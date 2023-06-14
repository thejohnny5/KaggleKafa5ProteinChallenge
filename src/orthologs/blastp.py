from Bio.Blast.Applications import NcbimakeblastdbCommandline
import subprocess
import pandas as pd
from io import StringIO
import os
import time

TRAIN_FASTA = "./Data/cafa5protein/Train/train_sequences.fasta"
PATH_TO_NCBI_BIN = "./ncbi-blast-2.14/bin/"

def make_blastp_db(input_file: str) -> None:
    """Create blastdb in same directory as input file"""
    command = NcbimakeblastdbCommandline(cmd=PATH_TO_NCBI_BIN + "makeblastdb", dbtype="prot", input_file=input_file)
    subprocess.run(str(command), shell=True)


def query_db(db: str, fasta_query: str, tol: float = 1e-30, max_results: int = 100) -> pd.DataFrame:
    """Query the db for a fasta sequence."""
    headers = ["query acc.ver", "subject acc.ver", "Percent identity", "alignment length", "mismatches", "gap opens",
               "q.start", "q.end", "s.start", "s.end", "evalue", "bit score"]
    threads = os.cpu_count()
    command = f"{PATH_TO_NCBI_BIN}blastp -db {db} -query {fasta_query} -evalue {tol} -max_target_seqs {max_results} -outfmt 6 -num_threads {threads}"
    stdout = subprocess.run(command, shell=True, stdout=subprocess.PIPE)
    outstring = stdout.stdout.decode("utf-8")
    df = pd.read_csv(StringIO(outstring), sep="\t", header=None)
    df.columns = headers

    def filter_df(df: pd.DataFrame) -> pd.DataFrame:
        """Remove same values in query and subject"""
        mask = (df["query acc.ver"] == df["subject acc.ver"])
        return df[~mask]
    
    return filter_df(df)


class ChunkBlastP:
    """Break and run blastP in small chunks"""
    def __init__(self, fasta_file: str, chunck_size: int=100, chunk_directory: str="./tmp/", results_directory: str="./tmp_results/"):
        """Read fasta file. Break up into chunks of files."""
        self.fasta_file = fasta_file
        self.chunk_size = chunck_size
        self.fastas: list[str] = []
        self.chunk_directory = chunk_directory
        self.results_dir = results_directory
        self.results: list[str] = []
        
    
    def run(self):
        from_directory = self.chunk_directory
        #Try load from temp if no files - prepare fasta. Then call _run
        if (not os.path.isdir(from_directory)) or os.listdir(from_directory)==[]:
            print("Splitting...")
            self.prepare_fasta(tmp_directory=from_directory)
            print("Finished Splitting!")
        
        #Load fastas from split directory
        self.fastas = self.load_from_temp(from_directory)
        #Load results from results directory
        self.results = self.load_from_temp(self.results_dir)
       
        to_run = [file for file in self.fastas if file not in self.results] #Run those which don't exist in results directory
        to_run_full_path = [f"{from_directory}/{i}" for i in to_run]

        self._run_all(to_run_full_path)
        

            

    def load_from_temp(self, directory: str) -> list[str]:
        
        """Properties of temp files should be that they can be converted to int"""
        try:
            return [int(i) for i in os.listdir(directory)]
        except:
            Exception("Type error")
        #return [f"{tmp_directory}/{file}" for file in files]

    def prepare_fasta(self) -> None:
        tmp_directory = self.chunk_directory
        """Read fasta and break into chunks"""
        directory_save = tmp_directory
        chunk_num: int = 0
        current_file_num: int = 0
        with open(self.fasta_file, "r") as f:
            """Logic for splitting file"""
            while line := f.readline():
                if line.startswith(">"):
                    chunk_num += 1
                if not chunk_num <= self.chunk_size:
                    chunk_num=1
                    current_file_num+=1
                self._write_line(f"{directory_save}/{current_file_num}", line)
        
    
    def _write_line(self, file_name: str, line: str) -> None:
        """Write a single line to specified file"""
        with open(file_name, "a") as f:
            f.write(line)
            f.close()
    
    def _run(self, fasta_file: str) -> None:
        df = query_db(db=TRAIN_FASTA, fasta_query=fasta_file)
        file_name = os.path.split(fasta_file)[-1]
        df.to_csv(os.path.join(self.results_dir, file_name), index=False)

    def _runall(self, full_path_fastas: list[str]) -> None:
        """Run all fasta querries from a list of paths to each fasta file"""
        for i, fasta in enumerate(full_path_fastas):
            print(f"Running {i+1}/{len(full_path_fastas)}")
            now = time.time()
            self._run(fasta)
            post = time.time()
            print(f"Time elapsed: {post-now}")

            
if __name__ == "__main__":
    make_blastp_db(TRAIN_FASTA)
    #df = query_db(TRAIN_FASTA, TRAIN_FASTA)
    #df.to_csv("query_results.csv", header=False)
  
    par = ChunkBlastP(TRAIN_FASTA)
    par.run()