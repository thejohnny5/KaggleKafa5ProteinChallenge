import pandas as pd
import numpy as np



def convert_2d_to_1d_cols(df: pd.DataFrame) -> pd.DataFrame:
    """Take a 2D columns and convert to a 1D column by using the second value only"""
    cols = [j for _, j in df.columns]
    df.columns = cols
    return df

def prepare_labels(terms: pd.DataFrame[str, str], labels_to_include: list[str]) -> pd.DataFrame:
    """Prepare a multilabel dataset based on certain labels to include
    Expects that the terms dataframe contains protein name in column 0 and go term in column 1
    returns shape of [x, N] where N is the number of labels to include and x is the number of unique protein IDs"""

    def create_dumby_df(terms: pd.DataFrame, protein_col_name: str, go_col_name: str) -> pd.DataFrame:
        """Create a dummy dataframe so that all values are present"""
        search_terms = terms[protein_col_name].unique() #Since we will create a wide data frame from a short data frame need to make dummy data
        df_filler = pd.DataFrame(search_terms, columns=[protein_col_name])
        df_filler[go_col_name] = "blank"
        return df_filler
    
    protein_col = terms.columns[0]
    go_col = terms.columns[1]
    terms_to_include =  terms[terms[go_col].isin(labels_to_include)] #Eliminate labels that are not of interest

    df_combined = pd.concat([terms_to_include, create_dumby_df(terms=terms, protein_col_name=protein_col, go_col_name=go_col)])
    df_combined["value"] = 1

    #Create Pivot Table so that the GO term is now a column. Proteins with that GO term will have a value of 1.
    df_wide = pd.pivot(df_combined, index=["EntryID"], columns=["term"], values=["value"])

    #Fix Columns
    df_wide = convert_2d_to_1d_cols(df_wide)

    df_wide.drop("blank", axis=1, inplace=True)
    df_wide.fillna(value=0, inplace=True)

    return df_wide
    #cols = df_wide["term"].unique()
    #return df_wide[cols]

if __name__ == "__main__":
    #Load terms
    terms = pd.read_csv("./Data/cafa5protein/Train/train_terms.tsv", sep="\t")

    labels_to_include = ["GO:0015267", "GO:0015318"]
    df = prepare_labels(terms=terms, labels_to_include=labels_to_include)
    print(df["GO:0015267"].sum())