import numpy as np
from src.data.data_load import GOTerm
from typing import List
from src.data.data_load import GO_NAMESPACE_MAP

class InvalidNamespace(Exception):pass

"""In future restructure dictionaries. Use structure dict[str, tuple[class, str]] 
where class specifies either string or list"""
GO_KEYWORDS: dict[str, str] = {
    "id":"id",
    "name":"name",
    "namespace":"namespace",
    "def":"definition",
    "comment":"comment",
    "is_obsolete":"is_obsolete"
    }
GO_KEYWORDS_FOR_LISTS: dict[str, str] = {
    "synonym":"synonym",
    "consider":"consider",
    "is_a":"is_a",
    "alt_id":"alt_id",
    "subset":"subset",
    "xref":"xref",
    "relationship":"relationship",
    "replaced_by":"replaced_by"
    }

def load_numpy(path: str) -> np.array:
    """Load from .npy file"""
    return np.load(path)





class OBOParser:
    """Parses an .obo file"""
    def __init__(self):
        self.lines: List(str) = []

    def add_line(self, new_line: str) -> None:
        """Add line to this parser"""
        self.lines.append(new_line)
    
    def parse_lines(self) -> GOTerm:
        """Parse all go lines"""
        go_dictionary: dict[str, str] = {}
        for line in self.lines:
            key, value = self._parse_go_line(line)
            """Add key to go_dicionary"""
            go_dictionary = self.add_key(go_dictionary, key, value)
        if go_dictionary["namespace"] not in GO_NAMESPACE_MAP.keys():
            raise InvalidNamespace(f"The namespace '{go_dictionary['namespace']}' does not exist")
        go_dictionary["namespace"] = GO_NAMESPACE_MAP[go_dictionary["namespace"]]
        return GOTerm(**go_dictionary)
    
    def add_key(self, go_dictionary: dict[str, list[str] or str], key, value):
        """Add key value pair to dictionary. Checks if the value should be of string or list first"""
        if key in GO_KEYWORDS.keys():
            go_dictionary[GO_KEYWORDS[key]] = value
            return go_dictionary
        elif key in GO_KEYWORDS_FOR_LISTS.keys():
            return self.make_list(go_dictionary=go_dictionary, new_key=GO_KEYWORDS_FOR_LISTS[key], value=value)
        raise InvalidID(f"The key '{key}' is not found in either GO_KEYWORDS or GO_KEYWORDS_FOR_LIST")
    
    def make_list(self, go_dictionary: dict[str, list[str] or str], new_key: str, value: str) -> dict[str, list[str] or str]:
        """Either initialize an empty list or add an element to an existing list"""
        if new_key in go_dictionary.keys():
            go_dictionary[new_key].append(value)
            return go_dictionary
        go_dictionary[new_key] = [value]
        return go_dictionary


    def _parse_go_line(self, line: str) -> tuple[str, str]:
        """Parse a single line for GO"""
        def process_line(line: str) -> str:
            line = line.strip("\n")
            line = line.strip()
            return line
        split_line = line.split(":")
        assert len(split_line) > 1
        key = process_line(split_line[0])
        value = process_line(":".join(split_line[1:]))
        return (key, value)



class InvalidID(Exception): pass

def parse_GO(obo_file: str) -> List[GOTerm]:
    """Parse OBO file by checking for [Term] keyword followed by ID. Breaks at \n
    :param v: specifies verbose printing of missing keywords"""
    with open(obo_file, "r") as f:
        data: List[dict[str, str]] = []
        parser = OBOParser()
        start_parse = False
        while line := f.readline():
            if not line.strip() and start_parse:
                """Add dictionary to data and empty parser"""
                data.append(parser.parse_lines())
                start_parse = False
                parser = OBOParser()
                
            elif start_parse:
                """Make sure dictionary is nonempty before adding fields"""
                """Handle exception for field is not contained with finally to close file"""
                try:
                    parser.add_line(line)
                except InvalidID:
                    InvalidID(f"ID: id not valid key for goTerm")
                    f.close()
                
            if line.startswith("[Term]"):
                """Initialize a new Parser"""
                start_parse = True
            
            
            
        """EOF add dictionary to data and close file."""
        f.close()
    return data

if __name__ == "__main__":
    file = "./Data/cafa5protein/Train/go-basic.obo"
    data = parse_GO(obo_file=file)
    print(data[0])