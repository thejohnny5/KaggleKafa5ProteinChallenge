import pandas as pd
from dataclasses import dataclass, field
from Bio.Seq import Seq
from typing import List, Optional
import numpy as np
from enum import Enum
from Bio import Entrez
import json

#PARAMS
ENTREZ_PARAMS = json.load(open("./src/data/entrez.params.json"))
Entrez.email = ENTREZ_PARAMS["email"]

class GONamespace(Enum):
    BIOLOGICAL_PROCESS="biological_process"
    MOLECULAR_FUNCTION="molecular_function"
    CELLULAR_COMPONENT="cellular_component"

GO_NAMESPACE_MAP: dict[str, GONamespace] = {
    "biological_process": GONamespace.BIOLOGICAL_PROCESS,
    "molecular_function": GONamespace.MOLECULAR_FUNCTION,
    "cellular_component": GONamespace.CELLULAR_COMPONENT
}

class EmptyVariablesException(Exception): pass

@dataclass
class Dataset:
    """Set up dataclass to be used for Train and Test data. Can also be used on validation data"""
    sequences: pd.DataFrame
    taxonomy: pd.DataFrame
    terms: Optional[pd.DataFrame]

@dataclass
class GOTerm:
    """Gets GO Terms to be used"""
    id: str
    name: str 
    namespace: GONamespace
    definition: Optional[str] = field(default_factory=str)
    comment: Optional[str] = field(default_factory=str)
    is_obsolete: Optional[str] = field(default_factory=str)
    synonym: Optional[List[str]] = field(default_factory=list)
    is_a: Optional[List[str]] = field(default_factory=list)
    consider: Optional[List[str]] = field(default_factory=list)
    alt_id: Optional[List[str]] = field(default_factory=list)
    subset: Optional[List[str]] = field(default_factory=list)
    xref: Optional[List[str]] = field(default_factory=list)
    relationship: Optional[List[str]] = field(default_factory=list)
    replaced_by: Optional[List[str]] = field(default_factory=list)

    def __str__(self) -> None:
        print(f"ID: {self.id}")
        print(f"Name: {self.name}")
        print(f"Namespace: {self.namespace}")

def load_embeds_to_df(path_to_embeds: str, path_to_embed_labels: str) -> pd.DataFrame:
    """Load embeds into a dataframe where the indexes are the labels"""
    embed = np.load(path_to_embeds)
    labels = np.load(path_to_embed_labels)
    embeds_df = pd.DataFrame(embed)
    embeds_df.index = labels
    return embeds_df

@dataclass
class Embeds:
    """Embeddings with ids for proteins"""
    EMBEDS: np.ndarray
    IDS: np.ndarray

def fetch_taxa(ids: List[int]) -> pd.DataFrame:
    """Get taxa information for a list of taxa ids"""
    ids = map(str, ids)
    handle = Entrez.efetch(id = ",".join(ids), db = "taxonomy", retmode = "xml")
    taxon_info = Entrez.read(handle)
    df = pd.DataFrame(taxon_info)
    df = df.astype({"TaxId": "int64"})
    return df[["TaxId", "ScientificName", "Rank", "Division", "GeneticCode", "MitoGeneticCode", "CreateDate"]]
