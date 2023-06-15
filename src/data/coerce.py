import pandas as pd
import numpy as np

def join_terms_and_embeds(terms_matrix: pd.DataFrame, embeds_matrix: pd.DataFrame, how='left') -> pd.DataFrame:
    """Indexes from both terms and embed must match
    Terms matrix is an [N, M] matrix where N is the protein term and M is the GO term. If the protein
    has a certain GO term, the value is 1, otherwise it is 0. Embeds Matrix is an [N,1024] matrix containing
    the protein term and the 1024 embeddings.
    
    Returns an [N, M+1024] Matrix """
    assert terms_matrix.shape[0] == embeds_matrix.shape[0]
    terms_matrix["key"] = terms_matrix.index.to_list()
    embeds_matrix["key"] = embeds_matrix.index.to_list()
    
    df = terms_matrix.merge(embeds_matrix, how="left", on="key")
    df.set_index("key", inplace=True)
    assert df.shape[0] == terms_matrix.shape[0]
    return df

def split_terms_and_embeds(joined_matrix: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Assumes that embed columns are named 0:1024
    Returns X, y of the joined matrix. For use after running train_test_split
    """
    cols = [i for i in range(0, 1024)]
    X = joined_matrix[cols]
    joined_matrix.drop(cols, inplace=True, axis=1)
    y = joined_matrix
    return X, y