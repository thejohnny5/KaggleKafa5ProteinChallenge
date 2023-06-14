import numpy as np
#from src.parsing import parse_GO
from src.data.data_load import Embeds

#####################################
#PARAMS
#####################################
path_to_ontology = "./Data/cafa5protein/Train/go-basic.obo"

#obo = parse_GO(path_to_ontology)
if __name__ == "__main__":
    embeds_root = "./Data/embeds/"
    train_embeds_path = embeds_root + "train_embeds.npy"
    train_labels_path = embeds_root + "train_ids.npy"
    test_embeds_path = embeds_root + "test_embeds.npy"
    test_labels_path = embeds_root + "test_ids.npy"

    train_set = Embeds(np.load(train_embeds_path), np.load(train_labels_path))
    test_set = Embeds(np.load(test_embeds_path), np.load(test_labels_path))
    print(len(train_set.IDS))
    print(len(test_set.IDS))
    print(len(set(train_set.IDS).intersection(set(test_set.IDS))))