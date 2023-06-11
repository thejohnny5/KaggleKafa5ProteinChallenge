from src.parsing import parse_GO


#####################################
#PARAMS
#####################################
path_to_ontology = "./Data/cafa5protein/Train/go-basic.obo"

obo = parse_GO(path_to_ontology)