import numpy as np
#from src.parsing import parse_GO

import os
import pandas as pd
from src.models.prepare_labels import prepare_labels
from src.data.data_load import load_embeds_to_df
from src.data.coerce import join_terms_and_embeds, split_terms_and_embeds
from sklearn.model_selection import train_test_split
from sklearn.linear_model import Ridge

#####################################
#PARAMS
#####################################
DATA_DIR = "./Competition_data/cafa5protein/Train/"
path_to_ontology = os.path.join(DATA_DIR, "go-basic.obo")
path_to_train_terms = os.path.join(DATA_DIR, "train_terms.tsv")
train_df = pd.read_csv(path_to_train_terms, sep="\t")
train_df = train_df[train_df.columns[[0, 1]]]
top_n_labels = train_df["term"].value_counts()[0:150].index.to_list()

label_matrix = prepare_labels(train_df, top_n_labels)
#print(label_matrix.head())

embeds_root = "./Competition_data/embeds/"

#################### Add function in data to handle
train_embeds_path = embeds_root + "train_embeds.npy"
train_labels_path = embeds_root + "train_ids.npy"
train_embeds_df = load_embeds_to_df(train_embeds_path, train_labels_path)
#print(train_embeds_df.head())
X_full = join_terms_and_embeds(label_matrix, train_embeds_df)
#############################################################

train, validate = train_test_split(X_full)
X_train, y_train = split_terms_and_embeds(train)
X_val, y_val = split_terms_and_embeds(validate)


print(X_train.shape)
print(y_train.shape)
print(X_val.shape)
print(y_val.shape)

##### To be handled by managed_models
#from sklearn.ensemble import RandomForestClassifier
#model = RandomForestClassifier(verbose=1)
import joblib
#model = Ridge(alpha=1)
#model.fit(X_train, y_train)
#joblib.dump(model, "test.pkl")
model = joblib.load("test.pkl")
import json
with open("test.json", "w") as outfile:
    json.dump(model.get_params(), outfile)
#print(time.time())
######################################

######### Need to be handled by managed_models and evaluator
y_proba = model.predict(X_val)
#from src.evaluation.eval_functions import ModelEvaluator
#evaluator = ModelEvaluator()
#evaluator.evaluate(y_val, y_proba)
#evaluator.print_results()
from sklearn.metrics import roc_auc_score
evaluation = []
"""print(y_val)
print(y_proba)"""
"""for i in range(y_proba.shape[1]):
    evaluation.append(roc_auc_score(y_val[:, i], y_proba[:, i]))
    print(evaluation)
    break"""
l = []
y_test=y_val.to_numpy()
for i in range(y_test.shape[1]):
    if len(np.unique(y_test[:, i]) ) > 1:
        s = roc_auc_score(y_test[:,i], y_proba[:,i]);
    else:
        s = 0.5
    l.append(s)        
    if i %10 == 0:
        print(i, s)

#evaluated = roc_auc_score(y_val, y_proba)
##################################################

#obo = parse_GO(path_to_ontology)
if __name__ == "__main__":
    
    test_embeds_path = embeds_root + "test_embeds.npy"
    test_labels_path = embeds_root + "test_ids.npy"

    #train_set = Embeds(np.load(train_embeds_path), np.load(train_labels_path))
    #test_set = Embeds(np.load(test_embeds_path), np.load(test_labels_path))
    #print(len(train_set.IDS))
    #print(len(test_set.IDS))
    #print(len(set(train_set.IDS).intersection(set(test_set.IDS))))