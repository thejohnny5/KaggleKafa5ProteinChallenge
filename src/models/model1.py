import tensorflow as tf
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.linear_model import Ridge

X_train = pd.DataFrame()
Y_train = pd.DataFrame()

model =Ridge(alpha=1.0)
model.fit(X_train, Y_train)