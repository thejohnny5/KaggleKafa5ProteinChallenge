import numpy as np
from sklearn.metrics import accuracy_score, precision_score, recall_score, log_loss
from sklearn.metrics import f1_score, hamming_loss, roc_auc_score, cohen_kappa_score, matthews_corrcoef
from enum import Enum

class LossType(Enum):
    HAMMING="hamming"

class AverageType(Enum):
    MICRO="micro"
    MACRO="macro"


class ModelEvaluator:
    def __init__(self, metric: AverageType=AverageType.MICRO, loss: LossType=LossType.HAMMING):
        self.metric = metric
        self.loss = loss
        self.results: dict={}
    
    def evaluate(self, ground_truth, predictions) -> dict[str, str]:
        
        accuracy = accuracy_score(ground_truth, predictions)
        precision = precision_score(ground_truth, predictions, average=self.metric.value)
        recall = recall_score(ground_truth, predictions, average=self.metric.value)
        f1 = f1_score(ground_truth, predictions, average=self.metric.value)

        if self.loss == LossType.HAMMING:
            loss = hamming_loss(ground_truth, predictions)
        else:
            raise ValueError("Invalid loss type. Choose 'hamming'.")
        
        results = {
            'accuracy': accuracy,
            'precision': precision,
            'recall': recall,
            'f1': f1,
            'loss': loss
        }
        self.results = results
        
        return results
    def print_results(self) -> None:
        """Prints the result"""

        results = self.results
        try:
            print("Accuracy: {:.2f}%".format(results['accuracy'] * 100))
            print("Precision: {:.2f}%".format(results['precision'] * 100))
            print("Recall: {:.2f}%".format(results['recall'] * 100))
            print("F1 Score: {:.2f}%".format(results['f1'] * 100))
            print("Hamming Loss: {:.2f}%".format(results['loss'] * 100))
        except:
            Exception("Must run 'evaluate' first")


#evaluator = ModelEvaluator(metric=AverageType.MICRO, loss=LossType.HAMMING)
#results = evaluator.evaluate(ground_truth, predictions)
"""
print("Accuracy: {:.2f}%".format(results['accuracy'] * 100))
print("Precision: {:.2f}%".format(results['precision'] * 100))
print("Recall: {:.2f}%".format(results['recall'] * 100))
print("F1 Score: {:.2f}%".format(results['f1'] * 100))
print("Hamming Loss: {:.2f}%".format(results['loss'] * 100))"""