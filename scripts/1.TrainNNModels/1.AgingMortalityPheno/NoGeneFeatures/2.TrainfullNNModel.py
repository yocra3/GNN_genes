"""
Train model
docker run --gpus all -v $PWD:$PWD -w $PWD -it gnn_python:1.2 python
"""

import torch
from torch import nn
import pandas as pd
import sys
from torchinfo import summary

## Import custom modules
sys.path.append('./scripts/utils')
from dataloader_tsv import data_loader_tsv
from engine import train

sys.path.append('./scripts/1.TrainNNModels/1.AgingMortalityPheno')
from NNfullmodel import fullNN
from tensorboard_writer import create_writer

# Setup target device
device = "cuda" if torch.cuda.is_available() else "cpu"


BATCH_SIZE = 32

## Import data
train_data, test_data, n_genes = data_loader_tsv( "results/NNmodels/data/train.tsv", 
                                        "results/NNmodels/data/test.tsv",
                                        batch_size = 32)

gene_map =  pd.read_csv("results/NNmodels/data/gene_go_map.tsv", sep = "\t")
n_gos = len(gene_map["GO_ID"].unique())

## Define model
model0 = fullNN(n_genes, (n_gos, 10), 1).to(device)
model0 = torch.compile(model0)

## Check model

summary(model0, 
         verbose=0,
         input_size = (n_genes, ),
         col_names=["input_size", "output_size", "num_params", "trainable"],
         col_width=20,
         row_settings=["var_names"]
 )

## Define function loss and optimizer
loss_fn = nn.BCELoss()
optimizer = torch.optim.Adam(params=model0.parameters(), lr = 1e-3)

## Define accuracy
def accuracy(y_pred: torch.tensor,
             y: torch.tensor) -> int:
    """Compute accuracy for binary classification

    Args:
    y_pred: Model prediction
    y: Real values
    """   
    y_class = torch.round(y_pred)
    return (y_class == y).sum().item()

epochs = 10
# Start training with help from engine.py
train_res = train(model = model0,
             train_dataloader=train_data,
             test_dataloader=test_data,
             loss_fn=loss_fn,
             acc_fn = accuracy,
             optimizer=optimizer,
             epochs = epochs,
             writer=create_writer(experiment_name = "GO_genes_0.05",
                                       model_name = "FullNN",
                                       extra=f"{epochs}_epochs"),  
             device=device)

test_feats = []
test_labels = []
for feats, label in test_data:
    test_feats.append(feats)
    test_labels.append(label)

def make_predictions(model: torch.nn.Module, data: list, device: torch.device = device):
    pred_probs = []
    model.eval()
    with torch.inference_mode():
        for sample in data:
            # Prepare sample
            sample = sample.to(device) # Add an extra dimension and send sample to device
            # Forward pass (model outputs raw logit)
            pred_prob = model(sample)
            # Get pred_prob off GPU for further calculations
            pred_probs.append(pred_prob.cpu())           
    # Stack the pred_probs to turn list into a tensor
    return torch.cat(pred_probs)

pred_probs = make_predictions(model=model0, data=test_feats)


from sklearn.metrics import confusion_matrix

# Example true labels and predicted labels
predicted_labels = torch.round(pred_probs)

true_labels = torch.cat(test_labels)
# Compute confusion matrix
cm = confusion_matrix(true_labels, predicted_labels)
