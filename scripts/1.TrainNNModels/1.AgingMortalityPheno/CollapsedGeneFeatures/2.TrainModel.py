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
from tensorboard_writer import create_writer

sys.path.append('./scripts/1.TrainNNModels/1.AgingMortalityPheno')
from NNfullmodel import fullNN

# Setup target device
device = "cuda" if torch.cuda.is_available() else "cpu"
## Remove warnings
import torch._dynamo
torch._dynamo.config.suppress_errors = True


BATCH_SIZE = 32

experiments = {}

percentages = ["0.5", "1", "2", "5", "10"]

## Load data
for percen in percentages:
    train_data, test_data, n_feats = data_loader_tsv( "results/NNmodels/data/collapse_feats_" + percen + "_percen_train.tsv", 
                                        "results/NNmodels/data/collapse_feats_" + percen + "_percen_test.tsv",
                                        batch_size = 32)
    experiments[percen] = [train_data, test_data, n_feats]
    ## Reduced experiments
    train_data_red, test_data_red, n_feats = data_loader_tsv( "results/NNmodels/data/collapse_feats_" + percen + "_percen_reduc_train.tsv", 
                                        "results/NNmodels/data/collapse_feats_" + percen + "_percen_reduc_test.tsv",
                                        batch_size = 32)
    new_lab = percen + "_red"
    experiments[new_lab] = [train_data_red, test_data_red, n_feats]



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
## Train models
for experiment, exp_list in experiments.items():
    print(f"Training {experiment} model...")
    train_data, test_data, n_feats = exp_list
    ## Define model
    model = fullNN(n_feats, (50, 10), 1).to(device)
    model = torch.compile(model)
    ## Define function loss and optimizer
    loss_fn = nn.BCELoss()
    optimizer = torch.optim.Adam(params=model.parameters(), lr = 1e-3)
    # Start training with help from engine.py
    train_res = train(model = model,
             train_dataloader=train_data,
             test_dataloader=test_data,
             loss_fn=loss_fn,
             acc_fn = accuracy,
             optimizer=optimizer,
             epochs = epochs,
             writer=create_writer(experiment_name = experiment,
                                model_name = "FullNN_CollapsedFeatures",
                                extra=f"{epochs}_epochs"),  
             device=device)

## Recover best model (full dataset 1%)
train_data, test_data, n_feats = experiments["1"]
## Define model
model0 = fullNN(n_feats, (50, 10), 1).to(device)
model0 = torch.compile(model0)


## Define function loss and optimizer
loss_fn = nn.BCELoss()
optimizer = torch.optim.Adam(params=model0.parameters(), lr = 1e-3)
# Start training with help from engine.py
train_res = train(model = model0,
             train_dataloader=train_data,
             test_dataloader=test_data,
             loss_fn=loss_fn,
             acc_fn = accuracy,
             optimizer=optimizer,
             epochs = epochs,
             writer = None,
             device=device)
test_pd =  pd.read_csv("results/NNmodels/data/collapse_feats_1_percen_test.tsv", sep = "\t")
test_feats =  torch.tensor(test_pd.iloc[:, :-1].values,dtype=torch.float32)
pred_probs = model0(test_feats.to(device))

test_pd['Pred_probs'] = pred_probs.cpu().detach().numpy()
test_pd.to_csv("results/NNmodels/collapse_feats_1_percen_test_prediction.tsv", sep = "\t",
               index=False)