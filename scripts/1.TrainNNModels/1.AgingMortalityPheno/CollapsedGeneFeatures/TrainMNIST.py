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

train_data, test_data, n_feats = data_loader_tsv( "mnist_train.tsv",  "mnist_test.tsv",
                                        batch_size = 32)


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
             writer=None, 
             device=device)

