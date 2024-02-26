"""
Pytorch template for creating a single NN model
"""
import torch
from torch import nn 

class fullNN(nn.Module):
    """Create a full NN network with two hiddedn layers.

    Args:
    input_shape: An integer indicating number of input features (e.g. genes).
    hidden_units: A tuple indicating number of hidden units in the two layers.
    output_shape: An integer indicating number of output units.
    """
    def __init__(self, input_shape: int, hidden_units: tuple, output_shape: int) -> None:
        super().__init__()
        self.NN1 = nn.Linear(in_features = input_shape,
                    out_features=hidden_units[0])
        self.NN2 = nn.Linear(in_features = hidden_units[0],
            out_features=hidden_units[1])
        self.out = nn.Linear(in_features = hidden_units[1],
            out_features=output_shape)
        self.relu = nn.ReLU()
        self.sigmoid = nn.Sigmoid()  
       
    def forward(self, x: torch.Tensor):
        return self.sigmoid(self.out(self.relu(self.NN2(self.relu(self.NN1(x))))))
