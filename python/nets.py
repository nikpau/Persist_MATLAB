"""
In here, the neural networks for the lateral and longitudinal control policy are built

"""
# TODO What do the nets get as input, what do they output? 

import torch 
import torch.nn as nn
import numpy as np
import onnx
from onnx import numpy_helper


class LongitudinalNet(nn.Module):

    def __init__(self, state_dim, action_dim):
        super().__init__()

        # Definition of Network structure after agentData_lonConPol.mat
        self.net = nn.Sequential(
            nn.Linear(state_dim, 32),
            nn.ReLU(),
            nn.Linear(32,32),
            nn.ReLU(),
            nn.Linear(32,action_dim),
            nn.Tanh()
        )

    def forward(self,x):

        return self.net(x)


class LateralNet(nn.Module):

    def __init__(self, state_dim, action_dim):
        super().__init__()

        # Definition of Network structure after agentData_latConPol.mat
        self.net = nn.Sequential(
            nn.Linear(state_dim, 128),
            nn.ReLU(),
            nn.Linear(128,128),
            nn.ReLU(),
            nn.Linear(128,action_dim),
            nn.Tanh()
        )

    def forward(self,x):

        return self.net(x)

def _init_from_onnx(input_net, path_to_onnx):
    """
    This function will assign pretrained weights from an onnx
    file to a network generated with pytorch.
    The source and origin nets must be of same size,
    otherwise the function breaks
    """

    # Load model
    onnx_model = onnx.load(path_to_onnx)

    # Extract raw byte data for the weights
    raw_weights = onnx_model.graph.initializer

    # Init list for weights and biases
    wb = []

    i = 0
    while True:
        try:
            x = numpy_helper.to_array(raw_weights[i])
            x = np.squeeze(x)
            wb.append(x)
            i += 1
        except:
            break

    for idx, params in enumerate(input_net.parameters()):
        params.data = nn.parameter.Parameter(torch.tensor(wb[idx]))

    return input_net