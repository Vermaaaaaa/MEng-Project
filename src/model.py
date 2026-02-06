import torch
import torch.nn as nn

class Residual_Block(nn.Module):
    def __init__(self, n_inputs, n_outputs, kernel_size, stride, dilation, padding, dropout=0.2):
        super(Residual_Block, self).__init__()
        
        self.conv1 = nn.Conv1d(n_inputs, n_outputs, kernel_size, stride, padding, dilation)
        self.activation_function1 = nn.ReLU()
        self.dropout1 = nn.Dropout(dropout)

        self.conv2 = nn.Conv1d(n_inputs, n_outputs, kernel_size, stride, padding, dilation)
        self.activation_function2 = nn.ReLU()
        self.dropout2 = nn.Dropout(dropout)

        self.net = nn.Sequential(self.conv1, self.activation_function1, self.dropout1, 
                                 self.conv2, self.activation_function2, self.dropout2)

        self.residual_connection = nn.Conv1d(n_inputs, n_outputs, 1)

        self.init_weights()

    def init_weights(self):
        self.conv1.weight.data.normal(0, 0.01)
        self.conv2.weight.data.normal(0, 0.01)

        
        