import torch
import torch.nn as nn

class Residual_Block(nn.Module):
    def __init__(self, n_inputs, n_outputs, kernel_size, stride, dilation, dropout=0.2):
        super(Residual_Block, self).__init__()
        
        padding = (kernel_size - 1) * dilation
        
        self.conv1 = nn.Conv1d(n_inputs, n_outputs, kernel_size, stride, padding=0, dilation=dilation)
        self.pad =  nn.ConstantPad1d((padding, 0), 0.0)
        self.activation_function1 = nn.ReLU()
        self.dropout1 = nn.Dropout(dropout)

        self.conv2 = nn.Conv1d(n_outputs, n_outputs, kernel_size, stride, padding=0, dilation=dilation)
        self.activation_function2 = nn.ReLU()
        self.dropout2 = nn.Dropout(dropout)

        self.net = nn.Sequential(self.pad, self.conv1, self.activation_function1, self.dropout1, 
                                 self.pad, self.conv2 ,self.activation_function2, self.dropout2)
        
        if n_inputs != n_outputs:
            self.residual_connection = nn.Conv1d(n_inputs, n_outputs, 1)
        else:
            self.residual_connection = None

        self.init_weights()

    def init_weights(self):
        self.conv1.weight.data.normal(0, 0.01)
        self.conv2.weight.data.normal(0, 0.01)
        if self.residual_connection is not None:
            self.residual_connection.weight.data.normal(0, 0.01)

    def forward(self, x):
        out = self.net(x)
        if self.residual_connection is None:
            residual = x
        else:
            residual = self.residual_connection(x)
        return nn.ReLU(out + residual)
        


        
        