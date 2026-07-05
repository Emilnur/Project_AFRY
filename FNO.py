# -*- coding: utf-8 -*-
"""
Created on Sun Apr 27 18:38:42 2025

@author: Sankalp
"""

import torch
import torch.nn as nn
from neuralop.models import FNO1d

class FNO(nn.Module):
    def __init__(self):
        super(FNO,self).__init__()
        
        self.model=FNO1d(
            in_channels=1,
            out_channels=1,
            n_modes_height=16,
            hidden_channels=64,
            num_layers=4,
            lifting_channels=128,
            projection_channels=128
            )
        
    def forward(self,x):
        return self.model(x)