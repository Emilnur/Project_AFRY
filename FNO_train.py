# -*- coding: utf-8 -*-
"""
Created on Sun Apr 27 22:30:59 2025

@author: Sankalp
"""

import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, random_split

from Data_reader import VesselFNO
from FNO import FNO

import matplotlib.pyplot as plt
import os

batch_size = 1

epochs = 200
learning_rate = 1e-3


device = torch.device('cuda'if torch.cuda.is_available() else 'cpu')

dataset = VesselFNO("Training")
train_size = int(0.8*len(dataset))
val_size = len(dataset)-train_size

training_set,validation_set = random_split(dataset,[train_size, val_size])

train_loader = DataLoader(training_set,batch_size=batch_size,shuffle=True)
val_loader = DataLoader(validation_set ,batch_size=batch_size,shuffle=False)

model = FNO().to(device)
criterion = nn.MSELoss()
optimizer = optim.Adam(model.parameters(),lr=learning_rate)
scheduler = optim.lr_scheduler.ReduceLROnPlateau(
    optimizer,
    mode='min',
    factor=0.5,       # reduce by 50%
    patience=10,      # wait 10 epochs with no improvement
    min_lr=1e-5
)

train_losses = []
val_losses = []

for epoch in range(epochs):
    model.train()
    total_loss = 0.0
    
    for inputs,targets in train_loader:
        inputs=inputs.to(device)
        targets=targets.to(device)
        
        optimizer.zero_grad()
        
        outputs = model(inputs)
        
        loss = criterion(outputs,targets)
        loss.backward()
        optimizer.step()
        total_loss+=loss.item()
    
    avg_loss = total_loss / len(train_loader)
    train_losses.append(avg_loss)
    
    model.eval()
    
    total_val_loss = 0.0
    
    with torch.no_grad():
        for inputs,targets in val_loader:
            inputs,targets = inputs.to(device),targets.to(device)
            outputs = model(inputs)
            loss = criterion(outputs,targets)
            
            total_val_loss += loss.item()
    
    avg_val_loss = total_val_loss/len(val_loader)
    val_losses.append(avg_val_loss)
    
    scheduler.step(avg_val_loss)
    current_lr = scheduler.optimizer.param_groups[0]['lr']
    
    print(f"Epoch [{epoch+1}/{epochs}] - Train Loss: {avg_loss:.6f} | Validation Loss: {avg_val_loss:.6f}")
    
torch.save(model.state_dict(), "model_FNO_final2.pth")

plt.figure()
plt.plot(train_losses, label="Train Loss")
plt.plot(val_losses, label="Validation Loss")
plt.xlabel("Epoch")
plt.ylabel("MSE Loss")
plt.title("Training and Validation Loss")
plt.legend()
plt.grid(True)
plt.savefig("loss_curve4.eps")
plt.show()


