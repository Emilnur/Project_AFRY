# -*- coding: utf-8 -*-
"""
Created on Tue Apr 29 10:33:41 2025

@author: Sankalp
"""

import torch
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from FNO import FNO


torch.serialization.safe_globals([torch._C._nn.gelu])
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

model = FNO()    
state = torch.load("model_FNO.pth",map_location='cpu',weights_only=False)
state.pop('_metadata',None)
model.load_state_dict(state)
model.to(device)
model.eval()



# sample_path = 'New_training/doe_1/run_0001/pvessel_0001.csv'
# # sample_path = 'Runfiles/run_0001/pvessel_0001.csv'

# node_data = pd.read_csv(sample_path,delimiter=';')

# x = node_data['X'].values.astype(np.float32)
# y = node_data['Y'].values.astype(np.float32)
# z = node_data['Z'].values.astype(np.float32)
# section = node_data['Section'].str.lower()
# # section = node_data['Set Name'].str.lower()
# t = np.zeros_like(x,dtype=np.float32)

# n_nodes = len(x)

# t1 = 11
# t2 =11
# t3 = 11

# t[section=='pad_section']=t1
# t[section=='vessel_section']=t2
# t[section=='nozzle_section']=t3

# # t[section=='surf_pad_section']=t1
# # t[section=='surf_vessel_section']=t2
# # t[section=='surf_nozzle_section']=t3


# input_features = t[np.newaxis,:]
# input_tensor = torch.tensor(input_features,dtype=torch.float32).unsqueeze(0).to(device)

# with torch.no_grad():
#     output = model(input_tensor)
    
# predicted_stress = output.squeeze().cpu().numpy()



def input_generator(t1,t2,t3):
    sample_path = 'New_training/doe_1/run_0001/pvessel_0001.csv'
    # sample_path = 'Runfiles/run_0001/pvessel_0001.csv'

    node_data = pd.read_csv(sample_path,delimiter=';')

    x = node_data['X'].values.astype(np.float32)
    y = node_data['Y'].values.astype(np.float32)
    z = node_data['Z'].values.astype(np.float32)
    section = node_data['Section'].str.lower()
    # section = node_data['Set Name'].str.lower()
    t = np.zeros_like(x,dtype=np.float32)

    n_nodes = len(x)
    
    t[section=='pad_section']=t1
    t[section=='vessel_section']=t2
    t[section=='nozzle_section']=t3
    input_features = t[np.newaxis,:]
    input_tensor = torch.tensor(input_features,dtype=torch.float32).unsqueeze(0).to(device)
    
    return input_tensor,x,y,z

def stress_eval(t1,t2,t3):
    input_tensor,x,y,z=input_generator(t1, t2, t3)
    input_tensor= input_tensor.to(device)
    with torch.no_grad():
        output = model(input_tensor)
        
    predicted_stress = output.squeeze().cpu().numpy()
    
    Max_vM = max(predicted_stress)
    
    return Max_vM, predicted_stress,x,y,z

t1 = 9.75
t2 =9
t3 = 10.5

Max_vM,predicted_stress,x,y,z = stress_eval(t1, t2, t3)

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

sc = ax.scatter(x, y, z, c=predicted_stress, cmap='turbo', s=5)

#plt.colorbar(sc, label='Predicted von Mises Stress (MPa)')
ax.set_xlabel('X (mm)')
ax.set_ylabel('Y (mm)')
ax.set_zlabel('Z (mm)')
#ax.set_title('Predicted Stress Distribution (3D View)')
plt.tight_layout()
plt.savefig("FNO1D_c1.png")

plt.show()
Max_vM = max(predicted_stress)

print(f'{Max_vM:.2f}')

    
