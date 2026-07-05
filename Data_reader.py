# -*- coding: utf-8 -*-
"""
Created on Sun Apr 27 17:38:46 2025

@author: Sankalp
"""

import torch 
from torch.utils.data import Dataset
import pandas as pd
import numpy as np
import os
import re

class VesselFNO(Dataset):
    def __init__(self,result_dir):
        
        self.result_dir=result_dir
        self.jobs = sorted(os.listdir(result_dir))
        
    def __len__(self):
        return len(self.jobs)
    
    def thickness(self,inp_path):
        thickness_dict = {'Top_Flange':None,'Bottom_Flange':None,'Web':None}
        
        with open(inp_path,'r')as f:
            lines=f.readlines()
            
        for i,line in enumerate(lines):
            line = line.strip()
            
            if line.startswith('**Section') or line.startswith('** Section'):
                parts = line.split(':')
                if len(parts)>1:
                    current_section = parts[1].strip().upper()
                
            if line.startswith('*Shell Section'):
                if (i+1)<len(lines):
                    thickness_lines = lines[i+1].strip()
                    thickness_v = float(thickness_lines.split(',')[0])
                        
                    if current_section == 'TOPFLANGE':
                        thickness_dict['Top_Flange']=thickness_v
                            
                    elif current_section== 'BOTTOMFLANGE':
                        thickness_dict['Bottom_Flange']=thickness_v
                            
                    elif current_section=='WEB':
                        thickness_dict['Web']=thickness_v
                            
        return thickness_dict
               
    def __getitem__(self, idx):
        job_name = self.jobs[idx]
        
        if not job_name.startswith("run_"):
            raise ValueError(f"Unexpected job name:{job_name}")
            
        job_id = ''.join(filter(str.isdigit, job_name))
        job_folder=os.path.join(self.result_dir, job_name)
        csv = os.path.join(job_folder,f"pvessel_{job_id}.csv")
       # inp_path = os.path.join(job_folder,f"{job_id}.inp")
        
        # row = self.df.iloc[idx]
        # s_id,t1,t2,t3 = row['sample_id'],row['t1'],row['t2'],row['t3']
        
        # abaqus_data = os.path.join(self.nodes_dir,f'{s_id}.csv')
        # data = pd.read_csv(csv,delimiter=';')
        
        try:
            data = pd.read_csv(csv, delimiter=';')
        except FileNotFoundError:
            #print(f"[Warning] Missing file: {csv}, skipping sample.")
        # Option 1: Skip by picking another valid sample (recursive call)
           # new_idx = (idx + 1) % len(self)
           return self.__getitem__((idx + 1) % len(self))
        
        x=data['X'].values
        y=data['Y'].values
        z=data['Z'].values
        #node_id = data['Node label'].values
        stress = data['VonMises'].values
        thickness = data['Thickness'].values
        
        x = x.astype(np.float32)
        y = y.astype(np.float32)
        z = z.astype(np.float32)
        stress = stress.astype(np.float32)        
        thickness = thickness.astype(np.float32)
        n_nodes = len(x)
       
        #input_features = np.stack([x,y,z,np.full(n_nodes,thickness,dtype=np.float32)])
        input_features = thickness[np.newaxis,:]
        output_features = stress[np.newaxis,:]
        
        input_tensor = torch.tensor(input_features,dtype=torch.float32)
        output_tensor = torch.tensor(output_features,dtype=torch.float32)

        return input_tensor,output_tensor        
    