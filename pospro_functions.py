from abaqus import *
from abaqusConstants import *
from caeModules import *
from odbAccess import openOdb
import csv
import os
import glob

def find_odb(parent_dir):
    odb_info = []
    # Loop through all entries in the parent directory
    for folder_name in os.listdir(parent_dir):
        folder_path = os.path.join(parent_dir, folder_name)
        # Check if it's a folder
        if os.path.isdir(folder_path):
            # Look for .odb files inside the folder
            for file_name in os.listdir(folder_path):
                if file_name.endswith('.odb'):
                    odb_path = os.path.join(folder_path, file_name)
                    odb_folder = os.path.dirname(odb_path)
                    odb_info.append((odb_path, odb_folder))
    return odb_info

def PostProcess(doe_path, purge=0):
    odb_list = find_odb(doe_path)

    for odb_f, odb_dir in odb_list:
        try:
            odb = openOdb(odb_f)
        except Exception as e:
            print(f"error opening {odb_f}: {e}")
            continue        # Skip file if there's an error opening it

        step_name = list(odb.steps.keys())[0]
        step = odb.steps[step_name]
        last_frame=step.frames[-1]
        
        stress=last_frame.fieldOutputs['S'].getSubset(position=NODAL)
        assembly=odb.rootAssembly
        
        for instance_name in assembly.instances.keys():
            print(instance_name)
            
        instance = assembly.instances['PART-1-1']
        
        base_filename = os.path.splitext(os.path.basename(odb_f))[0]
        # Create a custom CSV filename using your naming scheme.
        csv_filename = os.path.join(odb_dir, f"{base_filename}_output.csv")

        with open(csv_filename, mode='w', newline='', encoding='utf-8') as file:
            writer = csv.writer(file, delimiter=';')
            writer.writerow(['Node label','x', 'y', 'z', 'S11', 'S22','S33','Mises'])
            
            coord_dict = {node.label: node.coordinates for node in instance.nodes}
        
            
            written_nodes = set()
            for value in stress.values:
                node_label = value.nodeLabel
                if node_label not in written_nodes:
                    node_coords = coord_dict[node_label]
                    print(node_coords[0])
                    writer.writerow([node_label, node_coords[0], node_coords[1], node_coords[2], value.data[0], value.data[1], value.data[2], value.mises])
                    written_nodes.add(node_label)
                
    odb.close()
    
    # if purge==1:
    #     os.remove(odb_path)