import os
import pandas as pd

def summarize_vonmises_max(root_dir, output_file='sisummary.csv'):
    summary_data = []

    for subdir, _, files in os.walk(root_dir):
        for file in files:
            if file.lower().endswith('.csv'):
                file_path = os.path.join(subdir, file)
                try:
                    df = pd.read_csv(file_path, delimiter=';')
                    if 'VonMises' in df.columns:
                        max_vm = df['VonMises'].max()
                        summary_data.append({
                            'File': os.path.relpath(file_path, root_dir),
                            'Max VonMises': max_vm
                        })
                except Exception as e:
                    print(f"Error reading {file_path}: {e}")

    summary_df = pd.DataFrame(summary_data)
    summary_path = os.path.join(root_dir, output_file)
    summary_df.to_csv(summary_path, index=False)
    print(f"Summary saved to {summary_path}")

root_dir = r"C:\Users\librio\Documents\Clone_repo\Project_AFRY\doe_1"
summarize_vonmises_max(root_dir)