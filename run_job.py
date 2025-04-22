# run_job.py
import sys
from abaqus import mdb

inp_file = sys.argv[-2]
job_name = sys.argv[-1]

mdb.JobFromInputFile(name=job_name, inputFileName=inp_file)
mdb.jobs[job_name].submit()
mdb.jobs[job_name].waitForCompletion()
