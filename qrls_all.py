#!/usr/bin/env python3 

from pychemia.runner import get_jobs
import subprocess
import os

def qrls(job):
    job_id = job.split('.')[0]
    subprocess.call(f"qrls {job_id}", shell=True)

    
jobs = get_jobs(os.getenv("USER"))
for ijob in jobs:
    if jobs[ijob]['job_state'] == 'H':
        qrls(ijob)
