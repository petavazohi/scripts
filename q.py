#!/usr/bin/env python3 

from pychemia.runner import get_jobs
import os
import numpy as np


jobs = get_jobs(os.getenv("USER"))
job_queues = np.unique([jobs[x]['queue']for x in jobs])
job_queues = {x:{"count":0, 'Q':0, 'R':0} for x in job_queues}
for ijob in jobs:
    for iqueue in job_queues:
        if jobs[ijob]['queue'] == iqueue:
            job_queues[iqueue]['count'] += 1
            for status in ['Q', 'R']:
                if jobs[ijob]['job_state'] == status:
                    job_queues[iqueue][status] += 1
for iqueue in job_queues:
    print(
        "Total jobs in {: >10} : {: >5}, Running {: >5}, Queued {: >5}".format(iqueue,
                                                                               job_queues[iqueue]['count'],
                                                                               job_queues[iqueue]['R'],
                                                                               job_queues[iqueue]['Q']))
