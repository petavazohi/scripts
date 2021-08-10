#!/usr/bin/env python3
import os
import pychemia
import shutil
import time
import sys
import datetime
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-np" ,
                    dest="np",
                    type=int,
                    action="store",
                    help="Number of MPI processes for the code",
                    default = '40')
parser.add_argument('-i',
                    "--input",
                    dest='input',
                    type=str,
                    help='A file with list of the structures you want to be calculated')
parser.add_argument("-xc",
                    dest="xc",
                    nargs='+',
                    type=str,
                    help='exchange correlation functional',
                    default = ['PBE'])
parser.add_argument("--max_jobs",
                    dest="max_jobs",
                    type=int,
                    help='maximum jobs to be submited simultaniously',
                    default=700)
parser.add_argument("-t",
                    "--template",
                    type=str,
                    help="template for the submitting job",
                    default=None)


args = parser.parse_args()   
nparal = args.np
max_jobs = args.max_jobs
finished = 0
if args.template is None and not os.path.exists('template.pbs'):
    wf = open('template.pbs','w')
    wf.write("""
    source ~/.bashrc
    module load atomistic/vasp/5.4.4_intel18_seq
    vasp_relax.py -np {np} -xc {xc} --tags ENCUT 550 --structure {poscar}
    """)
    wf.close()
    args.template='template.pbs'
    
rf = open(args.input)
structures_init = rf.readlines()
rf.close()

to_submit=[]
for ifile in structures_init:
    ifile = ifile.replace('\n','')
    original_ifile = ifile
    path = ifile+'_relax'
    if not os.path.exists(path):
        os.mkdir(path)
        
    for ixc in args.xc:
        ifile = original_ifile
        jobname = (ixc+'-'+ifile).replace('.vasp','').replace('.ascii','')
        path = ifile+'_relax'
        if not os.path.exists(path+os.sep+ixc):
            os.mkdir(path+os.sep+ixc)
            
        path = ifile+'_relax'+os.sep+ixc
        if os.path.exists(path+os.sep+"relax_report.json"):
            finished +=1
            continue
        else :
            contcar = ''
            poscar = ''
            if os.path.exists(path+os.sep+'CONTCAR'):
                rf = open(path+os.sep+'CONTCAR','r')
                contcar = rf.read()
                rf.close()
            if os.path.exists(path+os.sep+'POSCAR'):
                rf = open(path+os.sep+'POSCAR','r')
                poscar = rf.read()
                rf.close()
            if len(contcar) != 0:
                ifile = 'CONTCAR'
            elif len(poscar) !=0:
                ifile = 'POSCAR'
            
            else:
                ifile = original_ifile
                if not os.path.exists(path+os.sep+ifile):
                    shutil.copy(ifile,path+os.sep+ifile)
                
            job = pychemia.runner.PBSRunner(template='template.pbs',jobname=jobname,workdir=path)
            job.template_text = job.template_text.replace('{poscar}','"'+ifile+'"')
            job.template_text = job.template_text.replace('{np}',str(nparal))
            job.template_text = job.template_text.replace('{xc}',ixc)
            job.set_pbs_params(nodes=1,
                               ppn=nparal,
                               walltime=[4,0,0],
                               mail='petavazohi.hpc@gmail.com',
                               queue='standby')
            job.write()
            to_submit.append(job)
    
counter = 1
print("%d of the jobs already finished"%finished)
while True:
    jobs_submitted = pychemia.runner.get_jobs(user=os.getenv('USER'))
    njobs = len(jobs_submitted)
    
    if njobs < max_jobs and len(to_submit) != 0 :
        job = to_submit[0]
        current_time = datetime.datetime.now()
        print("{: >20} {} {: >30} {: >20} {: >20} {: >20}".format(counter,current_time,job.jobname,job.queue,job.nodes,job.ppn))
        try :
            to_submit[0].submit()
        except :
            
            time.sleep(60)
            continue
        counter +=1 
        del to_submit[0]
    elif len(to_submit) == 0:
        print("Submitted All, Exiting")
        exit(0)
    else :
        # print("{} jobs submitted, sleeping for 3 minutes".format(max_jobs))
        time.sleep(180)
