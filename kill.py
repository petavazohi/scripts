#!/usr/bin/env python3

import os
import subprocess
import time
import psutil


while True:
    wf = open('kill.out','w+')
    ls = os.listdir('.')
    wf.write(os.getcwd()+'\n')
    for ifile in ls:
        wf.write(ifile+'\n')
        if  ifile == 'stop':
            wf.write('in the loop\n')
            rf = open(ifile)
            data = rf.read()
            rf.close()
            key_words = [x for x in data.split('\n')]
            for x in key_words:
                wf.write('{}\n'.format(x))
            listProc = [] 
            wf.write('checking procees\n')
            for proc in psutil.process_iter(): 
                pInfo = proc.as_dict(attrs=['pid','name','username']) 
                listProc.append(pInfo) 
                iproc = pInfo
                username = iproc['username']
                name     = iproc['name']
                pid      = iproc['pid']
                wf.write('{} , {} , {} \n'.format(username,name,pid))
         
            for iproc in listProc:
                username = iproc['username']
                name     = iproc['name']
                pid      = iproc['pid']
                wf.write('{} , {} , {} \n'.format(username,name,pid))
         
                if name in key_words and username == os.getenv('USER'):
                    p = psutil.Process(pid=pid)
                    p.kill()
            os.remove(ifile)
    wf.close()
    time.sleep(30)
