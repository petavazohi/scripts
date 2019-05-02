#!/usr/bin/env python3

import os
import subprocess
import time
import psutil

while True:
    ls = os.listdir('.')
    for ifile in ls:
        if 'stop' in ifile:
            rf = open(ifile)
            data = rf.read()
            rf.close()
            key_words = [x for x in data.split('\n')]
            listProc = [] 
            for proc in psutil.process_iter(): 
                pInfo = proc.as_dict(attrs=['pid','name','username']) 
                listProc.append(pInfo) 
            for iproc in listProc:
                username = iproc['username']
                name     = iproc['name']
                pid      = iproc['pid']
                if name in key_words and username == 'petavazohi':
                    p = psutil.Process(pid=pid)
                    p.kill()
    time.sleep(180)
