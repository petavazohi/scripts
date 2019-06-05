#!/usr/bin/env python3

import os
import subprocess
import time


while True:
    path = os.getcwd()+os.sep
    ls = os.listdir(path)
    checked = False
    for item in ls :
        if '.sh' in item and '~' not in item:
            if os.path.exists(path+'WAITING'):
                os.remove(path+'WAITING')
            wf = open(path+'RUNNING','w')
            wf.close()
            try :
                subprocess.call(path+item,shell=True)
                os.rename(path+item,path+item.replace('.sh','.finished'))
                os.remove(path+'RUNNING')
            except:
                wf = open(path+"ERROR"+item.replace('.sh',''))
                wf.close()
                os.rename(path+item,path+item.replace('.sh','.error'))
        else : 
            checked = True
    if checked :
        wf = open(path+'WAITING','w')
        wf.close()
        time.sleep(30)

