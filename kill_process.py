#!/usr/bin/env python3

import psutil
import sys
import time
import datetime
import os


user = os.getenv('USER')
def check_proccess(name):
    """
    
    """
    for proc in psutil.process_iter():
        try:
            # Check if process name contains the given name string.
            if name.lower() in proc.name().lower() and proc.username() == user:
                now = datetime.datetime.now()
                print(now.strftime("%Y-%m-%d %H:%M:%S"), proc.pid, proc.name())
                proc.kill()
                return True
        except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
            pass
    return False


if len(sys.argv) < 2:
    print("Please input a list of process names")
    exit()
while True:
    for arg in sys.argv[1:]:
        check_proccess(arg)
    time.sleep(30)
