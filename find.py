#!/usr/bin/env python

import subprocess
import argparse


def find(fname, address='.'):
    cmd = "find {} -name {}".format(address, fname)
    ret = subprocess.check_output(cmd, shell=True)
    ret = [x for x in  ret.decode('utf-8').split("\n") if len(x)!=0]
    return ret
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Find files and directories')
    parser.add_argument('name', type=str, nargs=1,
                                            help='name of the file or directory in search')
    parser.add_argument('--address', dest='address', default='.',
                                            help='address from where to start the search')

    args = parser.parse_args()
    found = find(args.name[0], args.address)
    for ifile in found:
        print(ifile)
        
