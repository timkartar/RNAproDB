import os, sys
import subprocess

def getSS(prefix):
    subprocess.run(["x3dna-dssr","-i=../rnaprodb_frontend/public/cifs/{}-assembly1.cif".format(prefix),"--nested"])
    l = open("dssr-2ndstrs.dbn","r").readlines()
    seq = l[1].strip().split("&")
    sec = l[2].strip().split("&")
    return seq, sec

