#!/usr/bin/env python
#-------------------------------------------------------------------------------
# Name: Convert.py
# Purpose: use Prepare_ligand4.py and prepare_receptor4.py to convert pdb file to pdbqt file
# Author:  BearShare
# Created:     04-17-2013
# Copyright:   (c) BearShare
# Licence:     <your licence>
# Last Modified: 
#-------------------------------------------------------------------------------

import os

def GetList():
    w=[]
    f=open("list.txt")
    for line in f.readlines():
        w.append(line.strip().lstrip())
    return w

def main():
    file_list=GetList()
    for each in file_list:
        print each
        ilig_file="%s_lig.pdb"%each
        olig_file="%s_lig.pdbqt"%each
        ipoc_file="%s_poc.pdb"%each
        opoc_file="%s_poc.pdbqt"%each
        os.system("pythonsh prepare_ligand4.py -l %s -o %s"%(ilig_file,olig_file))
        os.system("pythonsh prepare_receptor4.py -r %s -o %s"%(ipoc_file,opoc_file))
if __name__=="__main__":
    main()
