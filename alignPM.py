#----------------------------------------------------------------------------------
# Name: alignPM.py
# Author: Yolanda
# Instruction: Align PDBs. All the PDBs in dir_in will be aligned to the targetPDB,
#     by PyMol, and export the aligned PDBs to dir_exp.
#----------------------------------------------------------------------------------

from pymol import *
import os
from shutil import copy


def alignPymol(mobilePDB, targetPDB, dir_in, dir_exp):
    cmd.load(dir_in+mobilePDB)
    cmd.load(dir_in+targetPDB)
    cmd.align(mobilePDB[:-4], targetPDB[:-4])
    cmd.save(dir_exp+mobilePDB, selection=mobilePDB[:-4])
    cmd.delete('all')


def alignRun(dir_in, dir_exp, targetPDB):
    copy(dir_in+targetPDB, dir_exp+targetPDB)
    import __main__
    __main__.pymol_argv = ['pymol', '-qei'] # Run in quiet mode, but doesn't work
    finish_launching()
    pdbList = os.listdir(dir_in)
    pdbList.remove(targetPDB)
    for mobilePDB in pdbList: alignPymol(mobilePDB, targetPDB, dir_in, dir_exp)
    cmd.quit()
    

if __name__ == '__main__':
    dir_in = './original/'
    dir_exp = './aligned/'
    alignRun(dir_in, dir_exp, '3MXF.pdb')
