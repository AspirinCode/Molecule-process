#-------------------------------------------------------------------------
# Author: Yolanda
# Instructions: Vina docking results pdbqt file should be convert to sdf.
#      But title (name) lost in pdbqt file, so we match to original mol2.
#      Vina score should be added to sdf when converting.
#-------------------------------------------------------------------------

from pybel import *
import os

oridir = './ori_ligands_mol2/'
resdir = './out_ligands_pdbqt/'
oriLst = [i for i in os.listdir(oridir) if i.endswith('mol2')]
resLst = [i for i in os.listdir(resdir) if i.endswith('pdbqt')]
outdir = './out_ligands_sdf/'
if not os.path.exists(outdir): os.mkdir(outdir)

# Title (name) lost in pdbqt file, so we match to original mol2 file.
def matchFile(fres):
    num = fres.split('_')[-1].split('.')[0]
    a = '_active' in fres
    for i in oriLst:
        num1 = i.split('_')[-1].split('.')[0]
        b = '_active' in i
        if num == num1 and a==b:
            return i


def getName(fori):
    with open(oridir+fori) as f:
        f.readline()
        name = f.readline().strip()
    return name


# Vina score should be added to sdf when converting.
def getVinaScore(fres):
    s = ''
    with open(resdir+fres) as f:
        f.readline()
        s = f.readline().strip()
    sl = s.split()
    if len(sl)>4:
        vinaScore = float(sl[3])
        return vinaScore
    return s


if __name__ == '__main__':
    for i,r in enumerate(resLst):
        o = matchFile(r)
        name = getName(o)
        vinaScore = getVinaScore(r)
        active = int('_active' in r)
        mol = readfile('pdbqt',resdir+r).next()
        pdb_s = mol.write('pdb') # First convert to pdb, then to sdf.
        mol = readstring('pdb', pdb_s)
        mol.title = name
        mol.data['vina_score'] = vinaScore
        mol.data['Active'] = active
        mol.write('sdf',outdir+r[:-5]+'sdf')
        print i, name
        







