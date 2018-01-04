#----------------------------------------------------------------------------
# Name: smarts_filter.py
# Author: Yolanda
# Instruction: Match some SMARTS and filter to a single sdf. Advanced hit expansion.
# mainList() input: sdf database, list of patterns to match (any), and output
#     file name (matched file).
#----------------------------------------------------------------------------

#-------------- Examples of SMARTS patterns ------------------------------
macrocyle = '[C,N,O,S,P;!r6;!r5;!r4;!r3;R]'
def cycleN(n): return '[C,N,O,S,P;r'+str(n)+']'
#--------------------------------------------------------------------------

from pybel import *
import os


# Judge by SMARTS:
def judge(smiles, smarts):
    mol = readstring('smi', smiles)
    sm = Smarts(smarts)
    if len(sm.findall(mol))>0: return True
    return False
def judgeList(smiles, smartsList): # Match any.
    for sm in smartsList:
        if judge(smiles, sm): return True
    return False



def mainList(sdfile, SMARTS_list, matchFile):
    with open(matchFile, 'w') as fout:
        ferr = open('errorCpd.sdf', 'w')
        import time
        t = time.time()
        for i,mol in enumerate(readfile('sdf', sdfile)):
            if i%1000==0: print i, (time.time()-t)/60
            try:
                smls = mol.write(format='smi')
                if judgeList(smls, SMARTS_list): # Match any.
                    fout.write(mol.write(format='sdf'))
                    print i, smls.strip()
            except Exception, e:
                print e
                ferr.write(mol.write(format='sdf'))
                continue
        ferr.close()
            


if __name__ == '__main__':
    # Cpd.29
    pattern_list = ['[O,N]c1c([CH2,NH,O][C,N,O])c2[cH,n][c,n][cH,n][cH,n]c2[cH,n][cH,n]1',\
                    'C1=[CH,N]c2c([O,NH][CH2,O,NH]1)[cH,n][cH,n]c3[cH,n][cH,n][c,n][cH,n]c32',\
                    '[N,C,O][CH2,NH,O]c1c([O,N])[cH,n][cH,n]c2c1[c,n,o,s][c,n,o,s][c,n,o,s]2',\
                    'C1=[CH,N]c2c([O,NH][CH2,O,NH]1)[cH,n][cH,n]c3[c,n,o,s][c,n,o,s][c,n,o,s]c32']
    
    dataDir = 'e:/Database/chempartner/'
    DBList = ['Specs_Nov2015_nopains.sdf']
    DBList = [i for i in os.listdir(dataDir) if i.endswith('.sdf')]
    for i,db in enumerate(DBList):
        print db
        mainList(dataDir+db, pattern_list, '29_chempartner'+str(i+1)+'.sdf')
    
