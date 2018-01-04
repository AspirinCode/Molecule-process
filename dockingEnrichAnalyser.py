#------------------------------------------------------------------------------
# Name: dockingEnrichAnalyser.py
# Author: Yolanda
# Instruction: Calculate docking enrichment. For Glide docking results' sdf
#    file, this script will output Actives distribution of top % Results, also
#    output the AUC of ROC and a plot. User should input the sdf file and
#    implement 'telAct()' which can tell whether the molecule is active, if
#    necessary.
#------------------------------------------------------------------------------


# Information collection by reading input sdf file. Load molecules' titles and
#    docking scores.
def getInfo(sdfile):
    mol = []
    with open(sdfile) as f:
        ttl_on, score_on, ttl, score = True, False, '', 0.0
        for line in f:
            if ttl_on == True:
                ttl, ttl_on = line.strip(), False # get title
                continue
            if score_on == True:
                score, score_on = float(line.strip()), False # get docking score
                continue
            if line.strip() == '$$$$':
                ttl_on = True
                mol.append((score, ttl))
                continue
            if line.strip() == '> <r_i_docking_score>':
                score_on = True
    return sorted(mol)



# Tell the molecule is active or not. Should be RE-programmed for other circumstances.
def telAct(title):
    if title.startswith('ligand'): return False
    return True

def getActNum(mol):
    n = 0
    for score, ttl in mol:
        if telAct(ttl): n += 1
    return n



# Calculate True-positive and False-positive. Output the distribution and plot.
def percentTP(mol, p, actNum):
    rg, n, a = int(len(mol)*p), 0, 0
    for score, ttl in mol:
        if n == rg: break
        if telAct(ttl): a += 1
        n += 1
    return a, float(a)/actNum

def calRate(mol, actNum, dcyNum):
    tp, fp, rtNum, wrNum = 0.0, 0.0, 0, 0   
    x, y = [], []
    for score,ttl in mol:
        if telAct(ttl): rtNum += 1
        else: wrNum += 1
        tp, fp = float(rtNum)/actNum, float(wrNum)/dcyNum
        x.append(fp)
        y.append(tp)
    AUC = sum([(y[i]+y[i+1])*(x[i+1]-x[i])*0.5 for i in xrange(len(mol)-1)])
    return x, y, AUC

def pltROC(mol, actNum, dcyNum, filename):
    import matplotlib.pyplot as plt
    x, y, AUC = calRate(mol, actNum, dcyNum)
    print 'Area under ROC: '+ str(AUC)
    Fig = plt.figure()
    plt.title('ROC for '+filename)
    plt.xlabel('false positive')
    plt.ylabel('true positive')
    plt.plot(x, y)
    plt.plot([0,1],[0,1],'-.')    
    Fig.savefig('ROC_'+filename[:-4]+'.png')
    plt.show()
    return


import sys
filename = sys.argv[1]
#filename = 'enrich_????_dockingRes.sdf'
print 'Loading '+filename+' ...\n'
mol = getInfo(filename)
actNum = getActNum(mol)
dcyNum = len(mol)-actNum
print 'There are '+str(actNum)+' actives, and '+str(dcyNum)+' decoys.'
print
print 'Number of actives and true positive in top % results:'
print '1% : ', percentTP(mol, 0.01, actNum)
print '2% : ', percentTP(mol, 0.02, actNum)
print '5% : ', percentTP(mol, 0.05, actNum)
print '10% : ', percentTP(mol, 0.1, actNum)
print '15% : ', percentTP(mol, 0.15, actNum)
print '20% : ', percentTP(mol, 0.2, actNum)
print '50% : ', percentTP(mol, 0.5, actNum)
print
pltROC(mol, actNum, dcyNum, filename)

