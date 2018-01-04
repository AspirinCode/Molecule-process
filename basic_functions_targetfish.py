# -*- coding: utf-8 -*-

import numpy as np
from scipy.spatial import distance
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
import cPickle

# For data structuration.
def load_dict(struct_file):
    struct_dict = {}
    with open(struct_file) as f:
        for l in f:
            k, v  = eval(l.strip())
            struct_dict[k] = v
    return struct_dict

def load_bioinfo(bioinfo_file):
    go_dict, kegg_dict = {}, {}
    with open(bioinfo_file) as f:
        for l in f:
            tmp = eval(l.strip())
            go_dict[tmp[0]] = tmp[1]
            kegg_dict[tmp[0]] = tmp[2]
    return go_dict, kegg_dict

def load_by_group(dataset_file):
    buf = []
    with open(dataset_file) as f:
        for l in f:
            if l.startswith('Target') and len(buf)>1:
                yield buf
                buf = []
            buf.append(l)
        yield buf
def targetG_to_molEntry(target_group_file, dir_out):
    mol_fp_tars = {}
    for g in load_by_group(target_group_file):
        t_id = g[0].split('\t')[1]
        for l in g[1:]:
            mol_id, smi, b = l.strip().split('\t')
            if mol_fp_tars.has_key(mol_id):
                if b == '1': mol_fp_tars[mol_id][2].append(t_id)
                elif b == '0': mol_fp_tars[mol_id][3].append(t_id)
            else:
                entry = [mol_id, smi, [], []]
                if b == '1': entry[2].append(t_id)
                elif b == '0': entry[3].append(t_id)
                mol_fp_tars[mol_id] = entry
    with open(dir_out+'.txt', 'w') as f:
        for mol in mol_fp_tars:
            f.write(str(mol_fp_tars[mol]) + '\n')
    mol_smiles_tar_list = mol_fp_tars.values()
    cPickle.dump(mol_smiles_tar_list, open(dir_out+'.obj', 'w'))

#------------ Basic calculation ---------------------------------------
FP_SIZE = 2048
def get_ECFP4_bool(x, fp_size=FP_SIZE):
    mol = Chem.MolFromSmiles(x)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol,2, nBits=fp_size)
    return np.matrix([bool(int(i)) for i in fp.ToBitString()])
def get_ECFP4(x, fp_size=FP_SIZE):
    mol = Chem.MolFromSmiles(x)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol,2, nBits=fp_size)
    return fp
def get_topofp(x, fp_size=FP_SIZE, min_path=1, max_path=6):
    mol = Chem.MolFromSmiles(x)
    fp = Chem.RDKFingerprint(mol, minPath=min_path, maxPath=max_path, fpSize=fp_size)
    return fp
def get_MACCS(x):
    mol = Chem.MolFromSmiles(x)
    fp = MACCSkeys.GenMACCSKeys(mol)
    return fp
def cal_center(fp_matrix):
    center = np.mean(np.vectorize(int)(fp_matrix), axis=0)>=0.5
    return center
def cal_distance_bool(fps):
    return distance.jaccard(*fps)
def cal_distance_fp(fps):
    return 1 - DataStructs.FingerprintSimilarity(*fps)
def cal_c_score(fp_matrix_x, z=0.5):
    center = cal_center(fp_matrix_x)
    num = np.shape(fp_matrix_x)[0]
    tmp = zip(fp_matrix_x, [center for i in range(num)])
    distances_to_center = map(cal_distance_bool, tmp)
    return np.mean(distances_to_center) + z*np.std(distances_to_center)
def cal_query(fp_query, fp_matrix_x, center=None):
    num = np.shape(fp_matrix_x)[0]
    if center == None: center = cal_center(fp_matrix_x)
    query_vs_other = map(cal_distance_bool, zip(fp_matrix_x, [fp_query for i in range(num)]))
    return np.min(query_vs_other), np.mean(query_vs_other), np.std(query_vs_other)
def cal_query_sim(fp_query, fp_matrix_x, center=None):
    num = np.shape(fp_matrix_x)[0]
    if center == None: center = cal_center(fp_matrix_x)
    query_vs_other = map(cal_distance_bool, zip(fp_matrix_x, [fp_query for i in range(num)]))
    query_vs_other = [1-i for i in query_vs_other]
    return np.max(query_vs_other), np.mean(query_vs_other), np.std(query_vs_other)
def ClusterFps(fps, cutoff=0.8):
    from rdkit.ML.Cluster import Butina
    dists = []
    nfps = len(fps)
    for i in range(1, nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend([1-x for x in sims])
    cs = Butina.ClusterData(dists, nfps, cutoff, isDistData=True)
    return cs
def cal_BIO_distance(tarBIO1, tarBIO2):
    terms = set(tarBIO1.keys() + tarBIO2.keys())
    res = []
    for t in terms:
        if tarBIO1.has_key(t) and tarBIO2.has_key(t):
            res.append((tarBIO1[t]-tarBIO2[t])**2)
        elif tarBIO1.has_key(t) and (not tarBIO2.has_key(t)):
            res.append(tarBIO1[t]**2)
        elif (not tarBIO1.has_key(t)) and tarBIO2.has_key(t):
            res.append(tarBIO2[t]**2)
    return np.sqrt(np.sum(res))
def cal_BIO_distance_tc(tarBIO1, tarBIO2):
    terms_or = set(tarBIO1.keys()+tarBIO2.keys())
    terms_and = set(tarBIO1.keys()) & set(tarBIO2.keys())
    return 1-float(len(terms_and))/len(terms_or)
    
#-----------------------------------------------------------------------

def load_features(feature_data_file):
    X, y = [], []
    with open(feature_data_file) as f:
        for l in f:
            tmp = eval(l.strip())
            y.append(tmp[0])
            X.append(tmp[1:-1])
    X = np.matrix(X)
    y = np.matrix(y).T
    return X, y



def class_validation(y_prob, y, thr=0.5):
    tp, tn, fp, fn = 0, 0, 0, 0
    for a,b in zip(y_prob, y):
        if a>=thr and b==1: tp += 1
        elif a>=thr and b==0: fp += 1
        elif a<thr and b==1: fn += 1
        elif a<=thr and b==0: tn += 1
    precision = 0.0 if tp==0 and fp==0 else float(tp)/(tp+fp)
    recall = 0.0 if tp==0 and fn==0 else float(tp)/(tp+fn)
    accuracy = float(tn+tp)/(tn+fn+tp+fp)
    f1 = 0.0 if precision==0 and recall==0 else 2*precision*recall/(precision+recall)
    return precision, recall, accuracy, f1

# For single query molecule.
def rank_validation(top_pred, m_pred, act_tars):
    def cal_TP(n_pred, act_tars):
        return len([a for a in n_pred if a in act_tars])
    def cal_PR(TP_n, n):
        return float(TP_n)/n
    TOP = len(top_pred)
    M = len(act_tars)
    TP_n = cal_TP(top_pred, act_tars)
    PR_n = cal_PR(TP_n, TOP)
    RE_n = float(TP_n)/M
    if TP_n == 0.0: F_n = 0.0
    else: F_n = 2*PR_n*RE_n/float(PR_n+RE_n)
    PR1 = np.mean([cal_PR(cal_TP(m_pred[:i], act_tars),i) for i in range(1,M+1)])
    return PR_n, RE_n, F_n, PR1

def viz(result):
    result1 = []
    names = cPickle.load(open('names_dict.obj'))
    for r in result:
        l = list(r)
        l = [names[i] if str(i).startswith('CHEMBL') and i in names.keys() else i for i in l]
        result1.append(l)
    return result1
    












