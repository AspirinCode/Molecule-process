from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit import DataStructs
from rdkit.Chem import AllChem
import sascorer

def load_mols_all(molfile, form='smi'):
    if form in ['smi', 'txt']:
        return [Chem.MolFromSmiles(l.strip()) for l in open(molfile)]
    elif form == 'sdf':
        return [m for m in Chem.SDMolSupplier(molfile)]
    return None

def load_mol_single(molfile, form='smi'):
    if form in ['smi', 'txt']:
        with open(molfile) as f:
            for l in f:
                yield Chem.MolFromSmiles(l.strip())
    elif form == 'sdf':
        for m in Chem.SDMolSupplier(molfile):
            yield m


class MolProps:
    def __init__(self, molObj):
        self.molObj = molObj
        self.RO5_fail = 0
        self.HB_donor_num = Descriptors.NumHDonors(self.molObj)
        self.HB_acc_num = Descriptors.NumHAcceptors(self.molObj)
        self.logP = Descriptors.MolLogP(self.molObj)
        self.mol_weight = Descriptors.MolWt(self.molObj)
        self.SA_score = 0.0
        self.nearest_similarity = 0.0
        self.ecfp4 = None
    
    def get_SA_score(self):
        self.SA_score = sascorer.calculateScore(self.molObj)
    
    def get_RO5_fail(self):
        rules = [self.mol_weight < 600, self.HB_donor_num <= 5, \
                 self.HB_acc_num <= 10, 2.0 <= self.logP <= 4.0]
        self.RO5_fail = len(rules) - sum(rules)
    
    def get_and_return(self):
        self.get_SA_score()
        self.get_RO5_fail()
        return self.SA_score, self.RO5_fail
    
    def get_nearest_similarity(self, ref_fps):
        if self.ecfp4 == None:
            self.ecfp4 = AllChem.GetMorganFingerprint(self.molObj, 2)
        similarities = [DataStructs.DiceSimilarity(self.ecfp4, fp) for fp in ref_fps]
        self.ecfp4 = max(similarities)
        self.molObj.SetProp('Similarity', str(self.ecfp4))
        return self.ecfp4
        
    
    def set_mol_props(self):
        self.molObj.SetProp('logP', str(self.logP))
        self.molObj.SetProp('SA_score', str(self.SA_score))
        self.molObj.SetProp('HB_donor_number', str(self.HB_donor_num))
        self.molObj.SetProp('HB_acceptor_number', str(self.HB_acc_num))
        self.molObj.SetProp('MolWeight', str(self.mol_weight))
        

mol_file = 'D:/HAT/LSTM/lstm_fine_tune-180608/Character500000_Molecule500/S2/SP_p300_LSTM-180608-S2_lt-12.sdf'
mol_file_filtered = mol_file[:-4] + '_filtered'+ mol_file[-4:]
SIM_THR = 0.8
cal_similarity = True
ref_fps = None
if cal_similarity:
    ref_mols = [Chem.MolFromSmiles(l.strip()) for l in open('D:/HAT/LSTM/A485_analogues_active.smi')]
    ref_fps = [AllChem.GetMorganFingerprint(m, 2) for m in ref_mols]

if mol_file.endswith('.smi') or mol_file.endswith('.txt'):
    with open(mol_file_filtered, 'w') as f:
        for mol in load_mol_single(mol_file, 'smi'):
            mp = MolProps(mol)
            sa_score, ro5_fail = mp.get_and_return()
            if sa_score > 4.5 or ro5_fail > 1: continue
            if mp.get_nearest_similarity(ref_fps) >= SIM_THR: continue
            f.write(Chem.MolToSmiles(mol)+'\n')
elif mol_file.endswith('.sdf'):
    writer = Chem.SDWriter(mol_file_filtered)
    for mol in load_mol_single(mol_file, 'sdf'):
        mp = MolProps(mol)
        sa_score, ro5_fail = mp.get_and_return()
        if sa_score > 4.5 or ro5_fail > 1: continue
        if mp.get_nearest_similarity(ref_fps) >= SIM_THR: continue
        mp.set_mol_props()
        writer.write(mol)
    writer.close()

