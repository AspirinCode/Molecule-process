from rdkit import Chem
from rdkit.Chem import AllChem
import pybel
import os


def rmsd_inplace(refMol_dir, probMols_dirs):
    rmsd_result = []
    mols = Chem.MolFromPDBFile(refMol_dir)
    for i, probMol_dir in enumerate(probMols_dirs):
        probMol = Chem.MolFromPDBFile(probMol_dir)
        try:
            newconf = probMol.GetConformer(0)
        except Exception, e1:
            #print probMol_dir, e1
            continue
        mols.AddConformer(newconf, assignId=True)
        rmsd_result.append((probMol_dir, AllChem.GetConformerRMS(mols, 0, i+1, prealigned=True)))
    return rmsd_result

def batch_tag(pose_dir, ref_dir, tag_file, pose_format='pdb', ref_format='pdb'):
    refligs = os.listdir(ref_dir)
    f = open(tag_file, 'w')
    for code in os.listdir(pose_dir):
        probligs = [pose_dir+code+'/'+l for l in os.listdir(pose_dir+code) if 'ligand' in l]
        reflig = [ref_dir+r for r in refligs if code in r][0]
        try:
            rmsd_result = rmsd_inplace(reflig, probligs)
        except Exception, e:
            print 'WARNING:', code, e
            continue
        for lig,rmsd in rmsd_result:
            if rmsd <= 2.0:
                tag = 1
            elif rmsd > 4.0:
                tag = 0
            else: tag = None
            f.write('%s\t%s\t%s\t%s\n'%(tag, code, lig, rmsd))
        #print code, 'rmsd generated and taged.'
    f.close()

batch_tag('./set2/set2_converted/', \
          './set2/set2_ligand_native_pdb/', \
          './set2/tag_pose_rmsd_set2.txt')
            