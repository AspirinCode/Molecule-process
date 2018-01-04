import os
from pybel import *
from shutil import copy

def convert_pdb(pdbqt_file):
    mol = readfile('pdbqt', pdbqt_file).next()
    mol.write('pdb', pdbqt_file[:-2])
def convert_mol2(pdbqt_file):
    mol = readfile('pdbqt', pdbqt_file).next()
    mol.write('sdf', pdbqt_file[:-5]+'sdf')
    mol = readfile('sdf', pdbqt_file[:-5]+'sdf').next()
    mol.write('mol2', pdbqt_file[:-5]+'mol2')
    os.remove(pdbqt_file[:-5]+'sdf')
    
        

def batch_convert(data_dir, output_dir):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    res_dirs = os.listdir(data_dir)
    for d in res_dirs:
        if not os.path.exists(output_dir+d): os.mkdir(output_dir+d)
        print 'Converting', d, '...'
        for pdbqt_file in os.listdir(data_dir+d):
            convert_pdb(data_dir+d+'/'+pdbqt_file)
            copy(data_dir+d+'/'+pdbqt_file[:-2], output_dir+d+'/'+d+'_'+pdbqt_file[:-2])
            os.remove(data_dir+d+'/'+pdbqt_file[:-2])
            

batch_convert('./set2/vina_out_set2/', './set2/set2_converted/')
