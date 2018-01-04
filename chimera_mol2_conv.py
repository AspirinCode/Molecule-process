#import chimera
import os

def convert(dir_in, dir_out):
    chimera.runCommand('open '+dir_in)
    #chimera.runCommand('addh')
    chimera.runCommand('write format mol2 0 '+dir_out)
    chimera.runCommand('close session')


data_dir = 'd:/TIFP/set2/set2_converted/'
output_dir = 'd:/TIFP/set2/set2_chimera/'
if not os.path.exists(output_dir): os.mkdir(output_dir)
res_dirs = os.listdir(data_dir)
for d in res_dirs:
    os.mkdir(output_dir+d)
    print 'Converting', d, '...'
    for pdb_file in os.listdir(data_dir+d):
        convert(data_dir+d+'/'+pdb_file, output_dir+d+'/'+pdb_file[:-3]+'mol2')
