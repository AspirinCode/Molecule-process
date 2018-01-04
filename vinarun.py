import os

ligdir = './ligands_pdbqt/'
outdir = './out_pdbqt/'
logdir = './vina_logs/'
config_file = 'config.txt'
if not os.path.exists(outdir): os.mkdir(outdir)
if not os.path.exists(logdir): os.mkdir(logdir)

for i,l in enumerate(os.listdir(ligdir)):
    if not l.endswith('.pdbqt'): continue
    os.system('vina --ligand %s --out %s --log %s --config %s'%\
              (ligdir+l, outdir+'out_'+l, logdir+'log_'+l[:-6], config_file))
    os.system('echo '+str(i))
    
