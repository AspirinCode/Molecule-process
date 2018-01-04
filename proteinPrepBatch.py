dir_in = './'
dir_exp = './preped/'

import os
import time
pdbList = [i for i in os.listdir(dir_in) if i.endswith('.pdb')]
flog = open('proPrepLog.txt', 'w')
for i,pr in enumerate(pdbList):
    os.system('echo %s_%s...'%(str(i), pr))
    os.system('/home/dddc/users/hao/schrodinger/utilities/prepwizard \
-keepfarwat -nometaltreat -disulfides -mse -fillsidechains -s %s %s_prep.pdb'\
%(pr, pr[:-4]))
    # Run prepwizard for one protein a time, or the file name may mismatch.
    # Each run may take several minutes.
    t = time.time()
    while True:
        if time.time()-t > 600:
            with open(pr[:-4]+'_prep.pdb', 'w') as f: f.write('error')
            flog.write('%s error!\n'%pr)
        if os.path.exists('./%s_prep.pdb'%pr[:-4]):
            os.system('mv %s_prep.pdb %s_prep.pdb'%(pr[:-4], dir_exp+pr[:-4]))
            break

os.system('echo done.')
