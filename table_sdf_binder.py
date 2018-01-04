import pybel

result_src = 'bromoscan_tmp.txt'
struc_src = 'W-cpd_bromodomain_SARlist.sdf'
combine_out = 'W-cpd_bromoscan_result.sdf'
scan_titles = ['FALZ_Inh@20uM', 'FALZ_Inh@2uM', 'BRD4(1)_Inh@20uM_CP', 'BRD4(1)_Inh@2uM_CP', \
              'BRD4(2)_Inh@20uM_CP', 'BRD4(2)_Inh@2uM_CP', 'CECR2_Inh@20uM', 'CECR2_Inh@2uM', \
              'TAF1(2)_Inh@20uM', 'TAF1(2)_Inh@2uM', 'TRIM24_Inh@20uM', 'TRIM24_Inh@2uM', \
              'BRD7_Inh@20uM', 'BRD7_Inh@2uM', 'SMARCA2_Inh@20uM', 'SMARCA2_Inh@2uM', \
              'BRD4(1)_Inh@200nM', 'BRD4(2)_Inh@200nM', 'CECR2_Inh@200nM', 'FALZ_Inh@200nM', \
              'TAF1(2)_Inh@200nM', 'TRIM24_Inh@200nM']

restable = [l.strip().split('\t') for l in open(result_src)]
restable = dict([(i[0], i[2:]) for i in restable])

with open(combine_out, 'w') as f:
    for mol in pybel.readfile('sdf', struc_src):
        name = mol.title
        values = restable[name]
        for t,v in zip(scan_titles, values):
            tmp = float(v)/100
            mol.data[t] = tmp
        f.write(mol.write('sdf'))


