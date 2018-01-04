#-----------------------------------------------------------------------------
# Name: sdf2sheet.py
# Author: Yolanda
# Instruction: Output structure image and some information of a sdf file to
#     an Excel sheet. Titles of information are defined in OUTLIST.
# Todo: Adjust width and height.
#-----------------------------------------------------------------------------


from pybel import *
from pyExcelerator import *
import Image
import os

#------------- User set -------------------------------------------#
OUTLIST = ['IDNUMBER','MDLNUMBER','MOLWEIGHT','LOGP','LOGS','TPSA','INCHI',\
           'INCHIKEY','SMILES']
SDFLIST = ['./kinaseout2.sdf']
excelName = './kinases-out2.xls'
#------------------------------------------------------------------#


def readSDF(sdf):
    for mol in readfile('sdf', sdf):
        mol.draw(show=False, filename='tmpImage.png')
        Image.open('tmpImage.png').save('tmpImage.bmp', format='bmp')
        os.remove('tmpImage.png')
        outlist = {}
        for i in OUTLIST:
            if mol.data.has_key(i): outlist[i] = mol.data[i]
            else: outlist[i] = ' '
        #outlist['title'] = mol.title
        yield outlist

def readSDFlist(sdflist):
    for sdf in sdflist:
        for molout in readSDF(sdf): yield molout
        

w = Workbook()
ws = w.add_sheet('Compound')
for i,v in enumerate(['Structure']+OUTLIST): ws.write(0, i, v)
n = 1
for mol in readSDFlist(SDFLIST):
    ws.insert_bitmap('tmpImage.bmp', n, 0, scale_x=0.5, scale_y=0.5)
    for i,v in enumerate(OUTLIST): ws.write(n, i+1, mol[v])
    n += 1
    if n%10==1: print 'Output %s entries ...'%(n-1)
w.save(excelName)
os.remove('tmpImage.bmp')
print 'Output %s entries all!'%(n-1)
