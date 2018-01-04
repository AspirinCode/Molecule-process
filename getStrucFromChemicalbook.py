#casList = [i.strip() for i in open('./otherCAS.txt')]
casList = ['8001-54-5','7361-61-7','27299-12-3','130693-82-2','13009-99-9',\
           '122111-03-9','121268-17-5','11041-12-6','1094-08-2','1944-12-3']
notFound = []
out = ''

import urllib
from pybel import *
import time

for c in casList:
    try:
        a = urllib.urlopen('http://www.chemicalbook.com/CAS/MOL/'\
                           + c + '.mol').read()
        mol = readstring('mol', a)
        mol.title = c
        out += mol.write('sdf')
        print c
        time.sleep(0.5)
    except Exception, e:
        print 'not found',c
        notFound.append(c)
        print e[:20]

with open('webNotFound1.txt', 'w') as f:
    f.write('\n'.join(notFound))

with open('webFound1.sdf', 'w') as f:
    f.write(out)
        
