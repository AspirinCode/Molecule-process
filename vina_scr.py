import os
import time

def vinarun(l):
        print l
        t = time.time()
        os.system('vina --config conf.txt --ligand '+l+' --out res_'+l+' --log log_'+l)
        t = time.time() - t
        print 'It takes '+str(t)+'s to dock '+l

def scr():
        flst = [a for a in os.listdir('./') if a[:3]=='lig']
        map(vinarun, flst)


if __name__ == "__main__":
    scr()
