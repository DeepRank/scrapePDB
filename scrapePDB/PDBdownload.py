#!/usr/bin/env python
import pypdb
import pickle
import time
import multiprocessing
from functools import partial
from tqdm import tqdm
import os
import hdf5

class PDBdownload(object):

    def __init__(self,hdf5,outdir='./',nproc=1):

        f5 = h5py.File(hdf5)
        grp = f5['PDBunique']
        self.pdblist = [v.decode('utf-8') for v in grp['ids'].value]
        self.outdir = outdir
        self.nproc = nproc

    def download(self):

        if self.nproc == 1:
            for pdb in tqdm(self.pdblist):
                self._down(pdb,self.outdir)
        else:
            pool =  multiprocessing.Pool(self.nproc)
            part_process = partial(self._down,outdir=self.outdir)
            check = list(tqdm(pool.imap(part_process,self.pdblist), total=len(self.pdblist)))

    @staticmethod
    def _down(pdb,outdir):

            check,niter = True, 0
            while check and niter < 10:
                try:
                    data = pypdb.get_pdb_file(pdb,filetype='pdb',compression=False)
                    check = False
                except:
                    print('PDBDown -> Issue with :', pdb)
                    print('Trying again in 5 sec')
                    niter += 1
                    time.sleep(5)
            if check:
                print('PDBdown -> Could not download ', pdb)
                return False
            else:
                fname = os.path.join(outdir,pdb)
                fname += '.pdb'
                f = open(fname,'w')
                f.write(data)
                f.close()
                return True

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser('PDBdatabase scraper')
    parser.add_argument('hdf5',type = str, help='Pickle file containing the pdb ids')
    parser.add_argument('--outdir',type = str, default='./dataset/', help='Output directory for the dataset')
    parser.add_argument('--nproc',type = int, default=1, help='Number of procs')
    args = parser.parse_args()

    pdb = PDBdownload(args.hdf5,outdir=args.outdir,nproc=args.nproc)
    pdb.download()