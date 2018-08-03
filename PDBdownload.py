#!/usr/bin/env python
import pypdb
import pickle
from tqdm import tqdm
import os

class PDBdownload(object):

    def __init__(self,pklfile,outdir='./'):

        data = pickle.load(open(pklfile,'rb'))
        self.pdblist = data['ids']
        self.outdir = outdir

    def download(self):

        for pdb in tqdm(self.pdblist):

            data = pypdb.get_pdb_file(pdb,filetype='pdb',compression=False)
            fname = os.path.join(self.outdir,pdb)
            fname += '.pdb'
            f = open(fname,'w')
            f.write(data)
            f.close()

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser('PDBdatabase scraper')
    parser.add_argument('pklfile',type = str, help='Pickle file containing the pdb ids')
    parser.add_argument('--outdir',type = str, help='Output directory for the dataset')
    args = parser.parse_args()

    pdb = PDBdownload(args.pklfile,outdir=args.outdir)
    pdb.download()