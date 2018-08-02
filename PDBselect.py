#! /usr/bin/env python
import pypdb
import multiprocessing
from tqdm import tqdm
import numpy as np
from time import time
import pickle
from functools import partial

def _make_list(data):
    return [data] if not isinstance(data,list) else data

def screen_pdb(pdb,dict_cond=None):

    if dict_cond is None:
        dict_cond = {
        'method': 'xray', 'resolution': np.float('Inf'),
        'number_of_entity': 2,
        'type': ['protein'], 'len_min': 50, 'len_max': 5000
        }

    check = True
    info = pypdb.get_entity_info(pdb)

    # method
    check *= info['Method']['@name'] == dict_cond['method']
    if 'resolution' not in info:
        info['resolution'] = np.float('Inf')
    check *= float(info['resolution']) <= dict_cond['resolution']

    #number of entity
    entity = _make_list(info['Entity'])
    check *= len(entity) == dict_cond['number_of_entity']

    # number/type of chain
    bioAss = int(info['bioAssemblies'])
    types = _make_list(dict_cond['type'])
    for e in entity:
        check *= e['@type'] in types
        chain = _make_list(e['Chain'])
        check *= len(chain) == bioAss

    # lentgth
    if check == 1:
        polymer = pypdb.get_all_info(pdb)['polymer']
        for p in polymer:
            l = float(p['@length'])
            check *= (l > dict_cond['len_min'] and l < dict_cond['len_max'])

    return check

class PDBselect(object):

    def __init__(self, start = 0, size = -1, nproc = 1,tqdm=True,
                 method='xray', min_res = None, number_of_entity=2, types=['protein'],
                 len_min = 50, len_max = None, outfile='pdblist.pkl',outfailedfile='problist.pkl'):

        self.start = start
        self.size = size
        self.nproc = nproc
        self.tqdm = tqdm
        self.outfile = outfile
        self.out_failed_file = outfailedfile

        if self.nproc > 1:
            self._tmpfiles = ['_tmp_%03d' %i for i in range(self.nproc)]
            self._tmpfailedfiles = ['_tmp_failed_%03d' %i for i in range(self.nproc)]

        self.dict_cond = {
        'method': method, 'resolution': np.float('Inf') if min_res is None else min_res,
        'number_of_entity': number_of_entity, 'type': types,
        'len_min': 0 if len_min is None else len_min,
        'len_max': np.float('Inf') if len_max is None else len_max
        }

        self.results = self.dict_cond.copy()
        self.results['start'] = self.start
        self.results['size'] = self.size
        self.results['ids'] = []
        self.results['failed_ids'] = []

    def fetch(self):

        allpdbs = self.get_all_pdb(self.start,self.size)

        if self.nproc == 1:
            self.results['ids'] = self.process_pdb(allpdbs)

        else:
            print('PDBselect -> Processing the files on %d procs' %self.nproc)
            pool = multiprocessing.Pool(self.nproc)
            part_process = partial(self._select_pdb,dict_cond=self.dict_cond)
            if self.tqdm:
                self.results['ids'] = list(tqdm(pool.imap(part_process,allpdbs), total=len(allpdbs)))
            else:
                self.results['ids'] = pool.map(part_process,allpdbs)

        # remove the Nones
        self.results['ids'] = list(filter(None,self.results['ids']))

        # save results
        print('PDBselect -> %04d complexes found' %(len(self.results['ids'])))
        print('PDBselect -> Saving results in :', self.outfile)
        pickle.dump(self.results,open(self.outfile,'wb'))

    @staticmethod
    def get_all_pdb(start=0,size=-1):
            if size != -1:
                allpdbs = pypdb.get_all()[start:start+size]
            else:
                allpdbs = pypdb.get_all()[start:-1]
            print('PDBselect -> %d pdb retreived from the database' %len(allpdbs))
            return allpdbs

    def process_pdb(self,pdb_list):

        selected_pdbs, failed_pdbs  = [], []

        if self.tqdm:
            pdb_list = tqdm(pdb_list)

        # screen all the pdbs
        for pdb in pdb_list:
             selected_pdbs.append(self._select_pdb(pdb,self.dict_cond))

        return selected_pdbs


    @staticmethod
    def _select_pdb(pdb,dict_cond):
        try:
            if screen_pdb(pdb,dict_cond):
                    return pdb
        except :
            print('PDBselect -> Issue with ', pdb)
            #failed_pdbs.append(pdb)

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser('PDBdatabase scraper')

    parser.add_argument('--start',type = int, default=0, help='start index for the search')
    parser.add_argument('--size',type = int, default=-1, help='Number of pdbs to screen')
    parser.add_argument('--nproc',type = int, default=1, help='Number of concurrent procs to use')

    parser.add_argument('--method', type = str, default='xray',help='characterisation method')
    parser.add_argument('--min_res', type = float, default=None,help='minimum resolution')
    parser.add_argument('--number_of_entity', type = int, default=2, help='number of entities')
    parser.add_argument('--types', nargs = '+', type=str, default=['protein'], help='type of polymers')
    parser.add_argument('--len_min',type = int, default=50,help='Minimum number of residues')
    parser.add_argument('--len_max',type = int, default=None,help='Maximum number of residues')

    parser.add_argument('--tqdm',type = bool, default=True,help='use tqdm to monitor progress')

    args = parser.parse_args()

    pdbxt = PDBselect(start = args.start, size = args.size, nproc = args.nproc,tqdm=args.tqdm,
                 method=args.method, min_res=args.min_res,
                 number_of_entity=args.number_of_entity, types=args.types,
                 len_min = args.len_min, len_max = args.len_max)

    pdbxt.fetch()