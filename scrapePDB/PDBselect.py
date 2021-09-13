#! /usr/bin/env python
import os
import pypdb
import multiprocessing
from tqdm import tqdm
import numpy as np
from time import time
import pickle
from functools import partial
import h5py

from pypdb.clients.search.search_client import perform_search, perform_search_with_graph
from pypdb.clients.search.search_client import ReturnType
from pypdb.clients.search.search_client import QueryGroup, LogicalOperator
from pypdb.clients.search.operators import text_operators


def _make_list(data):
    return [data] if not isinstance(data, list) else data


def screen_pdb(pdb, dict_cond=None):

    if dict_cond is None:
        dict_cond = {
            'method': 'xray', 'resolution': np.float('Inf'),
            'number_of_entity': 2,
            'type': ['protein'], 'len_min': 50, 'len_max': 5000
        }

    check = True
    # info = pypdb.get_entity_info(pdb)
    info = pypdb.get_all_info(pdb)['rcsb_entry_info']

    # method
    check *= info['experimental_method'] == dict_cond['method']
    if not check:
        reason = 'Incorrect Method : %s' % info['experimental_method']
        return check, reason


    if 'diffrn_resolution_high' not in info:
        info['diffrn_resolution_high'] = {'provenance_source':None, 'value':np.Inf}

    check *= float(info['diffrn_resolution_high']['value']) <= dict_cond['resolution']
    if not check:
        reason = 'Low Resolution (%1.2f)' % float(info['resolution'])
        return check, reason

    # number of entity
    # entity = _make_list(info['Entity'])
    # check *= len(entity) == dict_cond['number_of_entity']
    # if not check:
    #     reason = 'Wrong number of entitites %d' % len(entity)
    #     return check, reason

    # number/type of chain
    types = info['selected_polymer_entity_types']
    check *= types in dict_cond['type']

    if not check:
        reason = 'Incorrect chain Type %s' % types
        return check, reason

        # chain = _make_list(e['Chain'])
        # check *= len(chain) == bioAss

        # if not check:
        #     reason = 'Incorrect Number of Chain %d/%d' % (
        #         len(chain), bioAss)
        #     return check, reason

    # lentgth
    l = info['deposited_polymer_monomer_count']

    check *= (l >= dict_cond['len_min']
                and l <= dict_cond['len_max'])
    if not check:
        reason = 'Incorrect chain length %d' % l
        return check, reason

    return check, 'Entries accepted'


class PDBselect(object):

    def __init__(self, start=0, size=-1, nproc=1, tqdm=True,
                 method='X-ray', min_res=None, number_of_entity=2, types=['Protein (only)'],
                 len_min=50, len_max=None, outfile='dataset.hdf5', outfailedfile='problist.pkl'):

        self.start = start
        self.size = size
        self.nproc = nproc
        self.tqdm = tqdm
        self.outfile = outfile
        self.out_failed_file = outfailedfile

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

        if self.dict_cond['type'] == ['protein']:

            search_operator =  text_operators.ComparisonOperator(
                attribute="rcsb_entry_info.assembly_count",
                value=2,
                comparison_type=text_operators.ComparisonType.EQUAL
            )

            return_type = ReturnType.ENTRY
            self.allpdbs = perform_search(search_operator, return_type)

        else:
            self.allpdbs = self.get_all_pdb(self.start, self.size)

        print('PDBselect -> %d pdb retreived from the database' %
              len(self.allpdbs))

        if 1:
            if self.nproc == 1:
                self.results['ids'] = self.process_pdb(self.allpdbs)

            else:
                print('PDBselect -> Processing the files on %d procs' %
                    self.nproc)
                pool = multiprocessing.Pool(self.nproc)
                part_process = partial(
                    self._select_pdb, dict_cond=self.dict_cond)
                if self.tqdm:
                    self.results['ids'] = list(
                        tqdm(pool.imap(part_process, self.allpdbs), total=len(self.allpdbs)))
                else:
                    self.results['ids'] = pool.map(part_process, self.allpdbs)

            # remove the Nones
            self.results['ids'] = list(filter(None, self.results['ids']))

            # save results
            print('PDBselect -> %04d complexes found' %
                (len(self.results['ids'])))
            print('PDBselect -> Saving results in :', self.outfile)
            self.save_hdf5()

    def save_pkl(self):

        if not os.path.isfile(self.outfile):
            with open(self.outfile, 'wb') as f:
                pickle.dump(self.results, f)

        else:
            with open(self.outfile, 'rb') as f:
                old_res = pickle.load(f)

            if self.compatible_out(self.results, old_res):

                self.results['ids'] = old_res['ids'] + \
                    self.results['ids']
                self.results['start'] = min(
                    self.results['start'], old_res['start'])
                self.results['size'] = self.results['size'] + \
                    old_res['size']

                with open(self.outfile, 'wb') as f:
                    pickle.dump(self.results, f)

            else:

                print(
                    'PDBselect -> New search not compatible with existing data in %s' % self.outfile)
                a, b = self.outfile.split('.')
                newname = a+'_new.'+b
                print('PDBselect -> Saving results in :', newname)
                with open(newname, 'wb') as f:
                    pickle.dump(self.results, f)

    def save_hdf5(self):

        if os.path.isfile(self.outfile):
            print('PDBselect -> Outifle %s already exists' %
                  self.outfile)
            a, b = self.outfile.split('.')
            self.outfile = a+'_new.'+b
            print('PDBselect -> Saving results in :', self.outfile)

        f5 = h5py.File(self.outfile, 'w')
        grp = f5.create_group('PDBselect')

        for k in ['method']:
            grp.create_dataset(k, data=bytes(
                self.results[k], 'utf-8'))
        for k in ['resolution', 'number_of_entity', 'len_min', 'len_max', 'start', 'size']:
            grp.create_dataset(k, data=self.results[k])
        for k in ['type', 'ids']:
            grp.create_dataset(
                k, data=[i.encode('utf-8') for i in self.results[k]])
        f5.close()

    @staticmethod
    def load_hdf5(outfile):

        if not os.path.isfile(outfile):
            raise FileNotFoundError(outfile)

        f5 = h5py.File(outfile, 'r')
        if 'PDBselect' in f5:
            grp = f5['PDBselect']
            res = {}
            for k in ['method', 'type']:
                res[k] = grp[k].value.decode('utf-8')
            for k in ['resolution', 'number_of_entity', 'len_min', 'len_max', 'start', 'size']:
                res[k] = grp[k].value.decode('utf-8')
            res['ids'] = [v.decode('utf-8') for v in grp['ids'].value]
        else:
            res = None
        f5.close()
        return res

    @staticmethod
    def compatible_out(new, old):
        check = []
        keys = ['method', 'resolution', 'number_of_entity',
                'type', 'len_min', 'len_max']
        for k in keys:
            check.append(new[k] == old[k])
        return np.all(check)

    @staticmethod
    def get_all_pdb(start=0, size=-1):
        if size != -1:
            allpdbs = pypdb.get_all()[start:start+size]
        else:
            allpdbs = pypdb.get_all()[start:-1]
        return allpdbs

    def process_pdb(self, pdb_list):

        selected_pdbs, failed_pdbs = [], []

        if self.tqdm:
            pdb_list = tqdm(pdb_list)

        # screen all the pdbs
        for pdb in pdb_list:
            selected_pdbs.append(
                self._select_pdb(pdb, self.dict_cond))

        return selected_pdbs

    @staticmethod
    def _select_pdb(pdb, dict_cond):
        try:
            check, reason = screen_pdb(pdb, dict_cond)
            if check:
                return pdb
            # else:
            #    print(pdb, reason)

        except Exception as e:
            print('PDBselect -> Issue with ', pdb)
            print(e)
            # failed_pdbs.append(pdb)


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser('PDBdatabase scraper')

    parser.add_argument('--start', type=int, default=0,
                        help='start index for the search')
    parser.add_argument('--size', type=int, default=-1,
                        help='Number of pdbs to screen')
    parser.add_argument('--nproc', type=int, default=1,
                        help='Number of concurrent procs to use')

    parser.add_argument('--method', type=str,
                        default='xray', help='characterisation method')
    parser.add_argument('--min_res', type=float,
                        default=None, help='minimum resolution')
    parser.add_argument('--number_of_entity', type=int,
                        default=2, help='number of entities')
    parser.add_argument('--types', nargs='+', type=str,
                        default=['protein'], help='type of polymers')
    parser.add_argument('--len_min', type=int,
                        default=50, help='Minimum number of residues')
    parser.add_argument('--len_max', type=int,
                        default=None, help='Maximum number of residues')
    parser.add_argument(
        '--hdf5', type=str, default='dataset.hdf5', help='Name of the output file')

    parser.add_argument('--tqdm', type=bool, default=True,
                        help='use tqdm to monitor progress')

    args = parser.parse_args()

    pdbxt = PDBselect(start=args.start, size=args.size, nproc=args.nproc, tqdm=args.tqdm,
                      method=args.method, min_res=args.min_res,
                      number_of_entity=args.number_of_entity, types=args.types,
                      len_min=args.len_min, len_max=args.len_max,
                      outfile=args.hdf5)

    pdbxt.fetch()
