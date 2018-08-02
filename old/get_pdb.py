import pypdb
import multiprocessing
from tqdm import tqdm
import numpy as np
from time import time


def chunks(full_list,nchunck):
	nitem = int(np.ceil(len(full_list)/nchunck))
	new_list = []
	for ic in range(nchunck):
		start, end = ic*nitem, (ic+1)*nitem
		new_list.append(full_list[start:end])
	return new_list

def check_structure(info):
	if 'polymer' in info:
		return len(info['polymer']) == 2
	else:
		return False

def check_polymer(info):
	check = True
	for p in info['polymer']:
		check *= _check_polymer(p)
	return check

def _check_polymer(p,min_length = 50):

	type_ = p['@type']
	if type_ != 'protein':
		return False

	num_chain = len(p['chain'])
	if num_chain != 1:
		return False

	length = int(p['@length'])
	if length < min_length:
		return False

	return True


def get_all_pdb(start=0,size=-1):
	if size != -1:
		allpdbs = pypdb.get_all()[start:start+size]
	else:
		allpdbs = pypdb.get_all()[start:-1]
	print('%d pdb retreived' %len(allpdbs))
	return allpdbs

def process_pdb(jobid,pdb_list):
	dimer = []
	pbl = []
	iadded = 0

	if jobid == 0:
		pdb_list = tqdm(pdb_list)

	for pdb in pdb_list:
		try:
			info = pypdb.get_all_info(pdb)
			if check_structure(info):
				if check_polymer(info):
					#print('[%04d] %s added to list' %(iadded,pdb))
					dimer.append(pdb)
					iadded += 1
		except :
			print('Issue with ', pdb)
			pbl.append(pdb)

	n = len(dimer)
	print('[proc%04d] %04d complexes found' %(jobid,n))
	f = open('pdblist%d.dat' %jobid,'w')
	for i in range(n):
		f.write('%s \n' %dimer[i])
	f.close()

	if len(pbl)>0:
		f = open('pbl%d.dat' %jobid,'w')
		for i in range(len(pbl)):
			f.write('%s \n' %pbl[i])
		f.close()

def dispatch_jobs(pdb_list,job_number):

	pdb_chuncks = chunks(pdb_list,job_number)
	jobs = []
	for i,pdb in enumerate(pdb_chuncks):
		j = multiprocessing.Process(target=process_pdb,args=(i,pdb))
		jobs.append(j)

	for j in jobs:
		j.start()

if __name__ == '__main__':
	pdbs = get_all_pdb(start=40000,size=10000)
	dispatch_jobs(pdbs,4)