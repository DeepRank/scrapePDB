import os, glob, shutil, re
from tqdm import tqdm

#from Bio.Blast.Applications import NcbipsiblastCommandline
#from pdb2sql.pdb2sqlcore import pdb2sql
#from .map_pssm2pdb import write_mapped_pssm_pdb

def get_bio_assembly(pdb):

    search_terms = 'APPLY THE FOLLOWING TO CHAINS'
    file = open(pdb,'r')
    chain_assemblies = []

    for line in file:

        if re.search(search_terms,line):
            l = line.split()
            chain_assemblies.append([l[-2].strip(','),l[-1].strip(',')])

        if line[0:4] == 'ATOM':
            break

    return chain_assemblies

def clean_pdb(pdb,chains,out='clean.pdb'):

    f = open(pdb,'r')
    fout = open(out,'w')
    for line in f:
        if line[0:4] == 'ATOM' and line[21] in chains:
            if line[21] == chains[0]:
                l = line[:21]+'A'+line[22:]
            elif line[21] == chains[1]:
                l = line[:21]+'B'+line[22:]
            fout.write(l)

def clean_dataset(input_dir,output_dir):

    pdb_list = list(filter(lambda x: x.endswith('.pdb'), os.listdir(input_dir)))

    for pdb in tqdm(pdb_list):

        inp = os.path.join(input_dir,pdb)
        outp = os.path.join(output_dir,pdb)

        chains = get_bio_assembly(inp)
        clean_pdb(inp,chains[0],outp)
