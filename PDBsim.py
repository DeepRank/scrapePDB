import pypdb
import multiprocessing
from tqdm import tqdm
import numpy as np
from time import time
import networkx as nx
import matplotlib.pyplot as plt
from plot_graph_plotly2 import plot_graph
from plot_graph_plotly import plot_graph_dimer
from collections import namedtuple
from tqdm import tqdm
import pickle

class SeqSimGraph(object):

    def __init__(self,xtfile,percent=95):

        self.xtfile = xtfile
        self.graphfile = graphfile
        self.percent = percent

        if self.graphfile is None:
            self.graphfile = 'SeqSimGraph_' + self.xtfile.split('.')[0] + '_' + str(percent) + '.pkl'

        # read the data
        self.pkl_xtdata = pickle.load(open(self.xtfile,'rb'))


    def get_graph(self, remove_self_loop = True):

        # get the pdbnames
        pdb_names = self.pkl_xtdata['ids']

        # init the graph
        self.seq_sim_graph = nx.Graph()
        self.seq_sim_graph.percent = self.percent
        self.seq_sim_graph.add_nodes_from(pdb_names)

        # go through all the pdbs
        for pdb in tqdm(pdb_names):

            # get info and chain labes
            info = pypdb.get_all_info(pdb)
            chain_labels = [info['polymer'][0]['chain']['@id'],info['polymer'][1]['chain']['@id']]

            # get all the neighbors for the all the chains
            neighbors = []
            for chain in chain_labels:

                check, niter = False, 0
                # get the cluster
                while not check and niter < 10:
                    try:
                        cluster,check = pypdb.get_seq_cluster_percent(pdb+'.'+chain, percent=percent)
                        cluster = cluster['pdbChain']
                    except :
                        print('Request failed for %s -> Trying again' %pdb)
                    niter += 1

                # add the pdb in the cluster to the
                # neighbor list
                if len(cluster)>0:
                    if not isinstance(cluster,list):
                        cluster = [cluster]
                    for n in cluster:
                        neighbors += [n['@name'].split('.')[0]]

            # get unique neighbors
            neighbors = list(set(neighbors))

            # add all the dges to the graph
            for n in neighbors:
                if n in pdb_names:
                    seq_sim_graph.add_edge(pdb,n)

    def save(self):
        nx.write_gpickle(self.seq_sim_graph,self.graphfile)


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser('Create the sequence similarity graph')

    parser.add_argument('--xtfile',type = str, default=None, help='PDBXtract file'
    parser.add_argument('--graphfile',type = str, default=None, help='sequence similarity graph file')
    parser.add_argument('--percent',type = int, default=95, help='Sequence similarity cutoff')
    parser.add_argument('--reomove_self_loop', type=bool, default=True, help='remove self edges'
    args = parser.parse_args()

    graph = SeqSimGraph(xtfile=args.xtfile,
                         graphfile=args.graphfile,
                         percent=args.percent)

    graph.get_graph(remove_self_loop = args.remove_self_loop)

    graph.save()

