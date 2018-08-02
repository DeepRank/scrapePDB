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


class PDBselect(object):

    def __init__(self,seqsimgraph):

        # load the graph
        self.ssg = nx.read_gpickle(seqsimgraph)
        self.percent = self.ssg.percent

        # identify clusters
        self.clusters = list(nx.algorithms.connected_components(self.ssg))

   def get_unique_entries(self):

        # number of custers
        num_cluster = len(self.clusters)
        unique_pbds = []

        # go through the list
        for icluster in range(num_cluster):

            c = list(cluster[icluster])

            # generate the protein cluster graph
            g = get_protein_cluster_graph(c,percent=self.percent)

            # extract single pdbID per edge
            for edge in g.edges():
                pdb_list = g.edges[edge[0],edge[1]]['txt'].split('<br>')
                unique_pbds.append(self.select_edge_pdb(pdb_list))

        return unique_pbds


    def select_edge_pdb(self,pdblist):

        if len(pdblist) == 1:
            return pdblist[0]
        else:
            length = []
            for pdb in pdb_list:
                length.append(np.sum(self._get_pdb_length(pdb)))
            index = np.argmin(length)
            return pdblist[index]

    @staticmethod
    def _get_pdb_length(pdb):
        polymer = pypdb.get_all_info(pdb)['polymer']
        l = []
        for p in polymer:
            l.append(float(p['@length']))
        return l

    @staticmethod
    def get_protein_cluster_graph(cluster,percent=95):

        edges, nodes,dict_chains = {},{},{}
        Edge = namedtuple('Edge',['weight','txt'])
        Node = namedtuple('Node',['number','txt'])

        for pdb in tqdm(cluster):

            # get the info and chain labels
            info = pypdb.get_all_info(pdb)
            chains = [info['polymer'][0]['chain']['@id'],info['polymer'][1]['chain']['@id']]
            names = [None,None]

            # enumerate chans
            for ic,chain in enumerate(chains):

                # pdb.chain ID
                c = pdb+'.'+chain

                # add the pdb.chain ID to the general dict
                # {pdb.chain: prot name}
                if c not in dict_chains:

                    # use the macromolecule name
                    if 'macroMolecule' in info['polymer'][ic]:
                        if isinstance(info['polymer'][ic]['macroMolecule'],list):
                            names[ic] = info['polymer'][ic]['macroMolecule'][0]['@name']
                        else:
                            names[ic] = info['polymer'][ic]['macroMolecule']['@name']

                    # or the polymer description name
                    elif 'polymerDescription' in info['polymer'][ic]:
                        if isinstance(info['polymer'][ic]['polymerDescription'],list):
                            names[ic] = info['polymer'][ic]['polymerDescription'][0]['@description']
                        else:
                            names[ic] = info['polymer'][ic]['polymerDescription']['@description']

                    # add the pdb.chain to the dict
                    dict_chains[c] = names[ic]

                    # get the seq similarity of the chain
                    cluster,_ = get_seq_cluster_percent(c,percent=percent)
                    cluster = cluster['pdbChain']

                    # add all the chains with similar seq
                    # to the dict_chain {pdb.chain: prot_name}
                    if len(cluster)>0:
                        if not isinstance(cluster,list):
                            cluster = [cluster]
                        for n in cluster:
                            dict_chains[n['@name']] = names[ic]

                # reuse the previously defined entry
                else:
                    names[ic] = dict_chains[c]

                # add the node to the dict of Node namedtuples
                key = names[ic]
                if key not in nodes:
                    nodes[key] = Node(number=1,txt=pdb)
                else:
                    nodes[key] = nodes[key]._replace(number=nodes[key].number+1)
                    nodes[key] = nodes[key]._replace(txt=nodes[key].txt+'<br>'+pdb)

            # add the edge to the dict of Edge namedtuples
            key = tuple(names)
            if key not in edges:
                edges[key] = Edge(weight=1,txt=pdb)
            else:
                edges[key] = edges[key]._replace(weight=edges[key].weight+1)
                edges[key] = edges[key]._replace(txt= edges[key].txt + '<br>' + pdb)

        # Create the graph
        g = nx.Graph()

        for node_key, node_val in nodes.items():
            g.add_node(node_key,number=node_val.number,txt=node_val.txt)

        for edge_key, edge_val in edges.items():
            g.add_edge(edge_key[0],edge_key[1],weight=edge_val.weight,txt=edge_val.txt)

        return g