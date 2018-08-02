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

def get_seq_cluster_percent(pdb_id_chain,percent=95):
    """Get the sequence cluster of a PDB ID plus a pdb_id plus a chain,
    Parameters
    ----------
    pdb_id_chain : string
        A string denoting a 4 character PDB ID plus a one character chain
        offset with a dot: XXXX.X, as in 2F5N.A
    Returns
    -------
    out : dict
        A dictionary containing the sequence cluster associated with the PDB
        entry and chain
    Examples
    --------
    >>> sclust = get_seq_cluster('2F5N.A')
    >>> print(sclust['pdbChain'][:10])
    [{'@name': '4PD2.A', '@rank': '1'},
     {'@name': '3U6P.A', '@rank': '2'},
     {'@name': '4PCZ.A', '@rank': '3'},
     {'@name': '3GPU.A', '@rank': '4'},
     {'@name': '3JR5.A', '@rank': '5'},
     {'@name': '3SAU.A', '@rank': '6'},
     {'@name': '3GQ4.A', '@rank': '7'},
     {'@name': '1R2Z.A', '@rank': '8'},
     {'@name': '3U6E.A', '@rank': '9'},
     {'@name': '2XZF.A', '@rank': '10'}]
    """

    url_root = 'http://www.rcsb.org/pdb/rest/sequenceCluster?cluster=%s&structureId=' %str(percent)
    out = pypdb.get_info(pdb_id_chain, url_root = url_root)
    out = pypdb.to_dict(out)
    check = True
    return pypdb.remove_at_sign(out['sequenceCluster']), check


def get_graph(pdb_names,percent=95):

    graph = nx.Graph()
    graph.add_nodes_from(pdb_names)

    ip = 0
    for pdb in tqdm(pdb_names):

        info = pypdb.get_all_info(pdb)
        chain_labels = [info['polymer'][0]['chain']['@id'],info['polymer'][1]['chain']['@id']]

        neighbors = []
        for chain in chain_labels:

            check, niter = False, 0
            while not check and niter < 25:
                try:
                    cluster,check = get_seq_cluster_percent(pdb+'.'+chain,percent=percent)
                    cluster = cluster['pdbChain']
                except :
                    print('Request failed for %s -> Trying again' %pdb)
                niter += 1

            if len(cluster)>0:
                if not isinstance(cluster,list):
                    cluster = [cluster]
                for n in cluster:
                    neighbors += [n['@name'].split('.')[0]]

        neighbors = list(set(neighbors))

        for n in neighbors:

            if n in pdb_names:
                jp = pdb_names.index(n)
                graph.add_edge(pdb,n)
        ip += 1
    return graph

def analyze_cluster(cluster,percent=95):

    g = nx.Graph()
    dict_chains = {}

    nodes = {}
    edges,nodes = {},{}
    Edge = namedtuple('Edge',['weight','txt'])
    Node = namedtuple('Node',['number','txt'])
    for pdb in tqdm(cluster):

        info = pypdb.get_all_info(pdb)
        chains = [info['polymer'][0]['chain']['@id'],info['polymer'][1]['chain']['@id']]
        names = [None,None]

        for ic,chain in enumerate(chains):

            c = pdb+'.'+chain

            if c not in dict_chains:
                if 'macroMolecule' in info['polymer'][ic]:
                    if isinstance(info['polymer'][ic]['macroMolecule'],list):
                        names[ic] = info['polymer'][ic]['macroMolecule'][0]['@name']
                    else:
                        names[ic] = info['polymer'][ic]['macroMolecule']['@name']
                elif 'polymerDescription' in info['polymer'][ic]:
                    if isinstance(info['polymer'][ic]['polymerDescription'],list):
                        names[ic] = info['polymer'][ic]['polymerDescription'][0]['@description']
                    else:
                        names[ic] = info['polymer'][ic]['polymerDescription']['@description']

                dict_chains[c] = names[ic]

                cluster,_ = get_seq_cluster_percent(c,percent=percent)
                cluster = cluster['pdbChain']

                if len(cluster)>0:
                    if not isinstance(cluster,list):
                        cluster = [cluster]
                    for n in cluster:
                        dict_chains[n['@name']] = names[ic]

            else:
                names[ic] = dict_chains[c]

            key = names[ic]
            if key not in nodes:
                nodes[key] = Node(number=1,txt=pdb)
            else:
                nodes[key] = nodes[key]._replace(number=nodes[key].number+1)
                nodes[key] = nodes[key]._replace(txt=nodes[key].txt+'<br>'+pdb)

        key = tuple(names)
        if key not in edges:
            edges[key] = Edge(weight=1,txt=pdb)
        else:
            edges[key] = edges[key]._replace(weight=edges[key].weight+1)
            edges[key] = edges[key]._replace(txt= edges[key].txt + '<br>' + pdb)

    for node_key, node_val in nodes.items():
        g.add_node(node_key,number=node_val.number,txt=node_val.txt)

    for edge_key, edge_val in edges.items():
        g.add_edge(edge_key[0],edge_key[1],weight=edge_val.weight,txt=edge_val.txt)

    return g

def read_pdb(fname):
    f = open(fname,'r')
    data = f.readlines()
    f.close()

    for i in range(len(data)):
        data[i] = data[i].split()[0]
    return data

if __name__ == "__main__":

    # create the graph
    #fname = 'PDBId5'
    #pdb_names = read_pdb(fname)
    #g = get_graph(pdb_names,percent=40)
    #nx.write_gpickle(g,'graph5_40.pkl')

    # or read it
    g = nx.read_gpickle('graph5_40.pkl')
    g.remove_edges_from([(p,p) for p in g.nodes])


    # compute the cluster
    cluster = list(nx.algorithms.connected_components(g))
    num_cluster = len(cluster)
    unique = len(list(nx.isolates(g)))
    print('There are %d isolated conformations' %unique)
    print('There are %d clusters' %num_cluster)

    # icluster = 37
    # c0 = list(cluster[icluster])
    # g0 = analyze_cluster(c0,percent='40')
    # pos = nx.spring_layout(g0)
    # plot_graph(g0,pos,'cluster'+str(icluster))

    # analyze one cluster
    for icluster in range(num_cluster):
        c0 = list(cluster[icluster])
        if len(c0) > 1:
            g0 = analyze_cluster(c0,percent='40')
            unique += g0.number_of_edges()
            print('Cluster %03d -> %d/%d pdb possible' %(icluster,g0.number_of_edges(),len(c0)))
    print('There are %d isolated conformations' %unique)