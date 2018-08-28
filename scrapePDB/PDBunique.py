#!/usr/bin/env python
import pypdb
import time
import multiprocessing
from tqdm import tqdm
import numpy as np
import networkx as nx
from collections import namedtuple
import pickle
import h5py

def print_id(item,id_,cond):
    if id_ == cond:
        print(item)


class PDBunique(object):

    def __init__(self,hdf5):

        # load the graph
        self.hdf5 = hdf5
        self.ssg, self.percent = self.load_seq_sim_graph()
        self.map = {}

        # identify clusters
        self.clusters = list(nx.algorithms.connected_components(self.ssg))


    def load_seq_sim_graph(self):

        f5 = h5py.File(self.hdf5,'r')
        grp = f5['PDBsim']
        nodes = grp['nodes'].value.astype('U')
        edges = [ tuple(e) for e in grp['edges'].value.astype('U') ]

        ssg = nx.Graph()
        ssg.add_nodes_from(nodes)
        ssg.add_edges_from(edges)
        percent = grp['percent'].value

        return ssg, percent

    def get_unique_entries(self,start,end):

        # number of custers
        num_cluster = len(self.clusters)
        if end != -1:
            num_cluster = end+1

        #self.unique_pbds = []
        #self.prot_graphs = []

        f5 = h5py.File(self.hdf5,'a')
        grp = f5.require_group('PDBunique')

        # go through the list
        for icluster in range(start,num_cluster):

            print('PDBUnique -> Cluster #%04d / %04d' %(icluster,num_cluster))
            c = list(self.clusters[icluster])

            subgrp_name = 'cluster_%04d' %icluster
            if subgrp_name not in grp:

                if len(c) > 1:

                    # generate the protein cluster graph
                    g = self.get_protein_cluster_graph(c,percent=self.percent)
                    #self.prot_graphs.append(g)

                    # extract single pdbID per edge
                    unique_pdbs = []
                    for edge in g.edges():
                        pdb_list = g.edges[edge[0],edge[1]]['txt'].split('<br>')
                        unique_pdbs.append(self.select_edge_pdb(pdb_list))

                else:
                    #self.prot_graphs.append(None)
                    unique_pdbs = c
                    g = None

                subgrp = grp.create_group(subgrp_name)
                self._save_cluster(subgrp,g,c,unique_pdbs)

        f5.close()

    @staticmethod
    def _save_cluster(subgrp,graph,cluster,unique):

        subgrp.create_dataset('unique',data=np.array(unique).astype('S'))
        subgrp.create_dataset('pdbids',data=np.array(list(cluster)).astype('S'))
        if isinstance(graph,nx.Graph):

            # nodes
            data = []
            for n,d in graph.nodes.data():
                data.append([n,d['number'],d['txt']])
            data = np.array(data).astype('S')
            subgrp.create_dataset('nodes',data = data)

            # nodes with data
            data = []
            for e1,e2,d in list(graph.edges.data()):
                data.append([e1,e2,d['weight'],d['txt']])
            data = np.array(data).astype('S')
            subgrp.create_dataset('edges',data = data)

    def get_unique_list(self,store=False):

        f5 = h5py.File(self.hdf5,'a')
        grp = f5['PDBunique']
        subgroup_list = filter(lambda x: x.startswith('cluster_'), list(grp.keys()))

        uniques = []
        for sub in subgroup_list:
            uniques += list(grp[sub+'/unique'].value.astype('U'))
        if store:
            grp.create_dataset('ids',data=np.array(uniques).astype('S'))
        f5.close()
        return uniques

    def save(self,outfile):
        results = {'graph':self.seqsimgraph,'percent':self.percent,'ids':self.unique_pbds,'map':self.map}
        f = open(outfile,'wb')
        pickle.dump(results,f)
        f.close()

    def select_edge_pdb(self,pdblist):

        if len(pdblist) == 1:
            self.map[pdblist[0]] = pdblist
            return pdblist[0]
        else:
            length = []
            for pdb in pdblist:
                length.append(np.sum(self._get_pdb_length(pdb)))
            index = np.argmin(length)
            self.map[pdblist[index]] = pdblist
            return pdblist[index]

    @staticmethod
    def _get_pdb_length(pdb):
        polymer = pypdb.get_all_info(pdb)['polymer']
        l = []
        for p in polymer:
            l.append(float(p['@length']))
        return l

    @staticmethod
    def get_protein_cluster_graph(cluster,percent):

        edges, nodes,dict_chains = {},{},{}
        Edge = namedtuple('Edge',['weight','txt'])
        Node = namedtuple('Node',['number','txt'])
        pdbid = None
        for pdb in tqdm(cluster):

            print_id(pdb,pdb,pdbid)

            # get the polymer infos
            check, niter = True, 0
            while check and niter < 10:
                try:
                    polymer = pypdb.get_all_info(pdb)['polymer']
                    check = False
                except:
                    print('PDBUnique -> Issue getting info for :', pdb)
                    print('PDBUnique -> Trying again in 5 sec')
                    time.sleep(5)
                    niter += 1
            if check:
                print('PDBUnique -> Entry %s ignored' %pdb)
                continue

            # get the chain labels
            chain_labels, chain_entity = [], []
            for ip,p in enumerate(polymer):

                chain = p['chain']

                # only conserve the first chain
                if isinstance(chain,list):
                    chain = chain[0]

                chain_labels.append(chain['@id'])
                chain_entity.append(ip)

            # init the names
            names = [None]*len(chain_labels)
            print_id(chain,pdb,pdbid)
            nup = 0
            # enumerate chans
            for ic,(chain,ip) in enumerate(zip(chain_labels,chain_entity)):

                # pdb.chain ID
                id_chain = pdb+'.'+chain

                # add the pdb.chain ID to the general dict
                # {pdb.chain: prot name}
                if id_chain not in dict_chains:

                    # use the macromolecule or polymer description name
                    for name_option,tag in zip(['polymerDescription','macroMolecule'],['@description','@name']):
                        if name_option in polymer[ip]:
                            if isinstance(polymer[ip][name_option],list):
                                names[ic] = polymer[ip][name_option][0][tag]
                            else:
                                names[ic] = polymer[ip][name_option][tag]
                            break

                        if names[ic] == 'Uncharacterized Protein':
                            names[ic] = 'UP_%03d' %nup
                            nup += 1

                    # add the pdb.chain to the dict
                    dict_chains[id_chain] = names[ic]

                    # get the seq similarity of the chain
                    check, niter = True, 0
                    while check and niter < 10:
                        try:
                            cluster,_ = pypdb.get_seq_cluster_percent(id_chain,percent=percent)
                            check = False
                        except:
                            print('PDBUnique -> Issue getting cluster for :', id_chain)
                            print('PDBUnique -> Trying again in 5 sec')
                            time.sleep(5)
                            niter += 1

                    if check:
                        print('PDBUnique -> Entry %s ignored' %id_chain)
                        cluster = []
                    else:
                        cluster = cluster['pdbChain']
                    print_id(cluster,pdb,pdbid)

                    # add all the chains with similar seq
                    # to the dict_chain {pdb.chain: prot_name}
                    if len(cluster)>0:
                        if not isinstance(cluster,list):
                            cluster = [cluster]
                        for n in cluster:
                            dict_chains[n['@name']] = names[ic]

                # reuse the previously defined entry
                else:
                    names[ic] = dict_chains[id_chain]

                # add the node to the dict of Node namedtuples
                key = names[ic]
                print_id(key,pdb,pdbid)
                if key not in nodes:
                    nodes[key] = Node(number=1,txt=id_chain)
                else:
                    nodes[key] = nodes[key]._replace(number=nodes[key].number+1)
                    if nodes[key].number < 35:
                        nodes[key] = nodes[key]._replace(txt=nodes[key].txt+'<br>'+id_chain)
                    elif nodes[key] == 35:
                        nodes[key] = nodes[key]._replace(txt=nodes[key].txt+'<br>'+'...')

            # add the edge to the dict of Edge namedtuples
            names.sort()
            key = tuple(names)
            print_id(key,pdb,pdbid)
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


    def map_find(self,pdbid):
        for k,v in self.map.items():
            if pdbid in v:
                return k,v

    def map_put(self,pdbid):

        key,vals = self.map_find(pdbid)
        index = self.unique_pdbs.index(key)
        print('Replacing %s with %s' %(key,pdbid))
        self.unique_pdbs[index] = pdbid


def load_prot_graph(f5name=None,index_cluster=None,grp=None):

    if f5name is None and grp is None:
        raise ValueError('f5name or grp must be specified')

    if f5name is not None:
        f5 = h5py.File(f5,'r')
        lgrp = f5['PDBunique/cluster_%03d' %index_cluster]
    else:
        lgrp = grp

    g = nx.Graph()
    nodes = grp['nodes'].value
    edges = grp['edges'].value
    #edges = [ tuple(e) for e in grp['edges'].value.astype('U') ]

    for n in nodes:
        node_id, number,txt = n
        node_id = node_id.decode('utf-8')
        g.add_node(node_id)
        g.nodes[node_id]['number'] = int(number.decode('utf-8'))
        g.nodes[node_id]['txt'] = txt.decode('utf-8')

    for e in edges:
        e1, e2, w, txt = e
        e1 = e1.decode('utf-8')
        e2 = e2.decode('utf-8')
        g.add_edge(e1,e2)
        g.edges[e1,e2]['weight'] = int(w.decode('utf-8'))
        g.edges[e1,e2]['txt'] = txt.decode('utf-8')

    return g


def plot_graph(G,fname,offline=False):


    if offline is False:
        import plotly.plotly as py
    else:
        import plotly.offline as py
    import plotly.graph_objs as go

    pos = nx.spring_layout(G)

    trace3_list = []
    middle_node_trace = go.Scatter(x=[],y=[],text=[],mode='markers',
                           hoverinfo='text',marker=go.Marker(opacity=0))
    for edge in G.edges():

        trace3 = go.Scatter(x=[],y=[],text=[],mode='lines',
                            line=go.Line(color='rgb(210,210,210)',width = min(G.edges[edge[0],edge[1]]['weight'],25)))

        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        trace3['x'] += [x0, x1, None]
        trace3['y'] += [y0, y1, None]
        trace3_list.append(trace3)

        middle_node_trace['x'].append((x0+x1)/2)
        middle_node_trace['y'].append((y0+y1)/2)
        middle_node_trace['text'].append(G.edges[edge[0],edge[1]]['txt'])

    node_trace = go.Scatter(
        x=[],
        y=[],
        text=[],
        mode='markers',
        hoverinfo='text',
        marker=dict(
            showscale=False,
            # colorscale options
            # 'Greys' | 'Greens' | 'Bluered' | 'Hot' | 'Picnic' | 'Portland' |
            # Jet' | 'RdBu' | 'Blackbody' | 'Earth' | 'Electric' | 'YIOrRd' | 'YIGnBu'
            colorscale='Electric',
            reversescale=True,
            color=[],
            size=[],
            colorbar=dict(
                thickness=15,
                title='Node Connections',
                xanchor='left',
                titleside='right'
            ),
            line=dict(color='rgb(50,50,50)',width=2)))

    for node in G.nodes():
        x, y = pos[node]
        node_trace['x'].append(x)
        node_trace['y'].append(y)
        node_trace['text'].append(''.join(node) + '<br>'+ G.nodes[node]['txt'])
        node_trace['marker']['color'].append(min(G.nodes[node]['number'],10))
        node_trace['marker']['size'].append(min(10+G.nodes[node]['number'],25))

    fig = go.Figure(data=[*trace3_list, middle_node_trace, node_trace],
                 layout=go.Layout(
                    title='<br>Network graph made with Python',
                    titlefont=dict(size=16),
                    showlegend=False,
                    hovermode='closest',
                    margin=dict(b=20,l=5,r=5,t=40),
                    annotations=[ dict(
                        text="Python code: <a href='https://plot.ly/ipython-notebooks/network-graphs/'> https://plot.ly/ipython-notebooks/network-graphs/</a>",
                        showarrow=False,
                        xref="paper", yref="paper",
                        x=0.005, y=-0.002 ) ],
                    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)))

    if not fname.endswith('.html'):
        fname += '.html'

    if offline is False:
        py.iplot(fig, filename=fname)
    else:
        py.plot(fig, filename=fname)



if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser('PDBdatabase scraper')
    parser.add_argument('hdf5',type = str, help='HDF5 file where the dataset is stored')
    parser.add_argument('-s','--start', type = int, default = 0, help='Index of the first clustr to analyze')
    parser.add_argument('-e','--end', type = int, default = -1, help='Index of the last cluster')
    parser.add_argument('--load', action="store_true", help='Only load the hdf5')

    args = parser.parse_args()

    pdb = PDBunique(args.hdf5)
    if not args.load:
        pdb.get_unique_entries(args.start,args.end)


