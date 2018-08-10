#! /usr/bin/env python
import pypdb
import numpy as np
import networkx as nx
from tqdm import tqdm
import pickle
import multiprocessing
from functools import partial
import h5py

class PDBsim(object):

    def __init__(self,hdf5,percent=95,nproc=1,tqdm=True):

        self.hdf5 = hdf5
        self.percent = percent
        self.nproc = nproc
        self.tqdm = tqdm

        if hdf5 is not None:
            f5 = h5py.File(self.hdf5,'r')
            self.pdbnames = [ n.decode('utf-8') for n in f5['PDBselect/ids'].value ]


    def get_graph(self, remove_self_loop = True):

        # get the pdbnames
        pdb_names = self.pdbnames

        # init the graph
        self.seq_sim_graph = nx.Graph()
        self.seq_sim_graph.percent = self.percent
        self.seq_sim_graph.add_nodes_from(pdb_names)

        if self.nproc == 1:
            print('PDBsim -> Processing the graph')
            edges = self.get_all_edges(pdb_names)

        else:
            print('PDBsim -> Processing the graph on %d procs' %self.nproc)
            pool = multiprocessing.Pool(self.nproc)
            part_process = partial(self._get_edge_cluster,pdb_names=pdb_names,percent=self.percent)
            if self.tqdm:
                edges = list(tqdm(pool.imap(part_process,pdb_names), total=len(pdb_names)))
                edges = [e for edge_pdb in edges for e in edge_pdb]
            else:
                edges = pool.map(part_process,pdb_names)

        # get unique edges
        edges = list(set(edges))

        # add all the dges to the graph
        for e in edges:
            self.seq_sim_graph.add_edge(e[0],e[1])

        #remove self loops
        if remove_self_loop:
            self.seq_sim_graph.remove_edges_from([(p,p) for p in self.seq_sim_graph.nodes])

    def get_all_edges(self,pdb_names):

        edges = []
        print(len(pdb_names))
        # go through all the pdbs
        for pdb in tqdm(pdb_names):
            edges += self._get_edge_cluster(pdb,pdb_names,self.percent)
        return edges

    @staticmethod
    def _get_edge_cluster(pdb,pdb_names,percent):

        edges = []

        # get info and chain labes
        polymer = pypdb.get_all_info(pdb)['polymer']
        chain_labels = []
        for p in polymer:

            chain = p['chain']
            if not isinstance(chain,list):
                chain = [chain]
            for c in chain:
                chain_labels.append(c['@id'])

        # get all the neighbors for the all the chains
        for chain in chain_labels:

            check, niter = False, 0
            # get the cluster
            while not check and niter < 10:
                try:
                    cluster,check = pypdb.get_seq_cluster_percent(pdb+'.'+chain, percent=percent)
                    cluster = cluster['pdbChain']
                except Exception as e:
                    print(str(e))
                    print('Request failed for %s.%s -> Trying again' %(pdb,chain))
                niter += 1

            # add the (pdb,pdbneighbor) to the edge list
            if len(cluster)>0:

                if not isinstance(cluster,list):
                    cluster = [cluster]

                for n in cluster:
                    pdbid = n['@name'].split('.')[0]

                    # make sure the neighbor is in the pdb_names
                    if pdbid in pdb_names:
                        edges.append((pdb,pdbid))

        return edges

    def save(self):
        nx.write_gpickle(self.seq_sim_graph,self.graphfile)

    def save_hdf5(self):

        f5 = h5py.File(self.hdf5,'a')
        grp = f5.create_group('PDBsim')

        # nodes
        data = np.array(list(self.seq_sim_graph.nodes)).astype('S')
        grp.create_dataset('nodes',data = data)

        # edges
        data = np.array(list(self.seq_sim_graph.edges)).astype('S')
        grp.create_dataset('edges',data=data)

        # percent
        grp.create_dataset('percent',data=self.percent)


    def load_hdf5(self,f5=None,grp=None):

        if f5 is None and grp is None:
            raise ValueError('f5 or grp must be specified')

        if f5 is not None:
            f5 = h5py.File(h5,'r')
            lgrp = f5['PDBsim']
        else:
            lgrp = grp

        self.seq_sim_graph = nx.Graph()
        nodes = lgrp['nodes'].value.astype('U')
        self.seq_sim_graph.add_nodes_from(nodes)

        edges = [tuple(e) for e in lgrp['edges'].value.astype('U')]
        self.seq_sim_graph.add_edges_from(edges)

        if f5 is not None:
            f5.close()


    def plot_graph(self,fname,offline=True,noedge=False,ind_cluster=None):

        if offline is False:
            import plotly.plotly as py
        else:
            import plotly.offline as py
        import plotly.graph_objs as go


        cluster = list(nx.algorithms.connected_components(self.seq_sim_graph))

        if ind_cluster is not None:
            self.seq_sim_graph = self.seq_sim_graph.subgraph(cluster[ind_cluster])

        pos = nx.spring_layout(self.seq_sim_graph)

        edge_trace = go.Scattergl(
            x=[],
            y=[],
            line=dict(width=0.5,color='#888'),
            hoverinfo='none',
            mode='lines')

        if not noedge:
            for edge in self.seq_sim_graph.edges():
                x0, y0 = pos[edge[0]]
                x1, y1 = pos[edge[1]]
                edge_trace['x'] += [x0, x1, None]
                edge_trace['y'] += [y0, y1, None]

        node_trace = go.Scattergl(
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
                colorscale='Portland',
                reversescale=True,
                color=[],
                size=15,
                colorbar=dict(
                    thickness=15,
                    title='Nodes',
                    xanchor='left',
                    titleside='right'
                ),
                line=dict(width=2)))

        for node in self.seq_sim_graph.nodes():
            x, y = pos[node]
            node_trace['x'].append(x)
            node_trace['y'].append(y)

            for ic in range(len(cluster)):
                if node in cluster[ic]:
                    index = ic
                    break
            node_trace['text'].append(node + '[%d]' %index)
            node_trace['marker']['color'].append(index)

        if noedge:
            data = [node_trace]
        else:
            data = [node_trace, edge_trace]

        fig = go.Figure(data=data,
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

    parser = argparse.ArgumentParser('Create the sequence similarity graph')

    parser.add_argument('hdf5',type = str, help='HDF5 file containing the selected pdbs')
    parser.add_argument('--percent',type = int, default=40, help='Sequence similarity cutoff (default 40)')
    parser.add_argument('--remove_self_loop', type=bool, default=True, help='remove self edges')
    parser.add_argument('--nproc',type = int, default=1, help='Number of concurrent procs to use')
    parser.add_argument('--tqdm',type = bool, default=True,help='use tqdm to monitor progress')

    # parser.add_argument('--load',type = str, default=None, help='Load a pre-existing graph')
    # parser.add_argument('--offline',action='store_true', help='Plot offline')
    # parser.add_argument('--noedge',action='store_true', help='Do not plot the edge')
    # parser.add_argument('--cluster',type = int, default = None, help='plot a single cluster')

    args = parser.parse_args()

    graph = PDBsim(hdf5=args.hdf5,
                         percent=args.percent,
                         nproc=args.nproc,
                         tqdm=args.tqdm)

    graph.get_graph(remove_self_loop = args.remove_self_loop)
    graph.save_hdf5()

    # else:

    #     graph = SeqSimGraph(load=args.load)
    #     graph.plot_graph('simseq_4',offline=args.offline,noedge=args.noedge,ind_cluster=args.cluster)