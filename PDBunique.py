#!/usr/bin/env python
import pypdb
import multiprocessing
from tqdm import tqdm
import numpy as np
import networkx as nx
from collections import namedtuple
import pickle


def print_id(item,id_,cond):
    if id_ == cond:
        print(item)


class PDBunique(object):

    def __init__(self,seqsimgraph,cluster=None,outfile='pdb_uniques.pkl'):

        # load the graph
        self.seqsimgraph = seqsimgraph
        self.ssg = nx.read_gpickle(self.seqsimgraph)
        self.percent = self.ssg.percent
        self.cluster_index = cluster
        self.outfile = outfile
        self.map = {}

        # identify clusters
        self.clusters = list(nx.algorithms.connected_components(self.ssg))

    def get_unique_entries(self):

        # number of custers
        num_cluster = len(self.clusters)
        self.unique_pbds = []

        # go through the list
        for icluster in range(num_cluster):

            print('DPBUnique -> Cluster #%03d / %03d' %(icluster,num_cluster))
            c = list(self.clusters[icluster])
            if len(c) > 1:

                # generate the protein cluster graph
                g = self.get_protein_cluster_graph(c,percent=self.percent)

                # extract single pdbID per edge
                for edge in g.edges():
                    pdb_list = g.edges[edge[0],edge[1]]['txt'].split('<br>')
                    self.unique_pbds.append(self.select_edge_pdb(pdb_list))

            else:
                self.unique_pbds.append(c[0])

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
    def get_protein_cluster_graph(cluster,percent=95):

        edges, nodes,dict_chains = {},{},{}
        Edge = namedtuple('Edge',['weight','txt'])
        Node = namedtuple('Node',['number','txt'])
        pdbid = None
        for pdb in tqdm(cluster):

            print_id(pdb,pdb,pdbid)

            # get the polymer ifos
            polymer = pypdb.get_all_info(pdb)['polymer']

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
                    cluster,_ = pypdb.get_seq_cluster_percent(id_chain,percent=percent)
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

def plot_graph(G,fname):


    import plotly.plotly as py
    import plotly.graph_objs as go

    pos = nx.spring_layout(G)

    trace3_list = []
    middle_node_trace = go.Scatter(x=[],y=[],text=[],mode='markers',
                           hoverinfo='text',marker=go.Marker(opacity=0))
    for edge in G.edges():

        trace3 = go.Scatter(x=[],y=[],text=[],mode='lines',
                            line=go.Line(color='rgb(210,210,210)',width = G.edges[edge[0],edge[1]]['weight']))

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
        node_trace['marker']['size'].append(min(10+G.nodes[node]['number'],50))

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

    py.iplot(fig, filename=fname)



if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser('PDBdatabase scraper')

    parser.add_argument('--ssg',type = str, help='Structure Similarity graph')
    parser.add_argument('--cluster', type = int, default = None, help='Index of the cluster to analyze')
    parser.add_argument('--outfile', type = str, default = 'pdb_unique.pkl', help='Name of the output file')
    parser.add_argument('--load', type = str,  default = None, help='Load a precomputed data file')
    args = parser.parse_args()



    if args.load is not None:

        data = pickle.load(open(args.load,'rb'))
        pdb = PDBunique(data['graph'],cluster = None,outfile=args.load)
        pdb.percent = data['percent']
        pdb.unique_pdbs = data['ids']
        pdb.map = data['map']

    else:
        pdb = PDBunique(args.ssg,cluster=args.cluster,outfile=args.outfile)

        if args.cluster is not None:
            graph = pdb.get_protein_cluster_graph(pdb.clusters[args.cluster],pdb.percent)
            plot_graph(graph,'protein_cluster%d' %(args.cluster))
        else:
            pdb.get_unique_entries()
            pdb.save(args.outfile)
