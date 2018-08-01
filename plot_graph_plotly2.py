import plotly.plotly as py
import plotly.graph_objs as go
import networkx as nx

def plot_graph(G,pos,fname):


    trace3_list = []
    middle_node_trace = go.Scatter(x=[],y=[],text=[],mode='markers',
                           hoverinfo='text',marker=go.Marker(opacity=0))

    # edge_trace = go.Scatter(
    #     x=[],
    #     y=[],
    #     line=dict(width=0.5,color='#888'),
    #     hoverinfo='text',
    #     mode='lines',
    #     text = [])

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
            showscale=True,
            # colorscale options
            # 'Greys' | 'Greens' | 'Bluered' | 'Hot' | 'Picnic' | 'Portland' |
            # Jet' | 'RdBu' | 'Blackbody' | 'Earth' | 'Electric' | 'YIOrRd' | 'YIGnBu'
            colorscale='YIGnBu',
            reversescale=True,
            color=[],
            size=10,
            colorbar=dict(
                thickness=15,
                title='Node Connections',
                xanchor='left',
                titleside='right'
            ),
            line=dict(width=2)))

    for node in G.nodes():
        x, y = pos[node]
        node_trace['x'].append(x)
        node_trace['y'].append(y)
        node_trace['text'].append(''.join(node) + '<br>'+ G.nodes[node]['txt'])
        node_trace['marker']['color'].append(G.nodes[node]['number'])

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