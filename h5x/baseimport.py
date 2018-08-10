#%matplotlib inline
from scrapePDB.PDBsim import PDBsim
from scrapePDB import PDBunique

def connection_graph(grp):

    import plotly.offline as py
    py.init_notebook_mode(connected=True)

    ssg = PDBsim(None)
    ssg.load_hdf5(grp=grp)
    ssg.plot_graph('PDBsim',offline=True,noedge=True,ind_cluster=None)


def protein_graph(grp):

    import plotly.offline as py
    py.init_notebook_mode(connected=True)

    g = PDBunique.load_prot_graph(grp=grp)
    PDBunique.plot_graph(g,'protgraph',offline=True)