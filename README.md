# PDB database scrapper for dataset creation

## Introduction

Clone the repository and install it with pip

```
git clone https://github.com/DeepRanl/scrapePDB
cd scrapePDB
pip install -e ./
```

If you want to be able to use the command line interface everywhere, add the repo to your Python path:

```
export $PYTHONPATH='/path/to/scrapePDB/scrapePDB/:$PYTHONPATH'
```

The creation of a data set is a four step process detailled in the following. The dataset is stored in a HDF5 file whose content can be vizualized using  `h5xplorer`. The four steps are:


#### 1 PDBselect : Select Potential PDB

Here we simly select all the PDBs that respect a set of conditions (number, types, length of entities in the strucrure)

#### 2 PDBsim : Identify PDBs with a common entity

Here we identify all the PDBs that share a common entity. This is done via the sequence similarity of the chains and is therefore defined by the sequence percentage below which two chains are considered idnetical

#### 3 PDBunique : Identify unique PDBs 

Here we go through all PDBs sharing a common entity and extract unique PDBs entries. This is done by creating the protein interaction graph and selecting a single PDB for each edge of the this graph.

#### 4 PDBdowload : download the dataset

Just dowload the unique PDBs we have selected


## PDBselect

PDBselect alows to find the PDb ID of all the entries that respect a set of predetermined conditions. These conditions are so far on the number, type and length of chains. PDB can be selected with the following options:

```
method           : method used to determine the structure (xray)
min_res          : minimum resolution desired (None)
number_of_entity : number of entities in the structure (2)
types            : type of the entity (protein)
len_min          : minimum number of residue in each entity (50)
len_max          : maximum number of residue in each entity (Inf)
```
For example:

```
./PDBselect --method xray --types protein --number_of_entity 2 --len_min 10 --len_max 500 --nproc 10 --hdf5 dataset.hdf5
```

selects all the PDB containging 2 proteins containing between 10 and 500 residues each and charaterized with X-ray. The program write a pickle file containing the search parameter and the the corresponding PDB IDs.

The last option `nproc` specfied that 10 processes will be used in parallel to process the data.

The selected PDB IDs are then stored in the HDF5 file specified together with all the parameter search.

## PDBsim

Many of the PDBs retruned by `PDBselect` contains chains that are somehow idnetical. `PDBsim` allows to identify al the entries that have one chain in common. The level of similarity can be specified. For example:

```
./PDBsim dataset.hdf5 --percent 40 --nproc 10
```

will identfy all the pdbs selected in `dataset.hdf5` that contains one common chain with a sequence similarity cutoff of 40 percent. Here as well `nproc` specifies that 10 processess will be used in parallel to process the data.

The corresponding connection graph will be stored in `dataset.hdf5`. In this graph each node corresponds to a given PDB ID and two nodes are connected if they share a common chain. This graph can be vizualized with the `h5xplorer` browser locate at `h5x/h5x.py`/. using thos browser open `dataset.hdf5` and right click on the `PDsim` folder. Click on `Connection Graph` to launch a plot.ly visualization in your web browser.

![alt-text](./seqsim.gif)

An online version of scuh a graph can be visualized here : https://plot.ly/~nicoreno/20/ as an example.

As you can see al the PDBs are grouped in clusters of entries that share at least one chain. These clusters can be further analyzed to extract unique PDB Ids.

## PDBunique

Now that we know which PDBs share a common cluster we can select unique pdbs with `PDBunique`. This program analyzes the cluster of PDB sharing a common chain and extract only the entries that are different from each other. To do so simply type:

```
./PDBunique dataset.hdf5
```

This will compute the protein interation graph of all the clusters, select unique PDB IDs and store all that in 'dataset.hdf5'. The result can be visualized with the `h5x.py` browser. Open `dataset.hdf5`, click on the PDBunique folder and right click on one of the cluster subfolder. Choose `Protein Graph` to lauch a plot.ly visualization in your web browser.

![alt-text](./protclust.gif)

In this graph each node represent a single chain and each edge shows the interaction of two distinct chains. Hovering over the nodes will show the name of the chains and the PDB entries where it is found. Hovering over the edges will show the PDB Ids that contains both chains the edge is linking. A copy of such a graph can be viewed here : https://plot.ly/~nicoreno/22/ for example.

The unique PDB Ids are then stored in `dataset.hdf5`, in the PDBunique folder and in the `ids` dataset.


## PDBdownload

We just need now to download the PDBs we have selected. This can be done via

```
./PDBdownload dataset.hdf5 --outdir ./dataset/
```