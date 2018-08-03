# PDB database scrapper for dataset creation

## Installation

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
./PDBselect --method xray --types protein --number_of_entity 2 --len_min 10 --len_max 500 --nproc 10
``

selects all the PDB containging 2 proteins containing between 10 and 500 residues each and charaterized with X-ray. The program write a pickle file containing the search parameter and the the corresponding PDB IDs.

The last option `nproc` specfied that 10 processes will be used in parallel to process the data.

## PDBsim

Many of the PDBs retruned by `PDBselect` contains chains that are somehow idnetical. `PDBsim` allows to identify al the entries that have one chain in common. The level of similarity can be specified. For example:

```
./PDBsim --xtfile pdblist.pkl --percent 40 --nproc 10
```

will identfy all the pdbs in `pdblist.pkl` that contains one common chain with a sequence similarity cutoff of 40 percent. Here as well `nproc` specifies that 10 processess will be used in parallel to process the data.

The program outputs a `networkX` graph file where each node corresponds to a given PDB ID and where nodes are connected if they share a common chain. This graph can be vizualized with the embedded `plotly` script.


![alt-text](./seqsim.gif)