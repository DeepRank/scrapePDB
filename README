# PDB database scrapper for dataset creation

## Installation

## PDBselect

PDBselect alows to find the PDb ID of all the entries that respect a set of predetermined conditions. These conditions are so far on the number, type and length of chains.

For example:

```
usage: PDBdatabase scraper [-h] [--start START] [--size SIZE] [--nproc NPROC]
                           [--method METHOD] [--min_res MIN_RES]
                           [--number_of_entity NUMBER_OF_ENTITY]
                           [--types TYPES [TYPES ...]] [--len_min LEN_MIN]
                           [--len_max LEN_MAX] [--tqdm TQDM]

optional arguments:
  -h, --help            show this help message and exit
  --start START         start index for the search
  --size SIZE           Number of pdbs to screen
  --nproc NPROC         Number of concurrent procs to use
  --method METHOD       characterisation method
  --min_res MIN_RES     minimum resolution
  --number_of_entity NUMBER_OF_ENTITY
                        number of entities
  --types TYPES [TYPES ...]
                        type of polymers
  --len_min LEN_MIN     Minimum number of residues
  --len_max LEN_MAX     Maximum number of residues
  --tqdm TQDM           use tqdm to monitor progress
```

To select the entries you can specify

```
method           : method used to determine the structure (xray)
min_res          : minimum resolution desired (None)
number_of_entity : number of entities in the structure (2)
types            : type of the entity (protein)
len_min          : minimum number of residue in each entity (50)
len_max          : maximu number of residue in each entity (Inf)
```

Other options allows to control the calculations

```
start            : index of the first entry considered (0)
size             : numner of entries processed (all)
nproc            : number of processes used (1)
```

The results will be stored in a pickle file containing the search parmaters and the resulting pdb IDs.