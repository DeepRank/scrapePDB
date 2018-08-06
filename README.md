# PDB database scrapper for dataset creation

## Introduction

No installation is required for the moment.

The creation of a dataset is a 4 step process detailled in the following:

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
./PDBselect --method xray --types protein --number_of_entity 2 --len_min 10 --len_max 500 --nproc 10
```

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

This graph can be visualized here : https://plot.ly/~nicoreno/20/

## PDBunique

Now that we know which PDBs share a common cluster we can select unique pdbs with `PDBunique`. This program analyzes the cluster of PDB sharing a common chain and extract only the entries that are different from each other. This is achieved by computing the protein cluster graphs as shown below.

![alt-text](./protclust.gif)

In this graph each node represent a single chain and each edge shows the interaction of two distinct chains. Hovering over the nodes will show the name of the chains and the PDB entries where it is found. Hovering over the edges will show the PDB Ids that contains both chains the edge is linking. To generate this graph you can use:

```
./PDBunique graph40.pkl --cluster 0
```


The graph can be viewed here : https://plot.ly/~nicoreno/22/

If you do not specify the cluster index the program will go through all the clusters and extract unique PDB IDs corresponding to single edges in the protein graph of each cluster. This will result in a `pickle` file containing the IDs of the unique PDBs. For example selecting the smallest PDB for each edge of the protein graph leads to the 778 following PDBs:

```
'5OCL', '1A7Q', '3K1K', '5E0Q', '3P9W', '5GXB', '4N9O', '5IP4', '5F1K', '5BOP', '1ZVH', '3ETB', '3K74' '5F72', '2XV6', '5IML', '4WEU', '4DK3', '1ZA6', '1MVF', '5F8Q', '5VXK', '4W6X', '1KXV', '5DHX', '5OMN' '3UX9', '5AAM', '4KSD', '6DBF', '4W2Q', '4LHQ', '4POU', '4NC0', '5DFW', '5IMM', '3H3B', '4NBZ', '4NIK' '4N1H', '4NC2', '5L21', '5LHR', '5H8O', '5C2U', '3K7U', '5M94', '2Z8W', '2I25', '2Z35', '4XVP', '5MA4' '5AQB', '1WQ9', '3V2A', '1FLT', '4JW3', '4ZV6', '4J4L', '2UZX', '2DZO', '3AJI', '1SVX', '4K5B', '4HRL' '1BI7', '2V5Q', '5MBM', '2TLD', '1D6R', '1ACB', '1SHY', '1K9O', '1HCG', '1BRC', '2QN5', '1C9T', '4AFU' '1CGJ', '1FLE', '5JSN', '5JSB', '1OC0', '5J7C', '4ML7', '3G3B', '4GLA', '1UUZ', '3TNI', '1BUH', '2W9Z' '5UQ2', '1JOW', '1H4L', '2JGZ', '1W98', '4UN0', '4KRD', '3CBK', '3K9M', '2OUL', '3C7U', '3QHY', '4O4B' '5ECJ', '4JE4', '6BYN', '5N7E', '1BHH', '1GGP', '3ZL7', '4PER', '1STF', '3OWE', '5JMC', '1R3H', '1CD1' '1PQZ', '2X83', '2WLV', '1X1X', '1AY7', '1B34', '1B8M', '2WWX', '2OT3', '6AXF', '1BKD', '5CM8', '2EFD' '3TW8', '1BPL', '1CAV', '1CPB', '2NU1', '1D4V', '1KFU', '5OVN', '1EAY', '1EM8', '4P5K', '1J7D', '2GRR' '5TUT', '2BF8', '5FQ2', '4DHI', '6CP0', '4IP3', '3ONG', '4WZ3', '3RCZ', '2XWU', '2O25', '1P3Q', '1YD8' '5HPL', '2OOB', '2HTH', '3TMP', '1XT9', '5ZQ3', '4I6L', '3K9O', '5LN1', '5KYD', '5J26', '3ZNH', '2QHO' '5JP3', '2D07', '2G4D', '1EUV', '5EQL', '2IO1', '1F02', '1US8', '1F3U', '1F3V', '2B7C', '1FM2', '3WSO' '1FQV', '5K35', '2E31', '5HYW', '1FX0', '3OUW', '1LUJ', '1QZ7', '1TH1', '1G3J', '1GH6', '1GO3', '1H6K' '1H3O', '1WQJ', '1H59', '1HL6', '1HWH', '1I1R', '1IAR', '1RY7', '4J23', '2FDB', '1ITB', '1IYJ', '3T5X' '1Y6N', '1JEN', '1JEQ', '1JKG', '1JQL', '1OZ7', '1KA9', '1KAC', '2NTS', '6AT6', '3EVS', '5HK5', '1KTZ' '1KVE', '2NZ8', '4XH9', '5FI1', '1XCG', '2YIN', '5UPL', '3GCG', '1UGH', '3WDG', '1M10', '1MBV', '1MHM' '3ME0', '1N0L', '1PDK', '1N1J', '1N8S', '1PGW', '1O6O', '1OC9', '1OEY', '1OQS', '1P5V', '2CO6', '1PK1' '1QA9', '1QAV', '1QDL', '1R8O', '4E18', '1RKE', '3H2V', '1SC5', '1SC4', '1VET', '1SPP', '1SV0', '1SYX' '1TFO', '1TNR', '1TTW', '1TUE', '2HKQ', '5EW5', '4UHP', '5D5K', '5CTT', '5T94', '3TPM', '1UN0', '3TJ3' '5TBK', '4UAD', '5H2W', '4XZR', '1US7', '1USV', '1V1P', '1VRS', '1WMH', '1XB2', '1XDT', '1XEW', '1XG2' '1XL3', '1XOU', '1Y8X', '3FN1', '1Y96', '1YKH', '1Z92', '1ZBX', '1ZM3', '2A1A', '2B42', '2BAP', '2BW3' '2C0L', '2C35', '2CKZ', '2D5R', '2DBU', '2FHZ', '2F9Z', '2FUN', '2G16', '2GAF', '3OED', '2GOX', '2GWF' '2H7Z', '2VSK', '3GXU', '3HEI', '2HSN', '2HSM', '2I32', '2I3T', '2J04', '2MAD', '3NVM', '2NPT', '2O2V' '2NXN', '2O4X', '2O8A', '5IOH', '2OEN', '4XWJ', '2OMZ', '2OZA', '3TG1', '2P24', '2PM9', '2PQA', '2PQN' '2Q2E', '2QK7', '2QKW', '3HGK', '2QSF', '2QVS', '5X3F', '4X6Q', '2R1V', '2RD0', '2RF4', '3D7U', '4I21' '2UY1', '2V1Y', '2V9T', '2VT1', '2W1T', '2W2K', '2W81', '2WD5', '2WUS', '2WWM', '2X0B', '2X11', '2X9A' '2XHE', '2Y9P', '2Y9Y', '2YJN', '2YSU', '5MS2', '2Z0D', '3VVW', '2Z3Q', '2Z5B', '2ZIX', '2ZNL', '2ZS6' '3A2F', '3A98', '3ABE', '3ANW', '3AQF', '3AU4', '3AXJ', '3B1S', '3BDW', '3BS5', '5N48', '3BX7', '3BZV' '4JEU', '3CQG', '3D3C', '3DGP', '3G9V', '3DLQ', '3E0J', '3E20', '5DMQ', '3EI4', '4A11', '3EYJ', '4H6J' '3F62', '3FBI', '3FHC', '3FXE', '3GB8', '3GC3', '3GFK', '3GQB', '3HVE', '3HW2', '3IEC', '3IF8', '3JRO' '3JV4', '3K1I', '3K1R', '3K51', '3K5B', '3K75', '3K8P', '3KBT', '3KF6', '3KF8', '3KMU', '3KYI', '3L0X' '3L10', '3LBX', '3LF4', '3LKX', '3M7F', '3MCA', '3MKR', '3MP7', '3MV3', '3NBH', '3O6B', '3OGI', '3ONA' '3OUR', '3PH0', '3PHF', '4DSS', '3PIN', '4KT6', '3Q4H', '3QBQ', '3QC8', '4KDI', '5X4L', '3QHE', '3QMZ' '3R0E', '3REA', '3RNQ', '3S4W', '4WNF', '3SF4', '4QRM', '3TBI', '3TL8', '3TLG', '4AKX', '3TUF', '3U82' '3U9T', '3UM3', '3V3K', '3V6B', '3V89', '3VLB', '3VLF', '3VPJ', '3VX7', '3VX8', '3VZ9', '3W03', '5H2V' '3W3W', '3W8I', '3WCY', '3WOE', '3WW0', '3WXE', '5OE7', '3X29', '3YGS', '4ILH', '3ZHE', '3ZVQ', '4A1G' '4A5U', '4A9A', '4ACV', '4AKR', '4APF', '5NLB', '4HXI', '4AXG', '4AYZ', '4B6H', '4B8A', '4B93', '4BH6' '4BI8', '4BIK', '5BMU', '4BVY', '5A5H', '4BSZ', '4C0H', '4C5H', '4CDK', '4CLQ', '4CT6', '4CZD', '4D73' '4DBG', '5X0W', '4DEX', '4DFC', '4EQ6', '4EQA', '4EYY', '4ZP4', '4F3L', '4F7G', '4FHM', '4FQ0', '4G8X' '4G6T', '4G7X', '4GI3', '4GVB', '4H5S', '4HDN', '4HGK', '4HPL', '4HR1', '4HWI', '4IRV', '4IW4', '4IYP' '4JE3', '4JHP', '4JOI', '4JSN', '4K0V', '4K12', '4KBM', '4KBQ', '4KF5', '4KT3', '4LI2', '4QXF', '4LLD' '4LLO', '4LO8', '5CHL', '4M6B', '4MMR', '4N7Z', '4NAF', '4NE3', '4NOO', '4NSO', '4NTQ', '4O9P', '4OH8' '4ONS', '4PJU', '4Q4J', '4QLP', '4RWS', '4S0S', '4TQ0', '4TU3', '4TXV', '4U0Q', '4U1C', '4U4P', '4WLQ' '4UF5', '4UZZ', '4V2C', '4W4K', '4WP5', '4WRM', '4WVM', '4WWI', '5KWG', '4XA9', '4XCI', '4XD9', '4XHR' '5FXY', '4XYI', '5WJC', '4Y1R', '4Y5O', '4Y66', '4YII', '4YK8', '4YVQ', '4YXC', '4ZGQ', '4ZGZ', '4ZOQ' '6FC0', '5ABV', '5NVK', '5AFQ', '5AJD', '5B2G', '5BXF', '5BY8', '5C6G', '5CSC', '5CTR', '5CW5', '5CY5' '5D6H', '5M5P', '5DCM', '5DHM', '5DMB', '5DMR', '5DO7', '5DOI', '5DQS', '5DSE', '5DUD', '5E6P', '5EN6' '5EQJ', '5ET0', '5F22', '5F3J', '5F3Y', '5F5S', '5F5T', '5G37', '5G4K', '5GLE', '5GNA', '5GPI', '5H3W' '5HFT', '5HKQ', '5HXG', '5HY3', '5I0Q', '5HYP', '5IM0', '5IWB', '5IXI', '5J3Y', '5J4A', '5J8Y', '5JFZ' '5JJA', '5JW9', '5K7M', '5KP6', '5L0W', '5LON', '5LSL', '5LUQ', '5LZ3', '5MAW', '5MMZ', '5MU7', '5N3U' '5N8A', '5NBT', '5NKZ', '5NT1', '5OEN', '5OQQ', '5OQR', '5OWU', '5T87', '5U1S', '5UN6', '5V8W', '5VIP' '5WEZ', '5WPA', '5WSU', '5WUR', '5WWL', '5WXL', '5WYL', '5WZF', '5X5W', '5XN6', '5XOC', '5XRW', '5XV1' '5YAH', '5YVI', '5ZRZ', '6ALX', '6AW2', '6B12', '6CH3', '6DXO', '6ET1', '6F3Z', '6F7S', '6F9N', '6FUD'
```

## PDBdownload

We just need now to download the PDBs we have selected. This can be done via

```
./PDBdownload pdb_unique.pkl --outdir ./dataset/
```