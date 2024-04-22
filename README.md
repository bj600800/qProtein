# qProtein
An automated platform for the analysis and quantification of protein structure-omics.

## Key features
:sparkles: **Batch Processing**: Transform multiple protein sequences simultaneously, saving you time and effort.
:sparkles: **Comprehensive Analysis**: Assess four major types of physical interactions: hydrogen bonding, electrostatic interactions, hydrophobic clusters, and disulfide bond.
:sparkles: **Local Structural Analysis**: Local structural features from specific regions of interest after performing multiple structural alignments.

## Requirements
- gcc (if missing)
  - for Ubuntu: sudo apt install gcc

## Installing and Building from source

```
git clone https://github.com/bj600800/qProtein.git

cd qProtein

chmod 775 -R ./

conda env create -f environment.yml

./run_qprotein.sh
```
