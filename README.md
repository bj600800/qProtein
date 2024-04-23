# qProtein: Leveraging evolutionary analysis of protein structures to explore the relationship between structure and function.

A platform for automated structural proteomics analysis.
![Interaction algorithm](https://github.com/bj600800/qProtein/blob/main/interaction_algorithm.png)

## Key features
:rocket: **Batch Processing**: Transform multiple protein sequences to structures.

:sparkles: **Interaction Analysis**: Compute four major types of physical interactions: hydrogen bonding, electrostatic interactions, hydrophobic clusters, and disulfide bond.

:heart: **Local Structural Analysis**: Analyze local structural features within targeted regions of interest following multiple structure alignment.

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
