# qProtein: Leveraging evolutionary analysis of protein structures to explore the relationship between structure and function.

**A platform for automated structural proteomics analysis.**

## Key features
:sparkles: **Batch Processing**: Transform multiple protein sequences to structures.

:sparkles: **Interaction Analysis**: Compute four major types of physical interactions: hydrogen bonding, electrostatic interactions, hydrophobic clusters, and disulfide bond.

:sparkles: **Overall and Local Structural Analysis**: Analyze local structural features within targeted regions of interest following multiple structure alignment.

## Additional requirements
For local structure analysis:
[USalign](https://zhanggroup.org/US-align/)

For prediction module:
[ESMFold](https://github.com/facebookresearch/esm)

## Installing and Building from source

```
git clone https://github.com/bj600800/qProtein.git

cd qProtein

conda env create -f environment.yml

python run_qprotein.py
```

## Schematic diagram of interaction algorithm
![Interaction algorithm](https://github.com/bj600800/qProtein/blob/main/interaction_algorithm.png)

