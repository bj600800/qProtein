# qProtein：Exploring physical interaction features of protein thermostability based on structural proteomics

**A platform for automated structural proteomics analysis.**

## Online web-server
https://qprotein.sdu.edu.cn

## Key features
:sparkles: **Batch Processing**: Transform multiple protein sequences to structures.

:sparkles: **Interaction Analysis**: Compute four major types of physical interactions: hydrogen bonding, electrostatic interactions, hydrophobic clusters, and disulfide bond.

:sparkles: **Overall and Local Structural Analysis**: Analyze local structural features within targeted regions of interest following multiple structure alignment.

## Install qProtein

```
git clone https://github.com/bj600800/qProtein.git

cd qProtein

conda create -n qprotein python=3.9
conda activate qprotein

pip install biotite==0.37.0 
pip install pdb2pqr
pip install tqdm
pip install numpy==1.23.5
pip install biopython==1.83
```

## Additional requirements
**Notice:**
Modify config params in run_qprotein.py for your installation.

### Structural alignment for local structure analysis ###
For Win
[USalign.exe](https://zhanggroup.org/US-align/bin/module/USalignWin64.zip)

For Linux
[USalign](https://zhanggroup.org/US-align/bin/module/USalignLinux64.zip)

### For prediction module ###
[ESMFold](https://github.com/facebookresearch/esm)


## Run qProtein
### Get structures ###
```
python run_qprotein.py --id test/id.txt --dir test
```

### Using user-prepared structures for overall analysis ###
```
python run_qprotein.py --dir test
```

### Local analysis ###
```
python run_qprotein.py --id test/id.txt --dir test --local --template_name P33557 --template_active_res 33,35,37,64,66,91,93,97,99,106,108,115,116,118,142,146,147,148,154,156,158,191,197,199,200 --dist1 12 --dist2 15
```

## Schematic diagram of interaction algorithm
![Interaction algorithm](https://github.com/bj600800/qProtein/blob/main/interaction_algorithm.png)

