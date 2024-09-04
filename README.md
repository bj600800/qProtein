<h1 align="center">
  <img src="logo-no-background.png" alt="Logo" width="200">
</h1>

**A platform for automated structural proteomics analysis.**

## Key features
:sparkles: **Ingrated platform**: An automated workflow integrating a variety of bioinformatics tools and a deep learning model.

:sparkles: **Batch Processing**: Transform multiple protein sequences to structures.

:sparkles: **Interaction Analysis**: Compute four major types of physical interactions: hydrogen bonds, electrostatic interactions, hydrophobic clusters, and disulfide bonds.

:sparkles: **Overall and Local Structural Analysis**: Analyze local structural features within targeted regions.

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

### Configure ESMFold for prediction module ###
[ESMFold](https://github.com/facebookresearch/esm)


## Run qProtein
### Getting structures ###
```
python run_qprotein.py --id test/id.txt --dir test
```

### Using user-prepared structures for overall analysis ###
**Notice:**
PDB structures shoud be included in the directory named "structure" under the working directory "--dir". Here it refers to the "test" fold.

So that qProtein will analyze the structures included in "structure/" fold, and output the results in the working directory "test/" fold.

```
python run_qprotein.py --dir test
```

### Local analysis ###
**Notice:**
Four params --template_name, --template_active_res --dist1, and --dist2 should be defined by users depending on the task requirement. 

For further details, please refer to the paper.

```
python run_qprotein.py --id test/id.txt --dir test --local --template_name P33557 --template_active_res 33,35,37,64,66,91,93,97,99,106,108,115,116,118,142,146,147,148,154,156,158,191,197,199,200 --dist1 12 --dist2 15
```

### Using ESMFold model for predicting protein structures ###
python run_qprotein.py --fasta test/sequence.fasta --dir test --local --template_name P33557_seq --template_active_res 33,35,37,64,66,91,93,97,99,106,108,115,116,118,142,146,147,148,154,156,158,191,197,199,200 --dist1 12 --dist2 15

## Citation
qProtein: Exploring the physical features of protein thermostability based on structural proteomics

## Online web-server (for using Alphafold database only)
http://qprotein.sdu.edu.cn:8888

## Schematic diagram of interaction algorithm
<img src="https://github.com/bj600800/qProtein/blob/main/interaction_algorithm.png" alt="algorithm" width="650" height="600">


