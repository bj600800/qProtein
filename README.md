<h1 align="center">
  <img src="logo-no-background.png" alt="Logo" width="200">
</h1>

**qProtein: A platform for structural proteomics quantitative analysis.** 

[![Static Badge](https://img.shields.io/badge/ACS_JCIM-10.1021%2Facs.jcim.4c01303-red)](https://pubs.acs.org/doi/full/10.1021/acs.jcim.4c01303)
[![DOI](https://zenodo.org/badge/608680173.svg)](https://doi.org/10.5281/zenodo.14214907)



## Key features
:sparkles: **Ingrated platform**: An automated workflow integrating a variety of bioinformatics tools and a deep learning model.

:sparkles: **Batch Processing**: Transform multiple protein sequences to structures.

:sparkles: **Interaction Analysis**: Compute four major types of physical interactions: hydrogen bonds, electrostatic interactions, hydrophobic clusters, and disulfide bonds.

:sparkles: **Overall and Local Structural Analysis**: Analyze local structural features within targeted regions.

:sparkles: **The useful 6 tools**: Fetch, Overall, Local, Visual, Surface, Landscape 

## Install qProtein

```Bash
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
### Fetch structures ###
```Bash
python run_qprotein.py --mode fetch --work_dir test --id test/id.txt 
```

### Overall analysis ###
```Bash
python run_qprotein.py --mode overall --work_dir test --pre_pdb pdb_dir
```

### Local analysis ###
```Bash
python run_qprotein.py --mode local --work_dir test --pre_pdb pdb_dir --template_name P33557 --template_active_res 33,35,37,64,66,91,93,97,99,106,108,115,116,118,142,146,147,148,154,156,158,191,197,199,200 --dist1 12 --dist2 15
```

### visual analysis ###
```Bash
python run_qprotein.py --mode visual --work_dir test --pre_pdb pdb_dir --pml
```

### surface analysis ###
```Bash
python run_qprotein.py --mode surface --work_dir test --pre_pdb pdb_dir --config_file config.ini
```

### function landscape analysis ###
```Bash
python run_qprotein.py --mode landscape --work_dir test --pre_pdb pdb_dir --label_file abs_label_txt --template_name P33557 --positions 141,142,143,144 --label_pos_thresold 50 --config_file config.ini
```

## Citation
```
@article{dou2024qprotein,
  title={qProtein: Exploring Physical Features of Protein Thermostability Based on Structural Proteomics},
  author={Dou, Zhixin and He, Jiaxin and Han, Chao and Wu, Xiuyun and Wan, Lin and Yang, Jian and Zheng, Yanwei and Gong, Bin and Wang, Lushan},
  journal={Journal of Chemical Information and Modeling},
  volume={64},
  number={20},
  pages={7885--7894},
  year={2024},
  publisher={ACS Publications}
  doi={10.1021/acs.jcim.4c01303}
}
```

## Online web-server (for ID mode only with seq2struct)
http://qprotein.sdu.edu.cn:8888

## Schematic diagram of interaction algorithm
<img src="https://github.com/bj600800/qProtein/blob/main/interaction_algorithm.png" alt="algorithm" width="650" height="600">


