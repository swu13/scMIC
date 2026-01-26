# scMIC

## Introduction

We propose an optimal transport¨Cbased framework, `scMIC`, for identifying metastasis-initiating cells (MICs) from primary tumor single-cell transcriptomic data.

# Repository Structure

The repository is organized as follows:

1. `Processing.R`£º Utility functions for preprocessing single-cell RNA-seq data.

2. `scMIC.py`: Identifies MICs in primary tumors using paired primary and metastatic single-cell samples.

3. `GeneProgram.py`: Infers MIC-related latent representations based on single-cell data and MIC labels.

4. `scMIC_unpaired.py`: Identifies MICs in query primary tumor cells using a reference primary tumor dataset with known MIC labels (unpaired setting).

5. `Figure2`: Code for data downloading, preprocessing, and analysis used to generate Figure 2 in the manuscript.

6. `Figure3`: Code for data downloading, preprocessing, and analysis used to generate Figure 3 in the manuscript.

7. `Figure4`: Code for data downloading, preprocessing, and analysis used to generate Figure 4 in the manuscript.

8. `Figure5`: Code for data downloading, preprocessing, and analysis used to generate Figure 5 in the manuscript.

9. `Figure6`: Code for data downloading, preprocessing, and analysis used to generate Figure 6 in the manuscript.

10. `Discussion`: Code for data downloading, preprocessing, and analysis used in the Discussion section of the manuscript.

## Overview of scMIC
![image](https://github.com/swu13/scMIC/blob/main/scMIC.png)

## Citation

If you use scMIC in your research, please cite:

[1] Sijia Wu, Jiangpeng Wei, Xin Liu, Jiajin Zhang, Jianguo Wen, Liyu Huang, Xiaobo Zhou, Identification and Characterization of Metastasis-initiating cells
