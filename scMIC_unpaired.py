# -*- coding: utf-8 -*-
"""
scMIC: MIC identification for unpaired primary-metastasis based on scFoundation

Author: xxx
Description:
This script identifies metastasis-initiating cells (MICs) by
integrating primary and metastatic tumor scRNA-seq data using
embedding learning and unbalanced optimal transport.

Input:
- Primary expression matrix (cell * embedding)
- Reference Primary expression matrix (cell * embedding)
- Reference Primary phenotype matrix (cell X group)

Output:
- Contribution score for each query primary tumor cell
"""

def MICmap(Pri_embedding,Ref_embedding,RefPhenotypeFile,outputFile,neighbork):
    import numpy as np
    import pandas as pd
    from sklearn.neighbors import NearestNeighbors
    Pri_emb = pd.read_csv(Pri_embedding, index_col=0)
    Ref_emb = pd.read_csv(Ref_embedding, index_col=0)
    RefPheno = pd.read_csv(RefPhenotypeFile, index_col=0)
    RefPheno = RefPheno.loc[Ref_emb.index]
    nn = NearestNeighbors(n_neighbors=neighbork, metric="euclidean")
    nn.fit(Ref_emb)
    dist, idx = nn.kneighbors(Pri_emb)
    labelgroup = RefPheno.loc[Ref_emb.index[idx.flatten()], "group"].values.reshape(idx.shape)
    labels = RefPheno.loc[Ref_emb.index[idx.flatten()], "score"].values.reshape(idx.shape)
    mic_counts = (labelgroup == "Y").sum(axis=1)
    nonmic_counts = (labelgroup == "N").sum(axis=1)
    vote1 = np.where(mic_counts > nonmic_counts, "Y", "N")
    vote2 = np.where(labels.mean(axis=1) > 1, "Y", "N")
    out = pd.DataFrame({
        "MIC_vote": vote1,
        "MIC_score": labels.mean(axis=1),
        "MIC_vote_score": vote2
    }, index=Pri_emb.index)
    out.to_csv(outputFile)
   
def main():
  import argparse
  parser = argparse.ArgumentParser(description="Run scMIC unpaired pipeline")
  parser.add_argument('--pri', type=str, required=True, help='Primary embedding file: cell*embedding')
  parser.add_argument('--Refpri', type=str, default='', help='Please input the reference primray embedding file: cell*embedding')
  parser.add_argument('--Refpripheno', type=str, default='', help='Please input the reference phenotype (MIC or non-MIC) of primray scRNA data')
  parser.add_argument('--out', type=str, default='', help='Please input the output file')
  parser.add_argument('--neighbork', type=int, default=5, help='K neighbors for query primary matching to reference primary')
  
  args = parser.parse_args()
  MICmap(
      Pri_embedding=args.pri,
      Ref_embedding=args.Refpri,
      RefPhenotypeFile=args.Refpripheno,
      outputFile=args.out,
      neighbork=args.neighbork,
  )
      

if __name__ == "__main__":
  main()

