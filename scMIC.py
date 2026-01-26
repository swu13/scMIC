# -*- coding: utf-8 -*-
"""
scMIC: MIC identification via unbalanced optimal transport

Author: xxx
Description:
This script identifies metastasis-initiating cells (MICs) by
integrating primary and metastatic tumor scRNA-seq data using
embedding learning and unbalanced optimal transport.

Input:
- Primary expression matrix (gene ¡Á cell)
- Metastatic expression matrix (gene ¡Á cell)

Output:
- Contribution score for each primary tumor cell
"""

import torch
import torch.nn as nn        

class MLP(nn.Module):
    def __init__(
        self,
        input_dim: int,
        output_dim: int,
        hidden_dims: list,
        act_fn=nn.SELU,
        dop: float = 0.1,
        out_dop: float | None = None
    ):
        super(MLP, self).__init__()
        layers = []

        layers.append(
            nn.Sequential(
                nn.Linear(input_dim, hidden_dims[0]),
                act_fn(),
                nn.Dropout(dop)
            )
        )

        for i in range(len(hidden_dims) - 1):
            layers.append(
                nn.Sequential(
                    nn.Linear(hidden_dims[i], hidden_dims[i + 1]),
                    act_fn(),
                    nn.Dropout(dop)
                )
            )

        self.net = nn.Sequential(*layers)

        if out_dop is None:
            self.output_layer = nn.Linear(hidden_dims[-1], output_dim)
        else:
            self.output_layer = nn.Sequential(
                nn.Linear(hidden_dims[-1], output_dim),
                act_fn(),
                nn.Dropout(out_dop)
            )

    def forward(self, x):
        return self.output_layer(self.net(x))



class AutoEncoder(nn.Module):
    def __init__(self, input_dim, latent_dim, hidden_dims):
        super(AutoEncoder, self).__init__()
        self.encoder = MLP(input_dim, latent_dim, hidden_dims, out_dop=0.1)
        self.decoder = MLP(latent_dim, input_dim, hidden_dims[::-1], out_dop=0.1)

    def forward(self, x):
        latent = self.encoder(x)
        reconstructed = self.decoder(latent)
        return reconstructed, latent

  
def ave(tau, u, u1):
    return tau * u + ( 1 - tau ) * u1

def lse(A):
    return np.log(np.sum(np.exp(A),axis=1)).reshape(-1,1)

def H(p):
    return -np.sum( p * np.log(p+1E-20)-1 )

def KL(h,p):
    return np.sum( h * np.log( h/p ) - h + p )

def KLd(u,p):
    return np.sum( p * ( np.exp(-u) - 1 ) )

def M(u,v,H1,H2,c,epsilon):
    y = -c + np.matmul(u.reshape(-1,1), H2.reshape(1,-1)) + \
        np.matmul(H1.reshape(-1,1), v.reshape(1,-1))
    return y/epsilon


import numpy as np
def uot(mu, nu, c, epsilon,niter=50, tau=-0.5, verb = 1, rho = np.inf, stopThr= 1E-7):
    import ot
    lmbda = rho / ( rho + epsilon )
    if np.isinf(rho): lmbda = 1
  
    mu = np.asarray(mu, float).reshape(-1,1)
    nu = np.asarray(nu, float).reshape(-1,1)
    N = [mu.shape[0], nu.shape[0]]
    H1 = np.ones([N[0],1]); H2 = np.ones([N[1],1])
  
    errs = []; Wprimal = []; Wdual = []
    u = np.zeros([N[0],1], float)
    v = np.zeros([N[1],1], float)
    for i in range(niter):
        u1 = u
        u = ave(tau, u, \
            lmbda * epsilon * np.log(mu) \
            - lmbda * epsilon * lse( M(u,v,H1,H2,c,epsilon) ) \
            + lmbda * u )
        v = ave(tau, v, \
            lmbda * epsilon * np.log(nu) \
            - lmbda * epsilon * lse( M(u,v,H1,H2,c,epsilon).T ) \
            + lmbda * v )
        gamma = np.exp( M(u,v,H1,H2,c,epsilon) )
    
        if np.isinf(rho):
            Wprimal.append(np.sum(c * gamma) - epsilon*H(gamma) )
            Wdual.append(np.sum(u*mu) + np.sum(v*nu) - epsilon*np.sum(gamma) )
            err = np.linalg.norm( np.sum(gamma,axis=1) - mu )
            errs.append( err )
        else:
            Wprimal.append(np.sum(c*gamma) - epsilon*H(gamma) \
                           + rho*KL(np.sum(gamma,axis=1), mu) \
                           + rho*KL(np.sum(gamma,axis=0), nu) )
            Wdual.append(- rho*KLd(u/rho,mu) - rho*KLd(v/rho,nu) \
                         - epsilon*np.sum(gamma) )
            err = np.linalg.norm(u-u1,1)
            errs.append( err )
        if err < stopThr or i == niter-1:
            print(err, i, sep="\t", flush=True)
            break
      
    return gamma

  
def set_seed(seed):
    import os, random
    import numpy as np
    import torch
    random.seed(seed)             
    np.random.seed(seed)          
    torch.manual_seed(seed)       
    torch.cuda.manual_seed_all(seed)  
    os.environ["PYTHONHASHSEED"] = str(seed)
    os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8" 
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    torch.use_deterministic_algorithms(True)
  
def setup_logger(log_path):
    import os
    import logging
    os.makedirs(os.path.dirname(log_path), exist_ok=True)    
    logger = logging.getLogger('training_logger')
    logger.setLevel(logging.INFO)   
    # prevent multiple handler
    if not logger.handlers:
        fh = logging.FileHandler(log_path)
        fh.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(message)s')
        fh.setFormatter(formatter)
        logger.addHandler(fh)
        # Optional: also log to console
        ch = logging.StreamHandler()
        ch.setFormatter(formatter)
        logger.addHandler(ch)   
    return logger
  
  
def MICot(PriFile, MetFile, outputFile, intergrate_flag, pcnum, AEdim,distance_flag, ot_flag, n_iter, epsilon, rho, top_k,seed,epochs,hidden_dims):
    import pandas as pd #pip install pandas
    Primatrix = pd.read_csv(PriFile, index_col=0) 
    Metmatrix = pd.read_csv(MetFile, index_col=0)
    
    common_index = Primatrix.index.intersection(Metmatrix.index)
    Primatrix = Primatrix.loc[common_index]
    Metmatrix = Metmatrix.loc[common_index]
    
    Primatrix_suffixed = Primatrix.add_suffix('-P')
    Metmatrix_suffixed = Metmatrix.add_suffix('-M') 
    
    matrixdata = pd.concat([Primatrix_suffixed, Metmatrix_suffixed], axis=1)
    
    if intergrate_flag == "PCA":
        from sklearn.decomposition import PCA
        pca = PCA(n_components=pcnum)
        pca_result = pca.fit_transform(matrixdata.T)
        
        Pri_embedding = pca_result[0:Primatrix.shape[1],]
        Met_embedding = pca_result[Primatrix.shape[1]:pca_result.shape[0],]
    elif intergrate_flag == "AE":
        import anndata as ad
        import scvi #pip install scvi-tools 
        import torch.optim as optim
        import anndata as ad     
        import os
        adata = ad.AnnData(X=matrixdata.T)
        adata.obs_names = matrixdata.columns.astype(str)
        adata.obs_names_make_unique()
        adata.var_names = matrixdata.index.astype(str)
        adata.var_names_make_unique()
        
        X = torch.tensor(matrixdata.T.values, dtype=torch.float32)
        input_dim = X.shape[1]  
        latent_dim = AEdim
        learning_rate = 1e-3   
        set_seed(seed)
        best_loss = float('inf')
        best_model = None
        log_path = os.path.join("./train.log")
        logger = setup_logger(log_path)
        model=AutoEncoder(input_dim, latent_dim, hidden_dims)
        optimizer = optim.Adam(model.parameters(), lr=learning_rate)
        loss_fn = nn.MSELoss()
        for epoch in range(epochs):
            optimizer.zero_grad()
            reconstruction, latent = model(X)
            loss = loss_fn(reconstruction, X)
            loss.backward()
            optimizer.step()
            if loss < best_loss:
                best_loss = loss
                best_model = model.state_dict().copy()
            logger.info(f"Epoch {epoch+1}/{epochs}, Loss: {loss.item():.4f}")    
        model.load_state_dict(best_model)
        with torch.no_grad():
            _, latent = model(X)    
        latent = latent.numpy()
        Pri_embedding = latent[0:Primatrix.shape[1],]
        Met_embedding = latent[Primatrix.shape[1]:latent.shape[0],]  
    else:
        Pri_embedding = Primatrix.T.to_numpy()
        Met_embedding = Metmatrix.T.to_numpy()
    
    if distance_flag == "Euclidean":
        distance = []
        for i in range(Pri_embedding.shape[0]):
            for j in range(Met_embedding.shape[0]):
                distance.append(np.linalg.norm(Pri_embedding[i,:] - Met_embedding[j,:]))
            
        distance_matrix = np.array(distance) 
        distance_matrix = distance_matrix.reshape(Pri_embedding.shape[0], Met_embedding.shape[0])  
        cost_matrix = (distance_matrix - np.min(distance_matrix)) / (np.max(distance_matrix) - np.min(distance_matrix))
        final_cost_matrix = cost_matrix
    elif distance_flag == "Pearson":
        from scipy.stats import pearsonr
        correlations = []
        for i in range(Pri_embedding.shape[0]):
            for j in range(Met_embedding.shape[0]):
                r, p = pearsonr(Primatrix.iloc[:, i], Metmatrix.iloc[:, j])
                correlations.append(r)
      
        cor_matrix = np.array(correlations) 
        cor_matrix = cor_matrix.reshape(Pri_embedding.shape[0], Met_embedding.shape[0])  
        cor_matrix[cor_matrix < 0] = 0
        final_cost_matrix = 1 - cor_matrix
    elif distance_flag == "Spearman":
        from scipy.stats import spearmanr
        correlations = []
        for i in range(Pri_embedding.shape[0]):
            for j in range(Met_embedding.shape[0]):
                r, p = spearmanr(Primatrix.iloc[:, i], Metmatrix.iloc[:, j])
                correlations.append(r)
        
        cor_matrix = np.array(correlations) 
        cor_matrix = cor_matrix.reshape(Pri_embedding.shape[0], Met_embedding.shape[0])  
        cor_matrix[cor_matrix < 0] = 0
        final_cost_matrix = 1 - cor_matrix
    else:
        raise ValueError("Please provide the method to calculate the distance: Euclidean, Pearson, Spearman")
    
    weight_matrix = np.exp(1-final_cost_matrix)
    w_a = np.sum(weight_matrix, axis=1)
    w_b = np.sum(weight_matrix, axis=0)
    w_a = w_a/np.sum(w_a)
    w_b = w_b/np.sum(w_b)
    
    if ot_flag:
        import ot #pip install POT
        gamma = uot(w_a, w_b, final_cost_matrix, epsilon, rho = rho, niter=n_iter)
    else:
        gamma = final_cost_matrix
    
    if top_k == np.inf:
        gamma_result = gamma / gamma.sum(axis=0)  
    else:
        top_k = top_k
        top_k_indices = np.argsort(gamma, axis=0)[-top_k:, :]
        gamma_result = np.zeros_like(gamma)
        for col in range(gamma.shape[1]):
            gamma_result[top_k_indices[:, col], col] = gamma[top_k_indices[:, col], col]
        gamma_result = gamma_result / gamma_result.sum(axis=0)   
        
    contri = gamma_result.sum(axis=1)
    with open(outputFile, 'w') as f:
        for value in contri:
            f.write(f"{value}\n")

  
  
def main():
  import argparse
  parser = argparse.ArgumentParser(description="Run scMIC pipeline")
  parser.add_argument('--pri', type=str, required=True, help='Primary matrix CSV file: gene*cellid')
  parser.add_argument('--met', type=str, default='', help='Metastatic matrix CSV file: gene*cellid')
  parser.add_argument('--out', type=str, required=True, help='Output file to save contribution scores: primary cellid*score')
  parser.add_argument('--integration', type=str, choices=['PCA', 'AE', 'None'], default='PCA', help='Integration method (PCA, AE, or None)')
  parser.add_argument('--pcnum', type=int, default=2, help='Number of PCs')
  parser.add_argument('--AEdim', type=int, default=10, help='hidden dimention of autoencoder')
  parser.add_argument('--distance', type=str, choices=['Euclidean', 'Pearson', 'Spearman'], default='Euclidean', help='Distance method (Euclidean, Pearson or Spearman)')
  parser.add_argument('--ot', action='store_true', help='Use Optimal Transport')
  parser.add_argument('--n_iter', type=int, default=5000, help='Number of OT iterations')
  parser.add_argument('--epsilon', type=float, default=0.5, help='Epsilon for OT')
  parser.add_argument('--rho', type=float, default=100, help='Rho for UOT')
  parser.add_argument('--top_k', type=int, default=1, help='Top k connections to keep')
  parser.add_argument('--seed', type=int, default=123, help='Random seed selection')
  parser.add_argument('--epochs', type=int, default=100)
  parser.add_argument('--hidden_dims', type=int, nargs='+', default=[512, 256, 128, 64, 32])
  
  args = parser.parse_args()
  MICot(
      PriFile=args.pri,
      MetFile=args.met,
      outputFile=args.out,
      intergrate_flag=args.integration,
      pcnum=args.pcnum,
      AEdim = args.AEdim,
      distance_flag=args.distance,
      ot_flag=args.ot,
      n_iter=args.n_iter,
      epsilon=args.epsilon,
      rho=args.rho,
      top_k=args.top_k,
      seed=args.seed,
      epochs=args.epochs,
      hidden_dims=args.hidden_dims           
  )
      

if __name__ == "__main__":
  main()

