# -*- coding: utf-8 -*-
"""
GeneProgram.py: Find gene program of MIC

Author: xxx
Description:
This script identifies MIC latent and MIC gene program

Input:
- Primary expression matrix (gene ¡Á cell)
- MIC and non-MIC labels

Output:
- Latent and F matrix
- ROC and Latent Heatmap
"""

import os
import copy
import torch
import logging
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from matplotlib.patches import Patch
from torch.utils.data import Dataset, DataLoader
from torch import nn, optim
import torch.backends.cudnn as cudnn
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, roc_curve

# ---------------------------- Dataset Class --------------------------------
class GroupedDataset(Dataset):
    def __init__(self, expr, pheno):
        self.expr = expr
        self.pheno = pheno
        self.labels = pheno.values

    def __len__(self):
        return len(self.expr)

    def __getitem__(self, idx):
        return (torch.tensor(self.expr[idx], dtype=torch.float32),
                torch.tensor(self.labels[idx], dtype=torch.float32))

# ---------------------------- Model Definition --------------------------------
class MLP(nn.Module):
    def __init__(self, input_dim, output_dim, hidden_dims, act_fn=nn.SELU, dop: float = 0.1, out_dop=None):
        super().__init__()
        self.dop = dop
        modules = [nn.Sequential(nn.Linear(input_dim, hidden_dims[0], bias=True), act_fn(), nn.Dropout(dop))]
        for i in range(len(hidden_dims) - 1):
            modules.append(nn.Sequential(nn.Linear(hidden_dims[i], hidden_dims[i + 1], bias=True), act_fn(), nn.Dropout(dop)))
        self.net = nn.Sequential(*modules)
        if out_dop is None:
            self.output_layer = nn.Linear(hidden_dims[-1], output_dim, bias=True)
        else:
            self.output_layer = nn.Sequential(nn.Linear(hidden_dims[-1], output_dim, bias=True), act_fn(), nn.Dropout(dop))

    def forward(self, x):
        embed = self.net(x)
        return self.output_layer(embed)

class AutoEncoder(nn.Module):
    def __init__(self, input_dim, latent_dim, hidden_dims):
        super(AutoEncoder, self).__init__()
        self.encoder = MLP(input_dim, latent_dim, hidden_dims, out_dop=0.1)
        self.decoder = MLP(latent_dim, input_dim, hidden_dims[::-1], out_dop=0.1)

    def forward(self, x):
        latent = self.encoder(x)
        reconstructed = self.decoder(latent)
        return reconstructed, latent

class PhenotypeClassifier(nn.Module):
    def __init__(self, latent_dim, phenotype_dim, hidden_dims,organ_flag):
        super(PhenotypeClassifier, self).__init__()
        self.organ_flag = organ_flag
        if organ_flag:
            self.net_organ = MLP(latent_dim, phenotype_dim - 1, hidden_dims)
        else:
            self.net_organ = None
        self.net_mic = MLP(latent_dim, 1, hidden_dims)

    def forward(self, latent):
        mic_logits = self.net_mic(latent)
        if self.organ_flag:
            organ_logits = self.net_organ(latent)
            return organ_logits, mic_logits
        else:
            return None, mic_logits


# ---------------------------- Utilities --------------------------------
def estimate_F(latent, x):
    latent_pinv = torch.linalg.pinv(latent)
    return latent_pinv @ x

def compute_pos_weight(labels):
    pos = labels.sum(dim=0)
    neg = labels.shape[0] - pos
    return neg / (pos + 1e-6)

def setup_logger(log_path):
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    logger = logging.getLogger('training_logger')
    logger.setLevel(logging.INFO)
    if not logger.handlers:
        fh = logging.FileHandler(log_path)
        fh.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(message)s')
        fh.setFormatter(formatter)
        logger.addHandler(fh)
        ch = logging.StreamHandler()
        ch.setFormatter(formatter)
        logger.addHandler(ch)
    return logger

def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)

    cudnn.deterministic = True
    cudnn.benchmark = False
    os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"
    if hasattr(torch, "use_deterministic_algorithms"):
        torch.use_deterministic_algorithms(True)

# ---------------------------- Training Functions --------------------------
# AutoEncoder Train
def pretrain_autoencoders(
    autoencoders,
    dataloader,
    epochs,
    lambda_recons,
    learning_rate,
    log_path,
    device
):
    autoencoders = autoencoders.to(device)
    autoencoders.train()
    optimizer = optim.Adam(autoencoders.parameters(), lr=learning_rate)
    criterion = nn.MSELoss()
    logger = setup_logger(log_path)
    best_loss = float("inf")
    best_autoencoders = copy.deepcopy(autoencoders)
    for epoch in range(epochs):
        epoch_recon_loss = []
        for batch_x, _ in dataloader:
            batch_x = batch_x.to(device)
            optimizer.zero_grad()
            reconstructed, _ = autoencoders(batch_x)
            reconstruction_loss = criterion(reconstructed, batch_x)
            total_loss = lambda_recons * reconstruction_loss
            total_loss.backward()
            optimizer.step()
            epoch_recon_loss.append(reconstruction_loss.item())
        mean_recon = np.mean(epoch_recon_loss)
        mean_total = lambda_recons * mean_recon
        if mean_total < best_loss:
            best_loss = mean_total
            best_autoencoders = copy.deepcopy(autoencoders)
        logger.info(
            f"Pre-Train AutoEncoder Epoch [{epoch+1}/{epochs}] | "
            f"Recons: {mean_recon:.4f} | Total: {mean_total:.4f}"
        )
    return best_autoencoders


# Classifer Train
def pretrain_classifier(
    autoencoders,
    classifier,
    dataloader,
    epochs,
    lambda_organ,
    lambda_mic,
    pos_weight,
    learning_rate,
    log_path,
    device
):
    autoencoders = autoencoders.to(device)
    classifier = classifier.to(device)
    autoencoders.eval()      # freeze AE
    classifier.train()
    optimizer = optim.Adam(classifier.parameters(), lr=learning_rate)
    logger = setup_logger(log_path)
    # one BCE per label (important for pos_weight)
    criterions = [
        nn.BCEWithLogitsLoss(pos_weight=pos_weight[i].to(device))
        for i in range(len(pos_weight))
    ]
    best_loss = float("inf")
    best_classifier = copy.deepcopy(classifier)
    for epoch in range(epochs):
        epoch_losses = []
        all_preds, all_labels = [], []
        for batch_x, batch_y in dataloader:
            batch_x = batch_x.to(device)
            batch_y = batch_y.to(device)
            # ---- AE frozen ----
            with torch.no_grad():
                _, latent = autoencoders(batch_x)
            # ---- classifier ----
            optimizer.zero_grad()
            organ_logits, mic_logits = classifier(latent)
            if organ_logits is not None:
                logits = torch.cat([organ_logits, mic_logits], dim=1)
                # ---- loss ----
    #            pos_weight = []
    #            for i in range(batch_y.shape[1]):
    #                n_MIC = batch_y[:, i].sum()
    #                n_nonMIC = batch_y.shape[0] - n_MIC
    #                pos_weight.append(n_nonMIC / (n_MIC + 1e-6))
    #            criterions = [
    #                nn.BCEWithLogitsLoss(pos_weight=pos_weight[i].to(device))
    #                for i in range(len(pos_weight))
    #            ]
                losses_per_label = [
                    criterions[i](logits[:, i], batch_y[:, i])
                    for i in range(batch_y.shape[1])
                ]
                organ_loss = torch.stack(losses_per_label[:-1]).sum()
                mic_loss = losses_per_label[-1]
                total_loss = lambda_organ * organ_loss + lambda_mic * mic_loss
            else:
                logits = mic_logits
                total_loss = criterions[0](logits[:, 0], batch_y[:, 0])    
            total_loss.backward()
            optimizer.step()
            epoch_losses.append(total_loss.item())
            all_preds.append(torch.sigmoid(logits).detach().cpu())
            all_labels.append(batch_y.detach().cpu())
        # ---- epoch summary ----
        mean_loss = np.mean(epoch_losses)
        if mean_loss < best_loss:
            best_loss = mean_loss
            best_classifier = copy.deepcopy(classifier)
        # ---- metrics ----
        all_preds = torch.cat(all_preds).numpy()
        bin_preds = (all_preds > 0.5).astype(int)
        all_labels = torch.cat(all_labels).numpy()
        n_labels = all_labels.shape[1]
        acc, precision, recall, f1, auc = [np.zeros(n_labels) for _ in range(5)]
        for i in range(n_labels):
            acc[i] = accuracy_score(all_labels[:, i], bin_preds[:, i])
            precision[i] = precision_score(all_labels[:, i], bin_preds[:, i], average="macro", zero_division=0)
            recall[i] = recall_score(all_labels[:, i], bin_preds[:, i], average="macro", zero_division=0)
            f1[i] = f1_score(all_labels[:, i], bin_preds[:, i], average="macro", zero_division=0)
            try:
                auc[i] = roc_auc_score(all_labels[:, i], all_preds[:, i])
            except ValueError:
                auc[i] = np.nan
        logger.info(
            f"Pre-Train Classifier Epoch [{epoch+1}/{epochs}] | "
            f"Loss: {mean_loss:.4f} | "
            f"Acc: {acc} | precision: {precision} | recall: {recall} | F1: {f1} | AUC: {auc}"
        )
    return best_classifier

# Joint Train
def joint_train(
    autoencoders,
    classifier,
    dataloader,
    epochs,
    lambda_recons,
    lambda_phenotype,
    lambda_organ,
    lambda_mic,
    pos_weight,
    learning_rate,
    log_path,
    device
):
    autoencoders = autoencoders.to(device)
    classifier = classifier.to(device)
    autoencoders.train()
    classifier.train()
    all_params = list(autoencoders.parameters()) + list(classifier.parameters())
    optimizer = torch.optim.Adam(all_params, lr=learning_rate)
    logger = setup_logger(log_path)
    recon_criterion = nn.MSELoss()
    criterions = [
        nn.BCEWithLogitsLoss(pos_weight=pos_weight[i].to(device))
        for i in range(len(pos_weight))
    ]
    best_loss = float("inf")
    best_model = {
        "autoencoders": copy.deepcopy(autoencoders),
        "classifier": copy.deepcopy(classifier)
    }
    for epoch in range(epochs):
        epoch_losses = []
        all_preds, all_labels = [], []
        for batch_x, batch_y in dataloader:
            batch_x, batch_y = batch_x.to(device), batch_y.to(device)
            optimizer.zero_grad()
            # ---- forward AE ----
            reconstructed, latent = autoencoders(batch_x)
            reconstruction_loss = recon_criterion(reconstructed, batch_x)
            # ---- forward classifier ----
            organ_logits, mic_logits = classifier(latent)
            if organ_logits is not None:
                logits = torch.cat([organ_logits, mic_logits], dim=1)
                # ---- multi-label BCE ----
    #            pos_weight = []
    #            for i in range(batch_y.shape[1]):
    #                n_MIC = batch_y[:, i].sum()
    #                n_nonMIC = batch_y.shape[0] - n_MIC
    #                pos_weight.append(n_nonMIC / (n_MIC + 1e-6))
    #            criterions = [
    #                nn.BCEWithLogitsLoss(pos_weight=pos_weight[i].to(device))
    #                for i in range(len(pos_weight))
    #            ]
                losses_per_label = [
                    criterions[i](logits[:, i], batch_y[:, i])
                    for i in range(batch_y.shape[1])
                ]
                organ_loss = torch.stack(losses_per_label[:-1]).sum()
                mic_loss = losses_per_label[-1]
                total_loss = lambda_recons * reconstruction_loss + \
                             lambda_phenotype * (lambda_organ * organ_loss + lambda_mic * mic_loss)
            else:
                logits = mic_logits
                phenotype_loss = criterions[0](logits[:, 0], batch_y[:, 0])
                total_loss = lambda_recons * reconstruction_loss + lambda_phenotype * phenotype_loss
            total_loss.backward()
            optimizer.step()
            epoch_losses.append(total_loss.item())
            all_preds.append(torch.sigmoid(logits).detach().cpu())
            all_labels.append(batch_y.detach().cpu())
        # ---- epoch metrics ----
        mean_loss = np.mean(epoch_losses)
        if mean_loss < best_loss:
            best_loss = mean_loss
            best_model = {
                "autoencoders": copy.deepcopy(autoencoders),
                "classifier": copy.deepcopy(classifier)
            }
        all_preds = torch.cat(all_preds).numpy()
        bin_preds = (all_preds > 0.5).astype(int)
        all_labels = torch.cat(all_labels).numpy()
        n_labels = all_labels.shape[1]
        acc, precision, recall, f1, auc = [np.zeros(n_labels) for _ in range(5)]
        for i in range(n_labels):
            acc[i] = accuracy_score(all_labels[:, i], bin_preds[:, i])
            precision[i] = precision_score(all_labels[:, i], bin_preds[:, i], average="macro", zero_division=0)
            recall[i] = recall_score(all_labels[:, i], bin_preds[:, i], average="macro", zero_division=0)
            f1[i] = f1_score(all_labels[:, i], bin_preds[:, i], average="macro", zero_division=0)
            try:
                auc[i] = roc_auc_score(all_labels[:, i], all_preds[:, i])
            except ValueError:
                auc[i] = np.nan
        logger.info(
            f"joint_train Epoch [{epoch+1}/{epochs}] | "
            f"Loss: {mean_loss:.4f} | "
            f"Acc: {acc} | precision: {precision} | recall: {recall} | F1: {f1} | AUC: {auc}"
        )
    return best_model


# ---------------------------- Main Training Entry -------------------------
def GeneProgram(
    scRNA_file,
    phenotype_file,
    latent_dim,
    epochs,
    Encoder_hidden_dims,
    Classifier_hidden_dims,
    learning_rate,
    lambda_recons,
    lambda_phenotype,
    lambda_organ,
    lambda_mic,
    output_path,
    batch_size,
    seed,
    organ_flag
):
    set_seed(seed)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    expr = pd.read_csv(scRNA_file, index_col=0).values
    phenos = pd.read_csv(phenotype_file, index_col=0)
    category_mappings = {}
    if organ_flag:
        for col in phenos.columns:
            phenos[col] = phenos[col].apply(
                lambda x: 1 if "MIC" in str(x) and "non" not in str(x) else 0
            )
            category_mappings[col] = {'MIC':1,'non-MIC':0}
    else:
        mic_col = phenos.columns[-1]
        phenos = phenos[[mic_col]]
        phenos[mic_col] = phenos[mic_col].apply(
            lambda x: 1 if "MIC" in str(x) and "non" not in str(x) else 0
        )
        category_mappings[mic_col] = {'MIC':1,'non-MIC':0}

    dataset = GroupedDataset(expr, phenos)
    dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True, num_workers=0,generator=torch.Generator().manual_seed(seed))
    pos_weight = compute_pos_weight(torch.tensor(dataset.labels))

    input_dim = expr.shape[1]
    autoencoders = AutoEncoder(input_dim, latent_dim, Encoder_hidden_dims).to(device)
    phenotype_dim = phenos.shape[1]
    classifier = PhenotypeClassifier(latent_dim, phenotype_dim, Classifier_hidden_dims, organ_flag).to(device)

    os.makedirs(output_path, exist_ok=True)
    log_path = os.path.join(output_path, "train.log")

    best_autoencoders = pretrain_autoencoders(autoencoders, dataloader, epochs, lambda_recons, learning_rate, log_path, device)
    best_classifier = pretrain_classifier(best_autoencoders, classifier, dataloader, epochs, lambda_organ, lambda_mic, pos_weight, learning_rate, log_path, device)
    best_model = joint_train(best_autoencoders, best_classifier, dataloader, epochs, lambda_recons, lambda_phenotype, lambda_organ, lambda_mic, pos_weight, learning_rate, log_path, device)

    save_dict = {"autoencoders": best_model["autoencoders"].state_dict(),
                 "classifier": best_model["classifier"].state_dict()}
    torch.save(save_dict, os.path.join(output_path, "full_model.pt"))
    print(f"Model saved to {os.path.join(output_path, 'full_model.pt')}")

# ---------------------------- Final Inference -----------------------------
def FinalInference(
    autoencoders,
    classifier,
    dataset,
    device,
    category_mappings,
    figure_path,
    percentile_ratio,
    std_threshold
):
    autoencoders = autoencoders.to(device).eval()
    classifier = classifier.to(device).eval()

    expr_tensor = torch.tensor(dataset.expr, dtype=torch.float32).to(device)
    with torch.no_grad():
        _, latent = autoencoders(expr_tensor)
        organ_logits, mic_logits = classifier(latent)
        if organ_logits is not None:
            logits = torch.cat([organ_logits, mic_logits], dim=1)
        else:
            logits = mic_logits

    all_preds = torch.sigmoid(logits).detach().cpu().numpy()
    all_labels = dataset.labels
    bin_preds = (all_preds > 0.5).astype(int)
    n_labels = all_labels.shape[1]

    # Metrics
    acc, precision, recall, f1, auc = [np.zeros(n_labels) for _ in range(5)]
    for i in range(n_labels):
        acc[i] = accuracy_score(all_labels[:,i], bin_preds[:,i])
        precision[i] = precision_score(all_labels[:,i], bin_preds[:,i], average="macro", zero_division=0)
        recall[i] = recall_score(all_labels[:,i], bin_preds[:,i], average="macro", zero_division=0)
        f1[i] = f1_score(all_labels[:,i], bin_preds[:,i], average="macro", zero_division=0)
        try:
            auc[i] = roc_auc_score(all_labels[:,i], all_preds[:,i])
        except:
            auc[i] = np.nan
    print(f"Acc: {acc}\nPrecision: {precision}\nRecall: {recall}\nF1: {f1}\nAUC: {auc}")

    # ROC Curve
    pheno_columns = list(category_mappings.keys())
    plt.figure(figsize=(8,6))
    for i in range(n_labels):
        try:
            fpr, tpr, _ = roc_curve(all_labels[:,i], all_preds[:,i])
            plt.plot(fpr, tpr, label=f"{pheno_columns[i]} (AUC={auc[i]:.2f})")
        except:
            continue
    plt.plot([0,1],[0,1],'k--'); plt.xlabel("FPR"); plt.ylabel("TPR"); plt.title("ROC Curve"); plt.legend()
    os.makedirs(figure_path, exist_ok=True)
    plt.savefig(os.path.join(figure_path,"ROC.pdf"), dpi=300)

    # Latent & labels sorting
    sort_keys = [all_labels[:,i] for i in range(len(pheno_columns))]
    order = np.lexsort(sort_keys[::-1])
    sorted_latent = torch.sigmoid(latent[order]).cpu().numpy()
    sorted_labels = all_labels[order]

    pd.DataFrame(sorted_latent).to_csv(os.path.join(figure_path,"sorted_latent.csv"), index=False)
    pd.DataFrame(sorted_labels).to_csv(os.path.join(figure_path,"sorted_labels.csv"), index=False)
    # F matrix
    F_all = estimate_F(latent, expr_tensor)
    pd.DataFrame(F_all.cpu().numpy()).to_csv(os.path.join(figure_path,"F_all.csv"), index=False)

    # Heatmap plotting 
    color_palettes = ["tab20", "Set3", "Pastel1", "Dark2", "hls", "husl"]
    legend_elements = []
    phenotype_colors = {}
    num=0
    for i, col in enumerate(pheno_columns):
        mapping = category_mappings[col]
        palette = sns.color_palette(color_palettes[i % len(color_palettes)], len(mapping))
        col_color_map = {}
        for name, code in mapping.items():
            col_color_map[code] = palette[code]
            legend_elements.append(Patch(facecolor=palette[code], label=f"{col}: {name}"))
            num=num+1
        phenotype_colors[col] = col_color_map
    label_color_rgb = np.zeros((*sorted_labels.shape, 3))
    for i, col in enumerate(pheno_columns):
        for j in range(sorted_labels.shape[0]):
            code = sorted_labels[j, i]
            label_color_rgb[j, i] = phenotype_colors[col][code]

    fig, (ax_label, ax_latent) = plt.subplots(
        ncols=2, figsize=(12, 6), gridspec_kw={'width_ratios': [len(pheno_columns), sorted_latent.shape[1]]}, sharey=True
    )
    fig = plt.figure(figsize=(10, 6))
    gs = fig.add_gridspec(1, 2, width_ratios=[len(pheno_columns) / 2, sorted_latent.shape[1] / 2])
    ax_label = fig.add_subplot(gs[0, 0])   
    ax_label.imshow(label_color_rgb, aspect='auto')
    ax_label.set_xticks(range(len(pheno_columns)))
    ax_label.set_xticklabels(pheno_columns, rotation=90)    
    # Normalize latent for color mapping
    norm = mpl.colors.Normalize(
        vmin=np.percentile(sorted_latent, percentile_ratio),
        vmax=np.percentile(sorted_latent, 100 - percentile_ratio)
    )
    latent_pd = pd.DataFrame(sorted_latent)
    stds = latent_pd.std()
    cols_to_keep = stds[stds > std_threshold].index
    latent_pd_keep = latent_pd[cols_to_keep]
    selected_sorted_latent = latent_pd_keep.values
    ax_latent = fig.add_subplot(gs[0, 1])
    sns.heatmap(selected_sorted_latent, ax=ax_latent, cmap='viridis', cbar=False, norm=norm)
    ax_latent.set_xlabel("Latent Dimension")
    ax_latent.set_xticklabels(cols_to_keep) 
    ax_latent.set_yticks([])
    ax_latent.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left')
    sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
    sm.set_array([])  # Required for ScalarMappable
    cbar_ax = fig.add_axes([0.85, 0.1, 0.02, 1-0.12*num])  # [left, bottom, width, height]
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label("Latent Value")
    plt.tight_layout()
    plt.savefig(os.path.join(figure_path, "Latent_Heatmap.pdf"), dpi=300)
    print(f"Latent heatmap saved to: {os.path.join(figure_path, 'Latent_Heatmap.pdf')}")


# ---------------------------- Output figures --------------------------------
def GeneProgramOutput(
    scRNA_file,
    phenotype_file,
    model_path,
    output_path,
    latent_dim,
    Encoder_hidden_dims,
    Classifier_hidden_dims,
    seed,
    percentile_ratio,
    std_threshold,
    organ_flag
):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    # -------------------- Load Data --------------------
    expr = pd.read_csv(scRNA_file, index_col=0).values
    phenos = pd.read_csv(phenotype_file, index_col=0)
    # Binarize MIC labels
    category_mappings = {}
    if organ_flag:
        for col in phenos.columns:
            phenos[col] = phenos[col].apply(
                lambda x: 1 if "MIC" in str(x) and "non" not in str(x) else 0
            )
            category_mappings[col] = {'MIC':1,'non-MIC':0}
    else:
        mic_col = phenos.columns[-1]
        phenos = phenos[[mic_col]]
        phenos[mic_col] = phenos[mic_col].apply(
            lambda x: 1 if "MIC" in str(x) and "non" not in str(x) else 0
        )
        category_mappings[mic_col] = {'MIC':1,'non-MIC':0}
    # -------------------- Dataset & Dataloader --------------------
    dataset = GroupedDataset(expr, phenos)
    # ------------------------- Load model -------------------------
    input_dim = expr.shape[1]
    phenotype_dim = phenos.shape[1]
    checkpoint = torch.load(os.path.join(model_path, "full_model.pt"), map_location=device)
    autoencoders =  AutoEncoder(input_dim, latent_dim, Encoder_hidden_dims).to(device)
    autoencoders.load_state_dict(checkpoint["autoencoders"])
    autoencoders.eval()
    classifier = PhenotypeClassifier(latent_dim, phenotype_dim, Classifier_hidden_dims,organ_flag).to(device)
    classifier.load_state_dict(checkpoint["classifier"])
    classifier.eval()
    # ------------------------- Inference -------------------------
    FinalInference(autoencoders, classifier, dataset,  device, category_mappings, output_path, percentile_ratio, std_threshold)


# ---------------------------- CLI Entry -----------------------------------
def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--Flag', type=str, required=True, choices=['Training','Inferencing'])
    parser.add_argument('--scRNA_file', type=str, required=True)
    parser.add_argument('--phenotype_file', type=str, required=True)
    parser.add_argument('--model_path', type=str, default=None)
    parser.add_argument('--output_path', type=str, required=True)
    parser.add_argument('--organ_flag', action='store_true',help='Whether to enable organ classifier (default: False)')
    parser.add_argument('--latent_dim', type=int, default=32)
    parser.add_argument('--epochs', type=int, default=200)
    parser.add_argument('--Encoder_hidden_dims', type=int, nargs='+', default=[512,256,128,64])
    parser.add_argument('--Classifier_hidden_dims', type=int, nargs='+', default=[16,4])
    parser.add_argument('--learning_rate', type=float, default=1e-3)
    parser.add_argument('--lambda_recons', type=float, default=1.0)
    parser.add_argument('--lambda_phenotype', type=float, default=1.0)
    parser.add_argument('--lambda_organ', type=float, default=1.0)
    parser.add_argument('--lambda_mic', type=float, default=1.0)
    parser.add_argument('--batch_size', type=int, default=64)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--percentile_ratio', type=float, default=5)
    parser.add_argument('--std_threshold', type=float, default=0.1)
    args = parser.parse_args()

    if args.Flag=="Training":
        GeneProgram(args.scRNA_file, args.phenotype_file, args.latent_dim, args.epochs,
                    args.Encoder_hidden_dims, args.Classifier_hidden_dims, args.learning_rate,
                    args.lambda_recons, args.lambda_phenotype, args.lambda_organ, args.lambda_mic,
                    args.output_path, args.batch_size, args.seed,args.organ_flag)
    elif args.Flag=="Inferencing":
        GeneProgramOutput(args.scRNA_file, args.phenotype_file, args.model_path, args.output_path,
                          args.latent_dim, args.Encoder_hidden_dims, args.Classifier_hidden_dims,
                          args.seed, args.percentile_ratio, args.std_threshold,args.organ_flag)

if __name__=="__main__":
    main()
