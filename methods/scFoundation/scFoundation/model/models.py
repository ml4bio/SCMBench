import argparse
import pandas as pd
import numpy as np
from sklearn.linear_model import Lasso
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LassoCV
from sklearn.model_selection import KFold
import json
import torch
from torch import nn
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import train_test_split
import anndata as ad
import matplotlib.pyplot as plt
from torch.optim.lr_scheduler import ReduceLROnPlateau
import sys
sys.path.append("../model/") # path to this folder
from load import *

class MLP(nn.Module):
    def __init__(self, input_dim, hidden_dims=[128, 64]):
        super().__init__()
        layers = []
        prev_dim = input_dim
        
        for hidden_dim in hidden_dims:
            layers.extend([
                nn.Linear(prev_dim, hidden_dim),
                nn.BatchNorm1d(hidden_dim),
                nn.ReLU(),
                nn.Dropout(0.3)
            ])
            prev_dim = hidden_dim
            
        # 添加一个线性层作为主要预测层
        layers.append(nn.Linear(prev_dim, 64))
        layers.append(nn.ReLU())
        layers.append(nn.Linear(64, 1))
        self.mlp = nn.Sequential(*layers)
        
    def forward(self, x):
        return self.mlp(x)

class BaseAgePredictor(nn.Module):
    def __init__(self, ckpt_path):
        super().__init__()
        self.ckpt_path = ckpt_path
        
    def build(self):
        model, model_config = load_model_frommmf(self.ckpt_path)
        self.token_emb = model.token_emb
        self.pos_emb = model.pos_emb
        self.encoder = model.encoder
        self.model_config = model_config
        
    def freeze_pretrained(self):
        for param in self.token_emb.parameters():
            param.requires_grad = False
        for param in self.pos_emb.parameters():
            param.requires_grad = False
        for param in self.encoder.parameters():
            param.requires_grad = False
            
    def unfreeze_last_encoder_layer(self):
        for param in self.encoder.transformer_encoder[-1].parameters():
            param.requires_grad = True
            
    def get_embeddings(self, sample_list):
        x = sample_list['x']
        value_labels = x > 0
        x, x_padding = gatherData(x, value_labels, self.model_config['pad_token_id'])
        data_gene_ids = torch.arange(19264, device=x.device).repeat(x.shape[0], 1)
        position_gene_ids, _ = gatherData(data_gene_ids, value_labels,
                                        self.model_config['pad_token_id'])
        
        x = self.token_emb(torch.unsqueeze(x, 2).float(), output_weight=0)
        position_emb = self.pos_emb(position_gene_ids)
        x += position_emb
        logits = self.encoder(x, x_padding)
        
        max_pool = torch.max(logits, dim=1)[0]  
        mean_pool = torch.mean(logits, dim=1)
        return torch.cat([max_pool, mean_pool], dim=1)

class AgePredictorLasso(BaseAgePredictor):
    def __init__(self, ckpt_path):
        super().__init__(ckpt_path)
        self.lasso = None
        self.scaler = StandardScaler()
        
    def train_lasso_cv(self, train_loader, val_loader, n_alphas=100, cv=5):
        print("Collecting training data...")
        embeddings_list = []
        ages_list = []
        self.eval()
        
        with torch.no_grad():
            for batch in train_loader:
                x = batch['x'].to(next(self.parameters()).device)
                age = batch['age'].numpy()
                emb = self.get_embeddings({'x': x}).cpu().numpy()
                embeddings_list.append(emb)
                ages_list.append(age)
        
        X = np.vstack(embeddings_list)
        y = np.vstack(ages_list).ravel()
        X_scaled = self.scaler.fit_transform(X)
        
        alphas = np.logspace(-4, 2, n_alphas)
        self.lasso = LassoCV(alphas=alphas, cv=cv, n_jobs=-1, max_iter=2000, verbose=True)
        self.lasso.fit(X_scaled, y)
        
        return self.lasso.alpha_
    
    def forward(self, sample_list):
        emb = self.get_embeddings(sample_list)
        
        if self.training:
            emb_scaled = self.scaler.transform(emb.detach().cpu().numpy())
            pred = self.lasso.predict(emb_scaled)
            pred_tensor = torch.FloatTensor(pred).view(-1, 1).to(emb.device)
            return pred_tensor + (emb.mean(dim=1, keepdim=True) * 0.0)
        else:
            emb_scaled = self.scaler.transform(emb.detach().cpu().numpy())
            pred = self.lasso.predict(emb_scaled)
            return torch.FloatTensor(pred).view(-1, 1).to(emb.device)

class AgePredictorMLP(BaseAgePredictor):
    def __init__(self, ckpt_path, hidden_dims=[512, 256, 128]):
        super().__init__(ckpt_path)
        self.hidden_dims = hidden_dims
        self.mlp = None
        
    def build(self):
        super().build()
        embedding_dim = self.encoder.transformer_encoder[0].self_attn.in_proj_weight.size(1) * 2
        self.mlp = MLP(embedding_dim, self.hidden_dims)
        
    def freeze_mlp(self):
        """冻结MLP层"""
        for param in self.mlp.parameters():
            param.requires_grad = False
            
    def unfreeze_mlp(self):
        """解冻MLP层"""
        for param in self.mlp.parameters():
            param.requires_grad = True
            
    def forward(self, sample_list):
        emb = self.get_embeddings(sample_list)
        return self.mlp(emb)