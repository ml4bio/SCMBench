import argparse
import pandas as pd
import numpy as np
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

class SingleCellDataset(Dataset):
    def __init__(self, data_path, age_labels, split_type='train'):
        # Load split data
        # self.data = pd.read_csv(data_path)
        self.data = pd.read_csv(data_path).values
        # Get corresponding ages from json
        self.ages = np.array(age_labels[split_type])
        
        # import pdb;pdb.set_trace()
        
    def __len__(self):
        return len(self.data)
        
    def __getitem__(self, idx):
        

        return {
            'x': torch.FloatTensor(self.data[idx].astype(np.float32)),
            'age': torch.FloatTensor([self.ages[idx]])
        }
    
class AgePredictor(nn.Module):
    def __init__(self, ckpt_path, frozenmore=True):
        super().__init__()
        self.ckpt_path = ckpt_path
        self.frozenmore = frozenmore
        
    def build(self):
        model, model_config = load_model_frommmf(self.ckpt_path)
        self.token_emb = model.token_emb
        self.pos_emb = model.pos_emb
        self.encoder = model.encoder
        
        if self.frozenmore:
            for _, p in self.token_emb.named_parameters():
                p.requires_grad = False
            for _, p in self.pos_emb.named_parameters():
                p.requires_grad = False
            print('self.pos_emb and self.token_emb frozen')
        
        for na, param in self.encoder.named_parameters():
            param.requires_grad = False
        for na, param in self.encoder.transformer_encoder[-1].named_parameters():
            print('self.encoder.transformer_encoder', na, 'have grad')
            param.requires_grad = True
            
        hidden_dim = model_config['encoder']['hidden_dim']
        self.fc1 = nn.Sequential(
            nn.LayerNorm(hidden_dim),
            nn.Linear(hidden_dim, 256),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(256, 128),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(128, 1)
        )
        self.model_config = model_config
        
    def forward(self, sample_list):
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
        
        logits, _ = torch.max(logits, dim=1)
        age_pred = self.fc1(logits)
        return age_pred

def plot_predictions(predictions, targets, title, lr=0.0001):
    plt.figure(figsize=(8, 6))
    plt.scatter(targets, predictions, alpha=0.5)
    plt.plot([min(targets), max(targets)], [min(targets), max(targets)], 'r--')
    plt.xlabel('True Age')
    plt.ylabel('Predicted Age')
    plt.title(f'{title}\nPearson r: {np.corrcoef(predictions, targets)[0,1]:.3f}')
    plt.savefig(f'{title.lower().replace(" ", "_")}_scatter_{lr}.png')
    plt.close()

def evaluate_model(model, data_loader, device='cuda', split_name=None, lr=0.0001):
    model.eval()
    predictions = []
    targets = []
    
    with torch.no_grad():
        for batch in data_loader:
            x = batch['x'].to(device)
            age = batch['age'].to(device)
            pred = model({'x': x})
            
            predictions.extend(pred.cpu().numpy().flatten())
            targets.extend(age.cpu().numpy().flatten())
    
    predictions = np.array(predictions)
    targets = np.array(targets)
    
    if split_name:
        plot_predictions(predictions, targets, f'{split_name} Predictions', lr)
    
    mae = np.mean(np.abs(predictions - targets))
    r2 = np.corrcoef(predictions, targets)[0,1] ** 2
    pearson_r = np.corrcoef(predictions, targets)[0,1]
    
    return {
        'mae': float(mae),
        'r2': float(r2),
        'pearson_r': float(pearson_r)
    }

def train_model(model, train_loader, val_loader, test_loader, epochs=20, lr=0.0001, device='cuda'):
    optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-5)
    scheduler = ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=5, verbose=True)
    criterion = nn.L1Loss()
    best_val_loss = float('inf')
    
    for epoch in range(epochs):
        # Training
        model.train()
        train_loss = 0
        for batch in train_loader:
            x = batch['x'].to(device)
            age = batch['age'].to(device)
            
            pred = model({'x': x})
            loss = criterion(pred, age)
            
            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()
            
            train_loss += loss.item()
        
        # Evaluation
        if (epoch + 1) % 1 == 0:
            results = {
                'train': evaluate_model(model, train_loader),
                'val': evaluate_model(model, val_loader)
            }
            
            print(f'Epoch {epoch+1}:')
            print(f'Train - MAE: {results["train"]["mae"]:.4f}, Pearson_r: {results["train"]["pearson_r"]:.4f}')
            print(f'Val - MAE: {results["val"]["mae"]:.4f}, Pearson_r: {results["val"]["pearson_r"]:.4f}')
            
            scheduler.step(results['val']['mae'])
            
            if results['val']['mae'] < best_val_loss:
                best_val_loss = results['val']['mae']
                torch.save({
                    'epoch': epoch,
                    'model_state_dict': model.state_dict(),
                    'optimizer_state_dict': optimizer.state_dict(),
                    'val_loss': best_val_loss,
                }, f'best_model_lr_{lr}.pth')
                with open(f'results_lr_{lr}.json', 'w') as f:
                    json.dump(results, f, indent=4)
    
    return best_val_loss


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--lr', type=float, default=0.0001, help='Learning rate')
    args = parser.parse_args()
    
    # Load age labels from json
    with open('/mnt/nas/user/yixuan/Bio_clock/cxg-clocks-v2/blood/age_splitted.json', 'r') as f:
        age_labels = json.load(f)
    
    # Create datasets
    print("Loading data...")
    print(f"Number of age labels - Train: {len(age_labels['train'])}, Val: {len(age_labels['val'])}, Test: {len(age_labels['test'])}")
    
    processed_data_path = '/mnt/nas/user/yixuan/Bio_clock/cxg-clocks-v2/blood/processed/'

    # Create datasets
    train_dataset = SingleCellDataset(processed_data_path+'B_cell_train_aligned_top2000_var.csv', age_labels, 'train')
    val_dataset = SingleCellDataset(processed_data_path+'B_cell_val_aligned_top2000_var.csv', age_labels, 'val')
    test_dataset = SingleCellDataset(processed_data_path+'B_cell_test_aligned_top2000_var.csv', age_labels, 'test')
    
    print(f"Dataset sizes - Train: {len(train_dataset)}, Val: {len(val_dataset)}, Test: {len(test_dataset)}")
    print(f"Input feature dimension: {train_dataset.data.shape[1]}")

    train_loader = DataLoader(train_dataset, batch_size=20, shuffle=True, drop_last=True)
    val_loader = DataLoader(val_dataset, batch_size=20, drop_last=True)
    test_loader = DataLoader(test_dataset, batch_size=20, drop_last=True)
    
    # Initialize and build model
    model = AgePredictor(ckpt_path='./models/models.ckpt')
    model.build()
    model = model.cuda()
    
    # Train model
    train_model(model, train_loader, val_loader, test_loader, lr=args.lr)
    
    # Final evaluation
    results = {
        'train': evaluate_model(model, train_loader,split_name='train', lr=args.lr),
        'val': evaluate_model(model, val_loader,split_name='val', lr=args.lr),
        'test': evaluate_model(model, test_loader,split_name='test', lr=args.lr)
    }
    
    # Save final results
    with open(f'final_results_lr_{args.lr}.json', 'w') as f:
        json.dump(results, f, indent=4)

   