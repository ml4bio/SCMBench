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

from datetime import datetime

class SingleCellDataset(Dataset):
    def __init__(self, data_path, age_labels, split_type='train'):
        self.data = pd.read_csv(data_path).values
        self.ages = np.array(age_labels[split_type])
        
    def __len__(self):
        return len(self.data)
        
    def __getitem__(self, idx):
        return {
            'x': torch.FloatTensor(self.data[idx].astype(np.float32)),
            'age': torch.FloatTensor([self.ages[idx]])
        }
    
class AgePredictorTwoStage(nn.Module):
    def __init__(self, ckpt_path):
        super().__init__()
        self.ckpt_path = ckpt_path
        self.lasso = None
        self.scaler = StandardScaler()
        
    def build(self):
        model, model_config = load_model_frommmf(self.ckpt_path)
        self.token_emb = model.token_emb
        self.pos_emb = model.pos_emb
        self.encoder = model.encoder
        self.model_config = model_config
        
    def freeze_pretrained(self):
        """冻结所有预训练参数"""
        for param in self.token_emb.parameters():
            param.requires_grad = False
        for param in self.pos_emb.parameters():
            param.requires_grad = False
        for param in self.encoder.parameters():
            param.requires_grad = False
            
    def unfreeze_last_encoder_layer(self):
        """解冻最后一个encoder层"""
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
        
        # 使用多种池化方式
        max_pool = torch.max(logits, dim=1)[0]  
        mean_pool = torch.mean(logits, dim=1)
        return torch.cat([max_pool, mean_pool], dim=1)
        
    
    def train_lasso_cv(self, train_loader, val_loader, n_alphas=100, cv=5):
        """训练LassoCV找到最优alpha"""
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
        
        # 标准化特征
        X_scaled = self.scaler.fit_transform(X)
        
        # 设置alphas范围
        alphas = np.logspace(-4, 2, n_alphas)
        
        # 训练LassoCV
        print("Training LassoCV...")
        self.lasso = LassoCV(
            alphas=alphas,
            cv=cv,
            n_jobs=-1,
            max_iter=2000,
            verbose=True
        )
        self.lasso.fit(X_scaled, y)
        
        print(f"Best alpha: {self.lasso.alpha_:.4f}")
        
        # 评估验证集
        val_embeddings = []
        val_ages = []
        with torch.no_grad():
            for batch in val_loader:
                x = batch['x'].to(next(self.parameters()).device)
                age = batch['age'].numpy()
                emb = self.get_embeddings({'x': x}).cpu().numpy()
                val_embeddings.append(emb)
                val_ages.append(age)
        
        X_val = np.vstack(val_embeddings)
        y_val = np.vstack(val_ages).ravel()
        X_val_scaled = self.scaler.transform(X_val)
        
        val_pred = self.lasso.predict(X_val_scaled)
        val_mae = np.mean(np.abs(val_pred - y_val))
        val_r = np.corrcoef(val_pred, y_val)[0,1]
        
        print(f"Validation MAE: {val_mae:.4f}")
        print(f"Validation Pearson r: {val_r:.4f}")
        
        return self.lasso.alpha_
    
    def forward(self, sample_list):
        """支持两阶段训练的前向传播"""
        emb = self.get_embeddings(sample_list)
        
        if self.training:  # 在训练模式下，使用原始embeddings以保持梯度
            emb_scaled = self.scaler.transform(emb.detach().cpu().numpy())
            pred = self.lasso.predict(emb_scaled)
            pred_tensor = torch.FloatTensor(pred).view(-1, 1).to(emb.device)
            
            # 添加一个小的线性层来保持梯度流
            final_pred = pred_tensor + (emb.mean(dim=1, keepdim=True) * 0.0)
            return final_pred
        else:  # 在评估模式下，直接使用LASSO预测
            emb_scaled = self.scaler.transform(emb.detach().cpu().numpy())
            pred = self.lasso.predict(emb_scaled)
            return torch.FloatTensor(pred).view(-1, 1).to(emb.device)

def evaluate_model(model, data_loader, device='cuda', split_name=None):
    """评估模型性能"""
    model.eval()
    predictions = []
    targets = []
    
    with torch.no_grad():
        for batch in data_loader:
            x = batch['x'].to(device)
            age = batch['age'].cpu().numpy()
            pred = model({'x': x}).detach().cpu().numpy()
            predictions.extend(pred.flatten())
            targets.extend(age.flatten())
    
    predictions = np.array(predictions)
    targets = np.array(targets)
    
    if split_name:
        plt.figure(figsize=(8, 6))
        plt.scatter(targets, predictions, alpha=0.5)
        plt.plot([min(targets), max(targets)], [min(targets), max(targets)], 'r--')
        plt.xlabel('True Age')
        plt.ylabel('Predicted Age')
        plt.title(f'{split_name}\nPearson r: {np.corrcoef(predictions, targets)[0,1]:.3f}')
        plt.savefig(f'{split_name.lower().replace(" ", "_")}_scatter.png')
        plt.close()
    
    mae = np.mean(np.abs(predictions - targets))
    r2 = np.corrcoef(predictions, targets)[0,1] ** 2
    pearson_r = np.corrcoef(predictions, targets)[0,1]
    
    return {
        'mae': float(mae),
        'r2': float(r2),
        'pearson_r': float(pearson_r)
    }

class CombinedLoss(nn.Module):
   def __init__(self, mae_weight=0.5, mse_weight=0.5):
       super().__init__()
       self.mae_loss = nn.L1Loss()
       self.mse_loss = nn.MSELoss()
       self.mae_weight = mae_weight
       self.mse_weight = mse_weight
       
   def forward(self, pred, target):
       mae = self.mae_loss(pred, target)
       mse = self.mse_loss(pred, target) 
       return self.mae_weight * mae + self.mse_weight * mse

def train_stage2(model, train_loader, val_loader, test_loader, epochs, lr, device='cuda'):
   optimizer = torch.optim.AdamW(
       filter(lambda p: p.requires_grad, model.parameters()),
       lr=lr, 
       betas=[0.9, 0.95],
       weight_decay=0.01
   )
   
#    scheduler = torch.optim.CosineAnnealingLR(
#         optimizer,
#         T_max=epochs,
#         eta_min=1e-7
#     )
   
   criterion = CombinedLoss(mae_weight=0.7, mse_weight=0.3)
   best_val_loss = float('inf')
   
   losses = []
   print("\nStarting stage 2 training...")
   
   for epoch in range(epochs):
       model.train()
       epoch_loss = []
       for batch in train_loader:
           x = batch['x'].to(device)
           age = batch['age'].to(device)
           
           pred = model({'x': x})
           loss = criterion(pred, age)
           epoch_loss.append(loss.item())
           
           optimizer.zero_grad()
           loss.backward()
           torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
           optimizer.step()
        #    scheduler.step()
           
       avg_loss = np.mean(epoch_loss)
       losses.append(avg_loss)
       
       results = {
           'train': evaluate_model(model, train_loader),
           'val': evaluate_model(model, val_loader)
       }
       
       print(f'Epoch {epoch+1}:')
       print(f'Loss: {avg_loss:.4f}')
       print(f'Train - MAE: {results["train"]["mae"]:.4f}, Pearson_r: {results["train"]["pearson_r"]:.4f}')
       print(f'Val - MAE: {results["val"]["mae"]:.4f}, Pearson_r: {results["val"]["pearson_r"]:.4f}')
       
       if results['val']['mae'] < best_val_loss:
           best_val_loss = results['val']['mae']
           torch.save({
               'epoch': epoch,
               'model_state_dict': model.state_dict(),
               'optimizer_state_dict': optimizer.state_dict(),
               'loss': losses,
               'val_loss': best_val_loss,
           }, 'best_model_stage2.pth')
           
   # 绘制loss曲线
   plt.figure(figsize=(10,5))
   plt.plot(losses)  
   plt.title('Training Loss')
   plt.xlabel('Epoch')
   plt.ylabel('Loss')
   plt.savefig('loss_curve.png')
   plt.close()

def save_checkpoint(model, optimizer, epoch, loss, val_loss, save_dir, stage, celltype):
    """
    保存模型检查点
    """
    os.makedirs(save_dir, exist_ok=True)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    filename = os.path.join(save_dir, f'{celltype}_stage{stage}_checkpoint_{timestamp}.pth')
    
    checkpoint = {
        'epoch': epoch,
        'model_state_dict': model.state_dict(),
        'optimizer_state_dict': optimizer.state_dict() if optimizer else None,
        'loss': loss,
        'val_loss': val_loss,
        'lasso': model.lasso if stage == 1 else None,
        'scaler': model.scaler if stage == 1 else None
    }
    
    torch.save(checkpoint, filename)
    return filename

def load_checkpoint(model, checkpoint_path, optimizer=None):
    """
    加载模型检查点
    """
    checkpoint = torch.load(checkpoint_path)
    model.load_state_dict(checkpoint['model_state_dict'])
    if optimizer and 'optimizer_state_dict' in checkpoint:
        optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
    return checkpoint['epoch'], checkpoint.get('loss', [])

def main():
    parser = argparse.ArgumentParser(description='Two-stage age predictor training')
    
    # Data parameters
    parser.add_argument('--data_dir', type=str, required=True,
                      help='Base directory for data')
    parser.add_argument('--celltype', type=str, required=True,
                      help='Cell type for analysis')
    parser.add_argument('--output_dir', type=str, default='./outputs',
                      help='Directory for saving outputs')
    
    # Training parameters
    parser.add_argument('--seed', type=int, default=42,
                      help='Random seed')
    parser.add_argument('--batch_size', type=int, default=20,
                      help='Batch size for training')
    parser.add_argument('--num_workers', type=int, default=4,
                      help='Number of workers for data loading')
    
    # Stage 1 parameters
    parser.add_argument('--n_alphas', type=int, default=100,
                      help='Number of alphas to test in LassoCV')
    parser.add_argument('--cv_folds', type=int, default=5,
                      help='Number of CV folds for LassoCV')
    parser.add_argument('--min_alpha', type=float, default=1e-4,
                      help='Minimum alpha value for LassoCV')
    parser.add_argument('--max_alpha', type=float, default=1e2,
                      help='Maximum alpha value for LassoCV')
    
    # Stage 2 parameters
    parser.add_argument('--lr2', type=float, default=0.1,
                      help='Learning rate for stage 2')
    parser.add_argument('--epochs2', type=int, default=5,
                      help='Number of epochs for stage 2')
    parser.add_argument('--weight_decay', type=float, default=0.01,
                      help='Weight decay for optimizer')
    parser.add_argument('--mae_weight', type=float, default=0.7,
                      help='Weight for MAE loss')
    parser.add_argument('--mse_weight', type=float, default=0.3,
                      help='Weight for MSE loss')
    
    # Model parameters
    parser.add_argument('--checkpoint_dir', type=str, default='./checkpoints',
                      help='Directory for saving model checkpoints')
    parser.add_argument('--resume', type=str, default=None,
                      help='Path to checkpoint to resume from')
    
    args = parser.parse_args()
    
    # Create output directories
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(args.checkpoint_dir, exist_ok=True)
    
    # Save configuration
    config_path = os.path.join(args.output_dir, f'{args.celltype}_config.json')
    with open(config_path, 'w') as f:
        json.dump(vars(args), f, indent=4)
    
    # Set random seed
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(args.seed)
        device = torch.device('cuda')
    else:
        device = torch.device('cpu')
    
    # Load data
    age_labels_path = os.path.join(args.data_dir, f'{args.celltype}_age_splitted.json')
    with open(age_labels_path, 'r') as f:
        age_labels = json.load(f)
    
    # Create datasets
    processed_data_path = os.path.join(args.data_dir, 'processed')
    train_dataset = SingleCellDataset(
        os.path.join(processed_data_path, f'{args.celltype}_train_aligned_top2000_var.csv'),
        age_labels, 'train'
    )
    val_dataset = SingleCellDataset(
        os.path.join(processed_data_path, f'{args.celltype}_val_aligned_top2000_var.csv'),
        age_labels, 'val'
    )
    test_dataset = SingleCellDataset(
        os.path.join(processed_data_path, f'{args.celltype}_test_aligned_top2000_var.csv'),
        age_labels, 'test'
    )
    
    # Create data loaders
    train_loader = DataLoader(
        train_dataset,
        batch_size=args.batch_size,
        shuffle=True,
        num_workers=args.num_workers
    )
    val_loader = DataLoader(
        val_dataset,
        batch_size=args.batch_size,
        num_workers=args.num_workers
    )
    test_loader = DataLoader(
        test_dataset,
        batch_size=args.batch_size,
        num_workers=args.num_workers
    )
    
    # Initialize model
    model = AgePredictorTwoStage(ckpt_path='./models/models.ckpt')
    model.build()
    model = model.to(device)
    
    if args.resume:
        start_epoch, losses = load_checkpoint(model, args.resume)
        print(f"Resuming from epoch {start_epoch}")
    
    # Stage 1: Train LassoCV
    print("Stage 1: Training LassoCV")
    model.freeze_pretrained()
    best_alpha = model.train_lasso_cv(
        train_loader,
        val_loader,
        n_alphas=args.n_alphas,
        cv=args.cv_folds
    )
    
    # Evaluate and save Stage 1
    stage1_results = {
        'train': evaluate_model(model, train_loader, device, 'Stage1 Train'),
        'val': evaluate_model(model, val_loader, device, 'Stage1 Val'),
        'test': evaluate_model(model, test_loader, device, 'Stage1 Test')
    }
    
    # Save Stage 1 results
    results_path = os.path.join(args.output_dir, f'{args.celltype}_stage1_results.json')
    with open(results_path, 'w') as f:
        json.dump(stage1_results, f, indent=4)
    
    # Save Stage 1 model
    stage1_checkpoint = save_checkpoint(
        model, None, 0, None, stage1_results['val']['mae'],
        args.checkpoint_dir, 1, args.celltype
    )
    
    # Stage 2: Fine-tune
    print("\nStage 2: Fine-tuning")
    model.unfreeze_last_encoder_layer()
    train_stage2(
        model,
        train_loader,
        val_loader,
        test_loader,
        args.epochs2,
        args.lr2,
        device=device
    )
    
    # Final evaluation
    final_results = {
        'train': evaluate_model(model, train_loader, device, 'Final Train'),
        'val': evaluate_model(model, val_loader, device, 'Final Val'),
        'test': evaluate_model(model, test_loader, device, 'Final Test')
    }
    
    # Save final results
    final_results_path = os.path.join(args.output_dir, f'{args.celltype}_final_results.json')
    with open(final_results_path, 'w') as f:
        json.dump(final_results, f, indent=4)

if __name__ == '__main__':
    main()
