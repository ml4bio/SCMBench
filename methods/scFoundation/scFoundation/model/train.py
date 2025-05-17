from typing import Optional
import torch
from torch import nn
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import os
import json
from torch.optim.lr_scheduler import LambdaLR
from transformers import get_scheduler, SchedulerType

class CombinedLoss(nn.Module):
    def __init__(self, mae_weight=0.7, mse_weight=0.3):
        super().__init__()
        self.mae_loss = nn.L1Loss()
        self.mse_loss = nn.MSELoss()
        self.mae_weight = mae_weight
        self.mse_weight = mse_weight
        
    def forward(self, pred, target):
        mae = self.mae_loss(pred, target)
        mse = self.mse_loss(pred, target)
        return self.mae_weight * mae + self.mse_weight * mse

from torch.optim.lr_scheduler import LambdaLR
from transformers import get_scheduler, SchedulerType

class CosineWarmupScheduler(LambdaLR):
    def __init__(
        self,
        optimizer,
        num_training_steps: int,
        num_warmup_steps: Optional[int] = None,
        warmup_ratio: Optional[float] = None,
        last_epoch: int = -1,
        verbose: bool = False
    ):
        """
        Args:
            optimizer: Optimizer to schedule learning rate
            num_training_steps: Total number of training steps
            num_warmup_steps: Number of warmup steps
            warmup_ratio: Ratio of warmup steps to total steps
            last_epoch: The index of last epoch
            verbose: If True, prints a message to stdout for each update
        """
        if num_warmup_steps is None and warmup_ratio is None:
            raise ValueError("Either num_warmup_steps or warmup_ratio must be provided.")
        
        if num_warmup_steps is None:
            num_warmup_steps = int(num_training_steps * warmup_ratio)
            
        scheduler = get_scheduler(
            SchedulerType.COSINE,
            optimizer=optimizer,
            num_warmup_steps=num_warmup_steps,
            num_training_steps=num_training_steps
        )
        
        super().__init__(optimizer, scheduler.lr_lambdas[0], last_epoch, verbose=verbose)
        self.num_warmup_steps = num_warmup_steps
        self.num_training_steps = num_training_steps

def evaluate_model(model, data_loader, device='cuda', split_name=None, output_dir=None, method=None, stage=None, celltype=None):
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
    
    if split_name and output_dir and method and stage and celltype:
        save_visualization(predictions, targets, output_dir, method, stage, split_name, celltype)
    
    mae = np.mean(np.abs(predictions - targets))
    r2 = np.corrcoef(predictions, targets)[0,1] ** 2
    pearson_r = np.corrcoef(predictions, targets)[0,1]
    
    return {
        'mae': float(mae),
        'r2': float(r2),
        'pearson_r': float(pearson_r)
    }

def save_visualization(predictions, targets, output_dir, method, stage, split_name, celltype):
    """保存可视化结果"""
    plt.figure(figsize=(8, 6))
    plt.scatter(targets, predictions, alpha=0.5)
    plt.plot([min(targets), max(targets)], [min(targets), max(targets)], 'r--')
    plt.xlabel('True Age')
    plt.ylabel('Predicted Age')
    pearson_r = np.corrcoef(predictions, targets)[0,1]
    plt.title(f'{method.upper()} - Stage {stage} - {split_name}\nPearson r: {pearson_r:.3f}')
    
    filename = f'{celltype}_{method}_stage{stage}_{split_name.lower()}_scatter.png'
    plt.savefig(os.path.join(output_dir, filename))
    plt.close()

def save_checkpoint(model, optimizer, epoch, loss, val_loss, save_dir, stage, method, celltype):
    """保存模型检查点"""
    os.makedirs(save_dir, exist_ok=True)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    filename = os.path.join(save_dir, f'{celltype}_{method}_stage{stage}_checkpoint_{timestamp}.pth')
    
    checkpoint = {
        'epoch': epoch,
        'model_state_dict': model.state_dict(),
        'optimizer_state_dict': optimizer.state_dict() if optimizer else None,
        'loss': loss,
        'val_loss': val_loss,
        'method': method,
        'stage': stage
    }
    
    if method == 'lasso':
        checkpoint.update({
            'lasso': model.lasso if stage == 1 else None,
            'scaler': model.scaler if stage == 1 else None
        })
    
    torch.save(checkpoint, filename)
    return filename

# 修改 train_stage2 函数：
def train_stage2(model, train_loader, val_loader, test_loader, args, device='cuda'):
    """第二阶段训练（适用于两种方法）"""
    optimizer = torch.optim.AdamW(
        filter(lambda p: p.requires_grad, model.parameters()),
        lr=args.lr2,
        weight_decay=args.weight_decay
    )
    
    # 计算总训练步数
    num_training_steps = len(train_loader) * args.epochs2
    
    # 创建scheduler
    scheduler = CosineWarmupScheduler(
        optimizer=optimizer,
        num_training_steps=num_training_steps,
        warmup_ratio=0.1,  # 10% steps for warmup
        verbose=True
    )
    
    criterion = CombinedLoss(mae_weight=args.mae_weight, mse_weight=args.mse_weight)
    best_val_loss = float('inf')
    
    # 记录训练过程
    train_losses = []
    val_losses = []
    learning_rates = []
    
    print(f"\nStarting stage 2 training ({args.method})...")
    
    for epoch in range(args.epochs2):
        model.train()
        epoch_losses = []
        
        for batch_idx, batch in enumerate(train_loader):
            x = batch['x'].to(device)
            age = batch['age'].to(device)
            
            pred = model({'x': x})
            loss = criterion(pred, age)
            epoch_losses.append(loss.item())
            
            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()
            scheduler.step()  # Update learning rate
            
            # Record learning rate
            learning_rates.append(optimizer.param_groups[0]['lr'])
        
        avg_train_loss = np.mean(epoch_losses)
        train_losses.append(avg_train_loss)
        
        # Validation
        val_results = evaluate_model(model, val_loader, device)
        val_losses.append(val_results['mae'])
        
        # Log progress
        log_training_progress(epoch, {
            'train': evaluate_model(model, train_loader, device),
            'val': val_results
        }, args.method, 2)
        
        # Save best model
        if val_results['mae'] < best_val_loss:
            best_val_loss = val_results['mae']
            save_checkpoint(
                model, optimizer, epoch, train_losses, best_val_loss,
                args.checkpoint_dir, 2, args.method, args.celltype
            )
    
    # Plot training curves
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5))
    
    # Loss curves
    ax1.plot(train_losses, label='Train Loss')
    ax1.plot(val_losses, label='Val Loss')
    ax1.set_title('Training and Validation Loss')
    ax1.set_xlabel('Epoch')
    ax1.set_ylabel('Loss')
    ax1.legend()
    ax1.grid(True)
    
    # Learning rate curve
    ax2.plot(learning_rates)
    ax2.set_title('Learning Rate Schedule')
    ax2.set_xlabel('Step')
    ax2.set_ylabel('Learning Rate')
    ax2.set_yscale('log')
    ax2.grid(True)
    
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, 
                f'{args.celltype}_{args.method}_stage2_training_curves.png'))
    plt.close()
    
    return best_val_loss

def log_training_progress(epoch, results, method, stage):
    """输出训练进度日志"""
    print(f"\n{method.upper()} - Stage {stage} - Epoch {epoch+1}:")
    print("-" * 50)
    for split in ['train', 'val']:
        print(f"{split.capitalize()}:")
        print(f"  MAE: {results[split]['mae']:.4f}")
        print(f"  Pearson r: {results[split]['pearson_r']:.4f}")
        print(f"  R²: {results[split]['r2']:.4f}")
    print("-" * 50)

def save_results(results, output_dir, method, stage, celltype):
    """保存训练结果
    
    Args:
        results: 包含评估指标的字典
        output_dir: 输出目录路径
        method: 使用的方法 ('lasso' 或 'mlp')
        stage: 训练阶段 (1 或 2)
        celltype: 细胞类型名称
    """
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 构建结果文件路径
    results_path = os.path.join(
        output_dir, 
        f'{celltype}_{method}_stage{stage}_results.json'
    )
    
    # 添加元数据到结果中
    results.update({
        'method': method,
        'stage': stage,
        'celltype': celltype,
        'timestamp': datetime.now().strftime('%Y%m%d_%H%M%S')
    })
    
    # 保存结果
    with open(results_path, 'w') as f:
        json.dump(results, f, indent=4)
        
    print(f"\nResults saved to: {results_path}")
    
    return results_path

def train_model(args, model, train_loader, val_loader, test_loader, device):
    """完整的训练流程"""
    if args.method == 'lasso':
        # Lasso方法的训练流程
        # Stage 1: 训练Lasso
        print("Stage 1: Training Lasso")
        model.freeze_pretrained()
        best_alpha = model.train_lasso_cv(
            train_loader,
            val_loader,
            n_alphas=args.n_alphas,
            cv=args.cv_folds
        )
        
        # 评估Stage 1
        stage1_results = {
            'train': evaluate_model(model, train_loader, device, 'Stage1 Train', 
                                  args.output_dir, args.method, 1, args.celltype),
            'val': evaluate_model(model, val_loader, device, 'Stage1 Val',
                                args.output_dir, args.method, 1, args.celltype),
            'test': evaluate_model(model, test_loader, device, 'Stage1 Test',
                                 args.output_dir, args.method, 1, args.celltype)
        }
        
        # 保存Stage 1结果
        save_results(stage1_results, args.output_dir, args.method, 1, args.celltype)
        
        # Stage 2: Fine-tuning
        print("\nStage 2: Fine-tuning last encoder layer")
        model.unfreeze_last_encoder_layer()
        train_stage2(model, train_loader, val_loader, test_loader, args, device)
        
    else:  # MLP方法
        # Stage 1: 冻结预训练模型，训练MLP
        print("Stage 1: Training MLP, freezing pretrained model")
        model.freeze_pretrained()
        model.unfreeze_mlp()
        
        best_val_loss = train_stage2(model, train_loader, val_loader, test_loader, args, device)
        
        # 评估Stage 1
        stage1_results = {
            'train': evaluate_model(model, train_loader, device, 'Stage1 Train', 
                                  args.output_dir, args.method, 1, args.celltype),
            'val': evaluate_model(model, val_loader, device, 'Stage1 Val',
                                args.output_dir, args.method, 1, args.celltype),
            'test': evaluate_model(model, test_loader, device, 'Stage1 Test',
                                 args.output_dir, args.method, 1, args.celltype)
        }
        
        # 保存Stage 1结果
        save_results(stage1_results, args.output_dir, args.method, 1, args.celltype)
        
        # Stage 2: 冻结MLP，微调预训练模型的最后一层
        print("\nStage 2: Fine-tuning pretrained model, freezing MLP")
        model.freeze_mlp()
        model.unfreeze_last_encoder_layer()
        args.lr2 = args.lr2 * 0.1
        train_stage2(model, train_loader, val_loader, test_loader, args, device)
    
    # 最终评估（对两种方法都一样）
    final_results = {
        'train': evaluate_model(model, train_loader, device, 'Final Train',
                              args.output_dir, args.method, 2, args.celltype),
        'val': evaluate_model(model, val_loader, device, 'Final Val',
                            args.output_dir, args.method, 2, args.celltype),
        'test': evaluate_model(model, test_loader, device, 'Final Test',
                             args.output_dir, args.method, 2, args.celltype)
    }
    
    save_results(final_results, args.output_dir, args.method, 2, args.celltype)
    
    return final_results