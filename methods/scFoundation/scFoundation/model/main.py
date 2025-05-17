import argparse
import os
import json
import torch
import numpy as np
from torch.utils.data import Dataset, DataLoader
from datetime import datetime
import pandas as pd
from models import AgePredictorLasso, AgePredictorMLP
from train import train_model

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

def add_common_args(parser):
    """添加通用命令行参数"""
    # 数据相关参数
    parser.add_argument('--data_dir', type=str, required=True,
                       help='数据根目录')
    parser.add_argument('--celltype', type=str, required=True,
                       help='细胞类型')
    parser.add_argument('--output_dir', type=str, default='./outputs',
                       help='输出目录')
    
    # 训练相关参数
    parser.add_argument('--method', type=str, choices=['lasso', 'mlp'],
                       required=True, help='使用的方法')
    parser.add_argument('--seed', type=int, default=42,
                       help='随机种子')
    parser.add_argument('--batch_size', type=int, default=32,
                       help='批次大小')
    parser.add_argument('--num_workers', type=int, default=4,
                       help='数据加载的工作进程数')
    
    # 模型相关参数
    parser.add_argument('--checkpoint_dir', type=str, default='./checkpoints',
                       help='模型检查点保存目录')
    parser.add_argument('--model_path', type=str, default='./models/models.ckpt',
                       help='预训练模型路径')
    parser.add_argument('--resume', type=str, default=None,
                       help='恢复训练的检查点路径')

def setup_environment(args):
    """设置环境和随机种子"""
    # 设置随机种子
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(args.seed)
        device = torch.device('cuda')
    else:
        device = torch.device('cpu')
    
    # 创建输出目录
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(args.checkpoint_dir, exist_ok=True)
    
    return device

def load_data(args):
    """加载数据集"""
    # 加载年龄标签
    age_labels_path = os.path.join(args.data_dir, f'{args.celltype}_age_splitted.json')
    with open(age_labels_path, 'r') as f:
        age_labels = json.load(f)
    
    # 创建数据集
    processed_data_path = os.path.join(args.data_dir, 'processed')
    datasets = {
        split: SingleCellDataset(
            os.path.join(processed_data_path, f'{args.celltype}_{split}_aligned_top2000_var.csv'),
            age_labels,
            split
        ) for split in ['train', 'val', 'test']
    }
    
    # 创建数据加载器
    dataloaders = {
        'train': DataLoader(
            datasets['train'],
            batch_size=args.batch_size,
            shuffle=True,
            num_workers=args.num_workers
        ),
        'val': DataLoader(
            datasets['val'],
            batch_size=args.batch_size,
            num_workers=args.num_workers
        ),
        'test': DataLoader(
            datasets['test'],
            batch_size=args.batch_size,
            num_workers=args.num_workers
        )
    }
    
    return dataloaders

def initialize_model(args):
    """初始化模型"""
    if args.method == 'lasso':
        model = AgePredictorLasso(args.model_path)
    else:
        model = AgePredictorMLP(args.model_path, args.hidden_dims)
    
    model.build()
    return model

def save_config(args):
    """保存配置"""
    config_path = os.path.join(args.output_dir, 
                              f'{args.celltype}_{args.method}_config.json')
    with open(config_path, 'w') as f:
        json.dump(vars(args), f, indent=4)

def main():
    # 首先创建解析器并添加通用参数
    parser = argparse.ArgumentParser(description='年龄预测模型训练')
    add_common_args(parser)
    
    # 添加所有训练相关的通用参数
    parser.add_argument('--lr2', type=float, default=0.001,
                       help='Stage 2的学习率')
    parser.add_argument('--epochs2', type=int, default=20,
                       help='Stage 2的训练轮数')
    parser.add_argument('--weight_decay', type=float, default=0.01,
                       help='权重衰减')
    parser.add_argument('--mae_weight', type=float, default=0.7,
                       help='MAE损失的权重')
    parser.add_argument('--mse_weight', type=float, default=0.3,
                       help='MSE损失的权重')
    
    # 首先解析以获取method参数
    temp_args, _ = parser.parse_known_args()
    
    # 然后根据method添加特定参数
    if temp_args.method == 'lasso':
        parser.add_argument('--n_alphas', type=int, default=100,
                           help='LassoCV中测试的alpha数量')
        parser.add_argument('--cv_folds', type=int, default=5,
                           help='交叉验证折数')
        parser.add_argument('--min_alpha', type=float, default=1e-4,
                           help='最小alpha值')
        parser.add_argument('--max_alpha', type=float, default=1e2,
                           help='最大alpha值')
    elif temp_args.method == 'mlp':
        parser.add_argument('--hidden_dims', nargs='+', type=int,
                           default=[512, 256, 128],
                           help='MLP隐藏层维度')
        parser.add_argument('--dropout', type=float, default=0.2,
                           help='Dropout比率')
        parser.add_argument('--lr_stage1', type=float, default=0.001,
                           help='Stage 1 (MLP训练) 的学习率')
        parser.add_argument('--lr_stage2', type=float, default=0.0001,
                           help='Stage 2 (预训练模型微调) 的学习率')
    
        # 最后解析所有参数
    args = parser.parse_args()
    
    # 设置环境
    device = setup_environment(args)
    
    # 保存配置
    save_config(args)
    
    # 加载数据
    dataloaders = load_data(args)
    
    # 初始化模型
    model = initialize_model(args)
    model = model.to(device)
    
    # 如果需要恢复训练
    if args.resume:
        checkpoint = torch.load(args.resume)
        model.load_state_dict(checkpoint['model_state_dict'])
        print(f"Resumed from checkpoint: {args.resume}")
    
    # 训练模型
    final_results = train_model(
        args,
        model,
        dataloaders['train'],
        dataloaders['val'],
        dataloaders['test'],
        device
    )
    
    # 打印最终结果
    print("\nTraining completed!")
    print("\nFinal Results:")
    print("-" * 50)
    for split in ['train', 'val', 'test']:
        print(f"{split.capitalize()}:")
        print(f"  MAE: {final_results[split]['mae']:.4f}")
        print(f"  Pearson r: {final_results[split]['pearson_r']:.4f}")
        print(f"  R²: {final_results[split]['r2']:.4f}")
    print("-" * 50)

if __name__ == '__main__':
    main()