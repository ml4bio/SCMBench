#!/bin/bash

# 创建存储最终结果的文件
echo "Learning Rate,Val MAE,Val R2,Train MAE,Train R2,Test MAE,Test R2" > lr_results.csv

# 定义要尝试的学习率数组
learning_rates=(0.00001 0.00005 0.0001 0.0005 0.001)

# 遍历每个学习率
for lr in "${learning_rates[@]}"
do
    echo "Training with learning rate: $lr"
    
    # 运行训练脚本
    python age-predictor.py --lr $lr
    
    # 从结果文件中提取指标
    val_mae=$(cat "final_results_lr_${lr}.json" | grep -o '"mae": [0-9.]*' | head -2 | tail -1 | cut -d' ' -f2)
    val_r2=$(cat "final_results_lr_${lr}.json" | grep -o '"r2": [0-9.]*' | head -2 | tail -1 | cut -d' ' -f2)
    train_mae=$(cat "final_results_lr_${lr}.json" | grep -o '"mae": [0-9.]*' | head -1 | cut -d' ' -f2)
    train_r2=$(cat "final_results_lr_${lr}.json" | grep -o '"r2": [0-9.]*' | head -1 | cut -d' ' -f2)
    test_mae=$(cat "final_results_lr_${lr}.json" | grep -o '"mae": [0-9.]*' | tail -1 | cut -d' ' -f2)
    test_r2=$(cat "final_results_lr_${lr}.json" | grep -o '"r2": [0-9.]*' | tail -1 | cut -d' ' -f2)
    
    # 将结果添加到CSV文件
    echo "$lr,$val_mae,$val_r2,$train_mae,$train_r2,$test_mae,$test_r2" >> lr_results.csv
done

# 找出最佳学习率（基于验证集R2）
best_lr=$(sort -t',' -k3 -n -r lr_results.csv | head -2 | tail -1 | cut -d',' -f1)
echo "Best learning rate (based on validation R2): $best_lr"

# 显示所有结果
echo "All results:"
cat lr_results.csv