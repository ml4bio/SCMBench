import os
for i in [0,0.5,1,1.5,2,3,5,10]:
    os.system("python3 simulate_process.py --path datasets/simulated_num_cell_5000_num_batch_3_num_gene_3000_effect_" + str(i))
        # os.system("bash run_embed.sh simulated_num_cell_5000_num_batch_3_num_gene_3000_effect_" + str(i) + "  " + j)