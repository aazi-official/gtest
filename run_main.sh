#!/bin/bash
#SBATCH --job-name=main_job            # 作业名称
#SBATCH --output=main_output.txt       # 标准输出文件
#SBATCH --error=main_error.txt         # 标准错误文件
#SBATCH --nodes=1                      # 使用的节点数量
#SBATCH --ntasks=4                     # 总任务数（核心数）
#SBATCH --time=01:00:00                # 预计运行时间（格式：hh:mm:ss）
#SBATCH --partition=batch              # 分区名称（默认是 batch）

# 加载必要的模块（根据集群环境调整）
module load openmpi//1.10.7

# 切换到作业提交的目录
cd $SLURM_SUBMIT_DIR

# 打印调试信息
echo "Running on node(s): $SLURM_NODELIST"
echo "Number of tasks: $SLURM_NTASKS"
echo "Working directory: $PWD"

# 运行已编译的程序
mpirun -np $SLURM_NTASKS ./main

