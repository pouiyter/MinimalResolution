#!/bin/bash
#SBATCH -J test                   # 作业名为 test
#SBATCH --cpus-per-task=127
#SBATCH -p node
#SBATCH --nodelist=node1
#SBATCH -o out1
#SBATCH -e err1
#module load intel/
# 输入要执行的命令，例如 ./hello 或 python test.py 
./mr_st 35 36
./BPtab 35
./mr_BP 35 35
