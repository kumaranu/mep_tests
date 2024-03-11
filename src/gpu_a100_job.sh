#!/bin/bash -l
#SBATCH --job-name=test
#SBATCH --account=pc_chippschao
#SBATCH --partition=es1
#SBATCH --qos=es_lowprio
#SBATCH --nodes=1
#sBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:V100:1
#SBATCH --time=24:00:00

source activate quacc

mkdir -p a100

cp test1.py a100

cd a100
  python test1.py
cd ..

