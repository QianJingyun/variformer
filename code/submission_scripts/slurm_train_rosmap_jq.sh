#!/bin/bash
#SBATCH -p encore-gpu
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2-0:00:00
#SBATCH -o /home/jqian54/sulab/enformer_fine_tuning/logs/train_rosmap/%j.out
#SBATCH -e /home/jqian54/sulab/enformer_fine_tuning/logs/train_rosmap/%j.err
#SBATCH --job-name=train_rosmap

source /sulab/users/jqian54/miniconda3/etc/profile.d/conda.sh
conda activate variformer_env2

export CUBLAS_WORKSPACE_CONFIG=:4096:8
export WANDB_MODE=online
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

cd /home/jqian54/sulab/enformer_fine_tuning/code


config_path=$1
fold=$2
model_type=$3

python ./train_rosmap.py --config_path $config_path --fold $fold --model_type $model_type