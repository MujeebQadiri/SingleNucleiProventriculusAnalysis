#!/bin/bash

#SBATCH --job-name=FC_08073_count%j  # Job name
#SBATCH --nodes=1                                    # Run all processes on a single node	
#SBATCH --ntasks=2                                   # Run a single task		
#SBATCH --cpus-per-task=8                            # Number of CPU cores per task
#SBATCH --mem=40gb                                   # Job memory request
#SBATCH --time=08:10:00                              # Time limit hrs:min:sec
#SBATCH --partition=short 
#SBATCH --output=FC_08073_count%j.log                # Standard output and error log

module load gcc/9.2.0
module use /n/app/lmod/lmod/modulefiles/Compiler/gcc/9.2.0/
module load cmake/3.22.2
module load gdal/3.1.4
module load udunits/2.2.28
module load geos/3.10.2
module load fftw
module load cellranger/7.0.0 

date 
cd FC_08073/Data
date
cellranger count --id=BEC2 \
  --transcriptome=/n/groups/flyrnai/yifang/Projects/With/Afroditi/2021-10-19_FC_07023/cellranger-6.1.1_Drosophila_melanogaster.BDGP6.32.104/Drosophila_melanogaster.BDGP6.32.104 \
  --fastqs=/n/groups/flyrnai/dropbox/raw_datasets/FC_08073/Unaligned_12_PF_TenX_mm1/cellranger/mkfastq.outs/Project_tsudhirg/sudhir_may17_pool_BEC2_CAATCCCGAC \
  --sample=BEC2 

cellranger count --id=BEC1 \
  --transcriptome=/n/groups/flyrnai/yifang/Projects/With/Afroditi/2021-10-19_FC_07023/cellranger-6.1.1_Drosophila_melanogaster.BDGP6.32.104/Drosophila_melanogaster.BDGP6.32.104 \
  --fastqs=/n/groups/flyrnai/dropbox/raw_datasets/FC_08073/Unaligned_12_PF_TenX_mm1/cellranger/mkfastq.outs/Project_tsudhirg/sudhir_may17_pool_BEC1_TGCGCGGTTT \
  --sample=BEC1

date