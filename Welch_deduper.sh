#!/usr/bin/bash

#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=bgmp               #REQUIRED: which partition to use
#SBATCH --cpus-per-task=8            #optional: number of cpus, default is 1
#SBATCH --mem=16GB                        #optional: amount of memory, default is 4G
#SBATCH --time=1-00:00:00


conda activate base

sam_file='/projects/bgmp/shared/deduper/C1_SE_uniqAlign.sam'

sorted_sam='C1_SE_uniqAlign.sorted.sam'
output='output.sam'
/usr/bin/time -v \
./Welch_deduper.py -f $sorted_sam -o $output -u STL96.txt
