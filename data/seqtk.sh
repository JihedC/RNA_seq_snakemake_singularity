#!/bin/bash
#SBATCH -J seqtk
#SBATCH --mail-type=ALL
#SBATCH --mail-user j.chouaref@lumc.nl
#SBATCH -t 1:0:0
#SBATCH --mem=4000
#SBATCH --ntasks-per-node=1


DATA="/exports/humgen/jihed/Morc3_RNA/RNA-seq-snakemake/data"

echo "start job SLURM `date`"

for file in $DATA/*.gz;
  do
    echo $file;
    name=`basename -s '.fq.gz' $file`;
    echo $name
    seqtk sample -s100 $file 5000 > 'sub_'$name'.fq'
  done 

echo "Finished `date`"
  
