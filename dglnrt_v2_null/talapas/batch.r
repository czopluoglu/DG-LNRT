

for(i in 2:100){
  
  a <- '#!/bin/bash'
  a <- rbind(a,'#SBATCH --account=edquant')
  a <- rbind(a,'#SBATCH --partition=longgpu')
  a <- rbind(a,'#SBATCH --job-name=rep1')
  a <- rbind(a,paste0('#SBATCH --output=/gpfs/projects/edquant/cengiz/dglnrt_null/rep',i,'.out'))
  a <- rbind(a,paste0('#SBATCH --error=/gpfs/projects/edquant/cengiz/dglnrt_null/rep',i,'.err'))
  a <- rbind(a,'#SBATCH --time=20160')
  a <- rbind(a,'#SBATCH --mem=4000M')
  a <- rbind(a,'#SBATCH --nodes=1')
  a <- rbind(a,'#SBATCH --ntasks-per-node=1')
  a <- rbind(a,'#SBATCH --cpus-per-task=4')
  a <- rbind(a,"")
  a <- rbind(a,'module load gcc')
  a <- rbind(a,'module load R/4.0.2')
  a <- rbind(a,"")
  a <- rbind(a,paste0('R --no-save < /gpfs/projects/edquant/cengiz/dglnrt_null/rep',i,'.r'))
  
  a <- noquote(a)
  
  write(a,paste0('batch',i,'.batch'))
}