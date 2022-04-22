#setwd(here('dglnrt_v2_simulation_partially_identified/talapas/'))



for(i in 1:100){
  
  a <- '#!/bin/bash'
  a <- rbind(a,'#SBATCH --account=edquant')
  a <- rbind(a,'#SBATCH --partition=long')
  a <- rbind(a,paste0('#SBATCH --job-name=rep',i))
  a <- rbind(a,paste0('#SBATCH --output=/gpfs/projects/edquant/cengiz/dglnrt_mis/rep',i,'.out'))
  a <- rbind(a,paste0('#SBATCH --error=/gpfs/projects/edquant/cengiz/dglnrt_mis/rep',i,'.err'))
  a <- rbind(a,'#SBATCH --time=10080')
  a <- rbind(a,'#SBATCH --mem=4000M')
  a <- rbind(a,'#SBATCH --nodes=1')
  a <- rbind(a,'#SBATCH --ntasks-per-node=1')
  a <- rbind(a,'#SBATCH --cpus-per-task=4')
  a <- rbind(a,"")
  a <- rbind(a,'module load gcc')
  a <- rbind(a,'module load R/4.0.2')
  a <- rbind(a,"")
  a <- rbind(a,paste0('R --no-save < /gpfs/projects/edquant/cengiz/dglnrt_mis/rep',i,'.r'))
  
  a <- noquote(a)
  
  write(a,paste0('/gpfs/projects/edquant/cengiz/dglnrt_mis/batch',i,'.batch'))
}
