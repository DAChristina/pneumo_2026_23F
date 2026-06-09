# things to do:
# 1. report NT back to Harry


# 2. compile all 23F

# 3. construct script for running gubbins as submitted job
Create an index file:
  <See https://github.com/nickjcroucher/gubbins/blob/master/docs/gubbins_tutorial.md>
  (gubbinsTRIAL) [dac23@login-e GPSC31_n739]$ for file in *.fasta; do echo -e "${file%.fasta}\t${file}"; done > AllTarget_isolates.list

Get the genome reference after creating an index list:
  For GPSC31:
  (gubbinsTRIAL) [dac23@login-e GPSC31_n739]$ wget -q -O - https://www.ebi.ac.uk/ena/browser/api/fasta/FQ312030.1?download=true > Strep_INV104_Ref.fasta

For GPSC2:
  (gubbinsTRIAL) [dac23@login-e GPSC2_n7]$ wget -q -O - https://www.ebi.ac.uk/ena/browser/api/fasta/CP000920.1?download=true > Strep_P1031_Ref.fasta

Run alignment (using SKA2 that has already contained inside gubbins instead of SKA2)
(gubbinsTRIAL) [dac23@login-e GPSC31_n739]$ generate_ska_alignment.py --reference Strep_INV104_Ref.fasta --input AllTarget_isolates.list --out GPSC31_n739.aln

(gubbinsTRIAL) [dac23@login-e GPSC31_n739]$ run_gubbins.py --prefix n739_ GPSC31_n739.aln

# FAILED because I run the file on HPC shell! I should have submitted a job instead:
(gubbinsTRIAL) [dac23@login-e GPSC31_n739]$ chmod +x submit.pbs
#!/bin/bash

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=8:mem=64gb

# Job sizing guidance:imperial college's Intro to HPC OR
# https://www.imperial.ac.uk/computational-methods/hpc/
source ~/anaconda3/etc/profile.d/conda.sh

cd $PBS_O_WORKDIR
#mkdir /trialFolder # Failed
conda run -n gubbinsTRIAL run_gubbins.py --prefix n739_ GPSC31_n739.aln > myresult.txt
#mv *.pdf myresult.txt $PBS_O_WORKDIR

printenv > myenv.txt
#mv myenv.txt $PBS_O_WORKDIR
(gubbinsTRIAL) [dac23@login-e GPSC31_n739]$ qsub submit.pbs
9668402.pbs
(gubbinsTRIAL) [dac23@login-e GPSC31_n739]$ qstat -u dac23








# age stratification due to weird Age_years column 