# directory guide

# /n/holyscratch01/hopkins_lab/webster/
# ├── ipyrad/                             (all shell/python scripts and instructions)
# └── radseq/                             (original zip file and (partially) unzipped fastq gzip files)
#     └── splits/                         (the 'aa, ab, etc' type fastq files, split for parallel processing)
#         └── zipped/                     (the files zipped using pigz and edited zip files after moving Ns to header)
#             └── reads_cat/              (part 2a.4, edited gzips concatenated)


### step 2a MoveNs
## part 1
# this unzips the GNU zip files and breaks each one down into files with a max of 30 million lines
# each of these is around ~2 GB

zcat Meghan_S6_L001_R1_001.fastq.gz | split -l 30000000 - Meghan_S6_L001_R1_001_ --additional-suffix=.fastq
zcat Meghan_S6_L001_R2_001.fastq.gz | split -l 30000000 - Meghan_S6_L001_R2_001_ --additional-suffix=.fastq
zcat Meghan_S6_L002_R1_001.fastq.gz | split -l 30000000 - Meghan_S6_L002_R1_001_ --additional-suffix=.fastq
zcat Meghan_S6_L002_R2_001.fastq.gz | split -l 30000000 - Meghan_S6_L002_R2_001_ --additional-suffix=.fastq

# ---------------------------------- #

## part 2
### splitting up the files and putting them into /n/holyscratch01/hopkins_lab/webster/radseq/splits

#!/bin/bash
#SBATCH -p serial_requeue	# Partition to submit to
#SBATCH -n 16               # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 2-0:00           # Runtime in days-hours:minutes
#SBATCH --mem 32000         # Memory in MB
#SBATCH -J pigz             # job name
#SBATCH -o pigz_%A.out         # File to which standard out will be written
#SBATCH -e pigz_%A.err         # File to which standard err will be written
#SBATCH --mail-type=ALL        				# Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=sophiewebster@college.harvard.edu	# Email to which notifications will be sent

source new-modules.sh
#module load GCC/8.2.0-2.31.1 pigz/2.4
module load GCCcore/8.2.0
module load pigz/2.4


cd /n/holyscratch01/hopkins_lab/webster/radseq/splits

# 'best' flag favors best compression over speed, 16 is number of processors used
# k flag keeps unzipped fastq files

pigz -k -p 16 --best Meghan_S6_L002_R1_001_ad.fastq Meghan_S6_L002_R2_001_aa.fastq Meghan_S6_L002_R2_001_ab.fastq Meghan_S6_L002_R2_001_ac.fastq Meghan_S6_L002_R2_001_ad.fastq

# ---------------------------------- #

## part 3

#!/bin/bash
#SBATCH -p shared			       			# Partition to submit to
#SBATCH -n 5	                   			# Number of cores
#SBATCH -N 1                   				# Ensure that all cores are on one machine
#SBATCH -t 7-00:00               			# Runtime in days-hours:minutes
#SBATCH --mem 64000              			# Memory in MB
#SBATCH -J move_Ns             				# job name
#SBATCH -o move_Ns_%A.out        			# File to which standard out will be written
#SBATCH -e move_Ns_%A.err        			# File to which standard err will be written
#SBATCH --mail-type=ALL        				# Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=[your email]			# Email to which notifications will be sent


source new-modules.sh
module load python/2.7.14-fasrc02

cd /n/holyscratch01/hopkins_lab/webster

# changed file names

# python path/to/Move_Ns_to_header_both_seqs.py path/to/split_raw_reads_R1_aa.fastq.gz path/to/split_raw_reads_R2_aa.fastq.gz &
python ipyrad/2a_Move_Ns_to_header_both_seqs.py radseq/splits/zipped/Meghan_S6_L001_R1_001_aa.fastq.gz radseq/splits/zipped/Meghan_S6_L001_R2_001_aa.fastq.gz &
python ipyrad/2a_Move_Ns_to_header_both_seqs.py radseq/splits/zipped/Meghan_S6_L001_R1_001_ab.fastq.gz radseq/splits/zipped/Meghan_S6_L001_R2_001_ab.fastq.gz &
python ipyrad/2a_Move_Ns_to_header_both_seqs.py radseq/splits/zipped/Meghan_S6_L001_R1_001_ac.fastq.gz radseq/splits/zipped/Meghan_S6_L001_R2_001_ac.fastq.gz &
python ipyrad/2a_Move_Ns_to_header_both_seqs.py radseq/splits/zipped/Meghan_S6_L001_R1_001_ad.fastq.gz radseq/splits/zipped/Meghan_S6_L001_R2_001_ad.fastq.gz &

python ipyrad/2a_Move_Ns_to_header_both_seqs.py radseq/splits/zipped/Meghan_S6_L002_R1_001_aa.fastq.gz radseq/splits/zipped/Meghan_S6_L002_R2_001_aa.fastq.gz &
python ipyrad/2a_Move_Ns_to_header_both_seqs.py radseq/splits/zipped/Meghan_S6_L002_R1_001_ab.fastq.gz radseq/splits/zipped/Meghan_S6_L002_R2_001_ab.fastq.gz &
python ipyrad/2a_Move_Ns_to_header_both_seqs.py radseq/splits/zipped/Meghan_S6_L002_R1_001_ac.fastq.gz radseq/splits/zipped/Meghan_S6_L002_R2_001_ac.fastq.gz &
python ipyrad/2a_Move_Ns_to_header_both_seqs.py radseq/splits/zipped/Meghan_S6_L002_R1_001_ad.fastq.gz radseq/splits/zipped/Meghan_S6_L002_R2_001_ad.fastq.gz &

wait

# -------------- #

# part 4
# concatenate edited gzip files

#!/bin/bash

#SBATCH -n 16                   # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 2-00:00              # Runtime in D-HH:MM
#SBATCH -p shared		        # Partition to submit to
#SBATCH --mem=32000             # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o cat_v11_%A.out     	# File to which STDOUT will be written
#SBATCH -e cat_v11_%A.err    	# File to which STDERR will be written
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=sophiewebster@college.harvard.edu     # Email to which notifications will be sent


cd /n/holyscratch01/hopkins_lab/webster/radseq/splits/zipped


cat ./Meghan_S6_L001_R1_001_aa.fastq.gzedited ./Meghan_S6_L001_R1_001_ab.fastq.gzedited ./Meghan_S6_L001_R1_001_ac.fastq.gzedited ./Meghan_S6_L001_R1_001_ad.fastq.gzedited > reads_cat/Meghan_S6_L001_edited_R1.fastq
cat ./Meghan_S6_L001_R2_001_aa.fastq.gzedited ./Meghan_S6_L001_R2_001_ab.fastq.gzedited ./Meghan_S6_L001_R2_001_ac.fastq.gzedited ./Meghan_S6_L001_R2_001_ad.fastq.gzedited > reads_cat/Meghan_S6_L001_edited_R2.fastq
cat ./Meghan_S6_L002_R1_001_aa.fastq.gzedited ./Meghan_S6_L002_R1_001_ab.fastq.gzedited ./Meghan_S6_L002_R1_001_ac.fastq.gzedited ./Meghan_S6_L002_R1_001_ad.fastq.gzedited > reads_cat/Meghan_S6_L002_edited_R1.fastq
cat ./Meghan_S6_L002_R2_001_aa.fastq.gzedited ./Meghan_S6_L002_R2_001_ab.fastq.gzedited ./Meghan_S6_L002_R2_001_ac.fastq.gzedited ./Meghan_S6_L002_R2_001_ad.fastq.gzedited > reads_cat/Meghan_S6_L002_edited_R2.fastq

# ---------- #
