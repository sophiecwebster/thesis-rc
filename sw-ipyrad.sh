# directory guide

# /n/holyscratch01/hopkins_lab/webster/
# ├── ipyrad/                             (all shell/python scripts and instructions)
# │    └── trimmomatic/                    (trimmomatic tools)
# └── radseq/                             (original zip file and (partially) unzipped fastq gzip files)
#     ├── demultiplex/                    (entirety of step 2b's results, all demultiplexed files in stacks and rem)
#     └── splits/                         (the 'aa, ab, etc' type fastq files, split for parallel processing)
#         └── zipped/                     (the files zipped using pigz and edited zip files after moving Ns to header)
#             └── reads_cat/              (part 2a.4, edited gzips concatenated)
#

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

### step 2b Demultiplexing

# SLURM job to run:

#!/bin/bash
#SBATCH -p shared       					# Partition to submit to
#SBATCH -n 1                   				# Number of cores
#SBATCH -N 1                   				# Ensure that all cores are on one machine
#SBATCH -t 3-0:00               			# Runtime in days-hours:minutes
#SBATCH --mem 16000              			# Memory in MB
#SBATCH -J stacks_dm             			# job name
#SBATCH -o Nstacks_example_%A.out     		# File to which standard out will be written
#SBATCH -e Nstacks_example_%A.err     		# File to which standard err will be written
#SBATCH --mail-type=ALL        				# Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=sophiewebster@college.harvard.edu		# Email to which notifications will be sent

# YOUR HOME DIRECTORY HERE
# this directory should include:
#	your raw reads files
#	a subdirectory named 'scripts' that contains '2a_run_stacks_1PR_example.sh' and '2a_run_stacks_2CF_example.sh'
#	a subdirectory named 'demultiplex' that contains your barcodes file
myhome="/n/holyscratch01/hopkins_lab/webster/ipyrad"

# YOUR PREFIX HERE
# this should match the beginning of the file names of your unzipped raw reads, which should end with '_R1.fastq' or '_R2.fastq'
prefix="Meghan_S6_L002"

# YOUR BARCODES FILE HERE
# column 1 = PstI barcode [tab], column 2 = MspI barcode [tab], column 3 = sample_name
barcodes="2b_sw_barcodes.txt"


cd $myhome

module purge
module load gcc/7.1.0-fasrc01 stacks/2.4-fasrc01


# call the two scripts to perform the two demultiplexing functions in STACKS (process radtags, clone filtering)
# make sure these scripts are named correctly, are in the correct directory ($myhome/scripts), and are executable

./2b_run_stacks_1PR_example.sh \
  demultiplex/"$prefix"_edited_R1.fastq \
  demultiplex/"$prefix"_edited_R2.fastq \
  demultiplex/$barcodes

./2b_run_stacks_2CF_example.sh \
  demultiplex/stacks_1PR

  # ----------------------------- #

  #Let's see what went wrong here!

  Processing paired-end data.
Using Phred+33 encoding for quality scores.
Found 1 paired input file(s).
Searching for single and paired-end, inlined barcodes.
Invalid barcode on line 1: '@A00794:10:H7TL2DRXX:2:2101:1217:1031 1:N:0:ATTCG'
mv: cannot stat ‘demultiplex/stacks_1PR/*.rem.*’: No such file or directory
mv: cannot stat ‘demultiplex/stacks_1PR/process_radtags*’: No such file or directory
mv: cannot stat ‘demultiplex/stacks_1PR/*.1.fq’: No such file or directory
mv: cannot stat ‘demultiplex/stacks_1PR/*.2.fq’: No such file or directory
Processing paired-end data.
Searching for index oligo (i7 Illumina read).
Found 1 paired input file(s).
Processing file 1 of 1 [*_R1_*]
  Reading data from:
  demultiplex/stacks_1PR/*_R1_* and
  demultiplex/stacks_1PR/*_R2_*
Error opening input file 'demultiplex/stacks_1PR/*_R1_*'
Error opening input file 'demultiplex/stacks_1PR/*_R2_*'
Aborted. (basic_string::substr: __pos (which is 18446744073709551615) > this->size() (which is 6))
mv: cannot stat ‘demultiplex/stacks_2CF/*.1.fq’: No such file or directory
mv: cannot stat ‘demultiplex/stacks_2CF/*.2.fq’: No such file or directory

# ------- #

# okay, running this again, this time on just Meghan_S6_L001_edited_R1.fastq and Meghan_S6_L001_edited_R2.fastq
# this is the output of 2b_run_stacks_1PR_example.sh ??? how'd we do

47753136 total sequences
 9792218 barcode not found drops (20.5%)
     465 low quality read drops (0.0%)
  887049 RAD cutsite not found drops (1.9%)
37073404 retained reads (77.6%)

# ------ #

# output of 2b_run_stacks_1PR_example.sh on Meghan_S6_L002_edited_R1.fastq and Meghan_S6_L002_edited_R2.fastq
# only thing that seemingly went wrong was this?
# mv: cannot move ‘demultiplex/stacks_1PR/process_radtags.demultiplex.log’ to
# 'logs/Meghan_S6_L002_edited_PR.log’: No such file or directory -- hopefully not a problem?
# this happened both times (with *_L001_* and *_L002_*)

50870570 total sequences
10510108 barcode not found drops (20.7%)
    1546 low quality read drops (0.0%)
  936557 RAD cutsite not found drops (1.8%)
39422359 retained reads (77.5%)

# -------- #

# step 2c: trimming
# using 2c_write_trimmomatic.xlsx to create giant series of commands for each demultiplexed file.

# renamed the files with these scripts, for the L001 and L002 files respectively:
# (just appending L001 or L002 on the end because otherwise their numbers are the same
# (because i ran them separately in step 2b!))

for f in *.fastq; do mv $f "${f%.fastq}_L001.fastq"; done
for f in *.fastq; do mv $f "${f%.fastq}_L002.fastq"; done

# using CF files (clone filtered), which have tossed out any PCR clones
# PR files are just demultipexed (standing for 'process radtags')
# so, using the contents of stacks_2CF_L001 and stacks_2CF_L002 in 2c_write_trimmomatic.xlsx
