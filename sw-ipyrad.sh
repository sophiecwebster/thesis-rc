# directory guide

# /n/holyscratch01/hopkins_lab/webster/
# ├── ipyrad/                             (all shell/python scripts and instructions)
# │    └── trimmomatic/                    (trimmomatic tools)
# └── radseq/                             (original zip file and (partially) unzipped fastq gzip files)
#     ├── demultiplex/                    (entirety of step 2b's results, all demultiplexed files in stacks and rem)
#     ├── trimmed/                        (products of step 2c: trimmomatic)
#     │   ├── trim_paired/
#     │   └── trim_unpaired/
#     ├── ipyRAD_assembly01/              (houses the assembly)
#     ├── sorted_fastqs/                  (trim_paired fastq files, ready for assembly!)
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

# TODO before you run: PUT ALL CF FILES IN SAME FOLDER AND RUN FROM THERE !!
# check what the names were ?

# to fix my foolish error of changing the names twice, i used this:

for f in A* ; do
    b=${f:0:15}${f:20:500}
    mv "$f" "$b"
done

# to rename the files for R1/R2 paired file column:

for f in 1* ; do
  mv $f "${f:0:8}${f:11:15}";
done

for f in 1* ; do
    mv $f "${f}_.fastq";
done

for f in 1* ; do
  mv $f ${f:0:12}${f:18:500};
done

# same thing, this time for files starting with "AA2"

for f in A* ; do
  mv $f ${f:0:7}${f:10:14};
done

for f in A* ; do
    mv $f "${f}_.fastq";
done

for f in A* ; do
  mv $f ${f:0:11}${f:17:500};
done

# fix names for R1/R2 unpaired columns

for f in *.fastq ; do
  mv $f "${f%.fastq}_unpaired.fastq";
done

# modifications to 2c_trimmomatic_example.#!/bin/sh

#!/bin/bash
#SBATCH -p shared 		      				# Partition to submit to
#SBATCH -n 16                 				# Number of cores
#SBATCH -N 1                  				# Ensure that all cores are on one machine
#SBATCH -t 1-0:00             				# Runtime in days-hours:minutes
#SBATCH --mem 32000           				# Memory in MB
#SBATCH -J trimmomatic            			# job name
#SBATCH -o trimmomatic_%A.out        		# File to which standard out will be written
#SBATCH -e trimmomatic_%A.err        		# File to which standard err will be written
#SBATCH --mail-type=ALL      				# Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=sophiewebster@college.harvard.edu     # Email to which notifications will be sent

# load a java module
module load jdk/10.0.1-fasrc01

# move to working directory
cd /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic

# make directories for paired and unpaired reads
mkdir trim_paired
mkdir trim_unpaired

# give R1 and R2 and then give output names for R1 paired, R1 unpaired, R2 paired and R2 unpaired.
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1011_R1_CF_L001.fastq 1011_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1011_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1011_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1011_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1011_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1012_R1_CF_L001.fastq 1012_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1012_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1012_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1012_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1012_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1013_R1_CF_L001.fastq 1013_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1013_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1013_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1013_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1013_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1014_R1_CF_L001.fastq 1014_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1014_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1014_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1014_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1014_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1015_R1_CF_L001.fastq 1015_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1015_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1015_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1015_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1015_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1016_R1_CF_L001.fastq 1016_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1016_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1016_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1016_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1016_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1017_R1_CF_L001.fastq 1017_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1017_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1017_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1017_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1017_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1018_R1_CF_L001.fastq 1018_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1018_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1018_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1018_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1018_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1019_R1_CF_L001.fastq 1019_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1019_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1019_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1019_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1019_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1020_R1_CF_L001.fastq 1020_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1020_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1020_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1020_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1020_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1021_R1_CF_L001.fastq 1021_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1021_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1021_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1021_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1021_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1022_R1_CF_L001.fastq 1022_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1022_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1022_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1022_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1022_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1023_R1_CF_L001.fastq 1023_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1023_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1023_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1023_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1023_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1025_R1_CF_L001.fastq 1025_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1025_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1025_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1025_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1025_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1026_R1_CF_L001.fastq 1026_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1026_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1026_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1026_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1026_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1027_R1_CF_L001.fastq 1027_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1027_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1027_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1027_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1027_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1028_R1_CF_L001.fastq 1028_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1028_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1028_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1028_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1028_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1029_R1_CF_L001.fastq 1029_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1029_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1029_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1029_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1029_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1030_R1_CF_L001.fastq 1030_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1030_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1030_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1030_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1030_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1031_R1_CF_L001.fastq 1031_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1031_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1031_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1031_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1031_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1032_R1_CF_L001.fastq 1032_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1032_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1032_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1032_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1032_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1033_R1_CF_L001.fastq 1033_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1033_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1033_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1033_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1033_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1034_R1_CF_L001.fastq 1034_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1034_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1034_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1034_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1034_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1035_R1_CF_L001.fastq 1035_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1035_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1035_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1035_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1035_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1036_R1_CF_L001.fastq 1036_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1036_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1036_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1036_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1036_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1037_R1_CF_L001.fastq 1037_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1037_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1037_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1037_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1037_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1038_R1_CF_L001.fastq 1038_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1038_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1038_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1038_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1038_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1039_R1_CF_L001.fastq 1039_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1039_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1039_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1039_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1039_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1040_R1_CF_L001.fastq 1040_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1040_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1040_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1040_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1040_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1041_R1_CF_L001.fastq 1041_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1041_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1041_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1041_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1041_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1042_R1_CF_L001.fastq 1042_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1042_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1042_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1042_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1042_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1043_R1_CF_L001.fastq 1043_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1043_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1043_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1043_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1043_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1044_R1_CF_L001.fastq 1044_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1044_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1044_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1044_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1044_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1045_R1_CF_L001.fastq 1045_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1045_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1045_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1045_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1045_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1046_R1_CF_L001.fastq 1046_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1046_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1046_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1046_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1046_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1047_R1_CF_L001.fastq 1047_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1047_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1047_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1047_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1047_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1048_R1_CF_L001.fastq 1048_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1048_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1048_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1048_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1048_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1049_R1_CF_L001.fastq 1049_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1049_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1049_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1049_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1049_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1050_R1_CF_L001.fastq 1050_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1050_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1050_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1050_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1050_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1051_R1_CF_L001.fastq 1051_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1051_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1051_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1051_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1051_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1052_R1_CF_L001.fastq 1052_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1052_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1052_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1052_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1052_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1053_R1_CF_L001.fastq 1053_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1053_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1053_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1053_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1053_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1054_R1_CF_L001.fastq 1054_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1054_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1054_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1054_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1054_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1055_R1_CF_L001.fastq 1055_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1055_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1055_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1055_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1055_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1056_R1_CF_L001.fastq 1056_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1056_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1056_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1056_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1056_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1057_R1_CF_L001.fastq 1057_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1057_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1057_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1057_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1057_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1058_R1_CF_L001.fastq 1058_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1058_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1058_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1058_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1058_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1059_R1_CF_L001.fastq 1059_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1059_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1059_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1059_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1059_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1060_R1_CF_L001.fastq 1060_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1060_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1060_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1060_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1060_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1061_R1_CF_L001.fastq 1061_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1061_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1061_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1061_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1061_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1062_R1_CF_L001.fastq 1062_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1062_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1062_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1062_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1062_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1063_R1_CF_L001.fastq 1063_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1063_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1063_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1063_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1063_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1064_R1_CF_L001.fastq 1064_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1064_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1064_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1064_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1064_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1065_R1_CF_L001.fastq 1065_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1065_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1065_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1065_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1065_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1066_R1_CF_L001.fastq 1066_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1066_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1066_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1066_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1066_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1067_R1_CF_L001.fastq 1067_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1067_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1067_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1067_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1067_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1068_R1_CF_L001.fastq 1068_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1068_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1068_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1068_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1068_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1069_R1_CF_L001.fastq 1069_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1069_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1069_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1069_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1069_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1070_R1_CF_L001.fastq 1070_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1070_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1070_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1070_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1070_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1071_R1_CF_L001.fastq 1071_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1071_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1071_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1071_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1071_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1072_R1_CF_L001.fastq 1072_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1072_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1072_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1072_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1072_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1073_R1_CF_L001.fastq 1073_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1073_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1073_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1073_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1073_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1074_R1_CF_L001.fastq 1074_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1074_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1074_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1074_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1074_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1075_R1_CF_L001.fastq 1075_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1075_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1075_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1075_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1075_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1076_R1_CF_L001.fastq 1076_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1076_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1076_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1076_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1076_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1077_R1_CF_L001.fastq 1077_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1077_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1077_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1077_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1077_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 AA2_R1_CF_L001.fastq AA2_R2_CF_L001.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/AA2_R1_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/AA2_R1_CF_L001_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/AA2_R2_L001_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/AA2_R2_CF_L001_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1011_R1_CF_L002.fastq 1011_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1011_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1011_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1011_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1011_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1012_R1_CF_L002.fastq 1012_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1012_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1012_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1012_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1012_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1013_R1_CF_L002.fastq 1013_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1013_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1013_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1013_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1013_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1014_R1_CF_L002.fastq 1014_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1014_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1014_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1014_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1014_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1015_R1_CF_L002.fastq 1015_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1015_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1015_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1015_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1015_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1016_R1_CF_L002.fastq 1016_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1016_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1016_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1016_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1016_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1017_R1_CF_L002.fastq 1017_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1017_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1017_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1017_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1017_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1018_R1_CF_L002.fastq 1018_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1018_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1018_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1018_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1018_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1019_R1_CF_L002.fastq 1019_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1019_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1019_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1019_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1019_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1020_R1_CF_L002.fastq 1020_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1020_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1020_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1020_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1020_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1021_R1_CF_L002.fastq 1021_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1021_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1021_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1021_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1021_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1022_R1_CF_L002.fastq 1022_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1022_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1022_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1022_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1022_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1023_R1_CF_L002.fastq 1023_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1023_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1023_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1023_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1023_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1025_R1_CF_L002.fastq 1025_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1025_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1025_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1025_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1025_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1026_R1_CF_L002.fastq 1026_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1026_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1026_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1026_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1026_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1027_R1_CF_L002.fastq 1027_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1027_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1027_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1027_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1027_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1028_R1_CF_L002.fastq 1028_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1028_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1028_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1028_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1028_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1029_R1_CF_L002.fastq 1029_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1029_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1029_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1029_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1029_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1030_R1_CF_L002.fastq 1030_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1030_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1030_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1030_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1030_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1031_R1_CF_L002.fastq 1031_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1031_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1031_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1031_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1031_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1032_R1_CF_L002.fastq 1032_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1032_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1032_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1032_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1032_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1033_R1_CF_L002.fastq 1033_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1033_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1033_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1033_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1033_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1034_R1_CF_L002.fastq 1034_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1034_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1034_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1034_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1034_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1035_R1_CF_L002.fastq 1035_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1035_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1035_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1035_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1035_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1036_R1_CF_L002.fastq 1036_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1036_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1036_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1036_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1036_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1037_R1_CF_L002.fastq 1037_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1037_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1037_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1037_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1037_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1038_R1_CF_L002.fastq 1038_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1038_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1038_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1038_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1038_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1039_R1_CF_L002.fastq 1039_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1039_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1039_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1039_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1039_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1040_R1_CF_L002.fastq 1040_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1040_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1040_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1040_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1040_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1041_R1_CF_L002.fastq 1041_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1041_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1041_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1041_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1041_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1042_R1_CF_L002.fastq 1042_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1042_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1042_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1042_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1042_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1043_R1_CF_L002.fastq 1043_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1043_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1043_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1043_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1043_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1044_R1_CF_L002.fastq 1044_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1044_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1044_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1044_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1044_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1045_R1_CF_L002.fastq 1045_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1045_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1045_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1045_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1045_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1046_R1_CF_L002.fastq 1046_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1046_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1046_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1046_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1046_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1047_R1_CF_L002.fastq 1047_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1047_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1047_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1047_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1047_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1048_R1_CF_L002.fastq 1048_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1048_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1048_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1048_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1048_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1049_R1_CF_L002.fastq 1049_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1049_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1049_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1049_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1049_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1050_R1_CF_L002.fastq 1050_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1050_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1050_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1050_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1050_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1051_R1_CF_L002.fastq 1051_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1051_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1051_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1051_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1051_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1052_R1_CF_L002.fastq 1052_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1052_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1052_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1052_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1052_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1053_R1_CF_L002.fastq 1053_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1053_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1053_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1053_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1053_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1054_R1_CF_L002.fastq 1054_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1054_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1054_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1054_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1054_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1055_R1_CF_L002.fastq 1055_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1055_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1055_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1055_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1055_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1056_R1_CF_L002.fastq 1056_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1056_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1056_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1056_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1056_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1057_R1_CF_L002.fastq 1057_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1057_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1057_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1057_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1057_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1058_R1_CF_L002.fastq 1058_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1058_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1058_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1058_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1058_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1059_R1_CF_L002.fastq 1059_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1059_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1059_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1059_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1059_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1060_R1_CF_L002.fastq 1060_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1060_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1060_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1060_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1060_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1061_R1_CF_L002.fastq 1061_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1061_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1061_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1061_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1061_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1062_R1_CF_L002.fastq 1062_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1062_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1062_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1062_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1062_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1063_R1_CF_L002.fastq 1063_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1063_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1063_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1063_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1063_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1064_R1_CF_L002.fastq 1064_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1064_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1064_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1064_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1064_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1065_R1_CF_L002.fastq 1065_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1065_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1065_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1065_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1065_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1066_R1_CF_L002.fastq 1066_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1066_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1066_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1066_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1066_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1067_R1_CF_L002.fastq 1067_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1067_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1067_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1067_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1067_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1068_R1_CF_L002.fastq 1068_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1068_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1068_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1068_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1068_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1069_R1_CF_L002.fastq 1069_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1069_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1069_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1069_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1069_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1070_R1_CF_L002.fastq 1070_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1070_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1070_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1070_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1070_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1071_R1_CF_L002.fastq 1071_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1071_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1071_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1071_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1071_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1072_R1_CF_L002.fastq 1072_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1072_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1072_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1072_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1072_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1073_R1_CF_L002.fastq 1073_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1073_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1073_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1073_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1073_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1074_R1_CF_L002.fastq 1074_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1074_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1074_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1074_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1074_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1075_R1_CF_L002.fastq 1075_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1075_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1075_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1075_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1075_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1076_R1_CF_L002.fastq 1076_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1076_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1076_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1076_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1076_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 1077_R1_CF_L002.fastq 1077_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1077_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1077_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/1077_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/1077_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 AA2_R1_CF_L002.fastq AA2_R2_CF_L002.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/AA2_R1_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/AA2_R1_CF_L002_unpaired.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_paired/AA2_R2_L002_.fastq /n/holyscratch01/hopkins_lab/webster/ipyrad/trimmomatic/trim_unpaired/AA2_R2_CF_L002_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# -------------------------------------------- #

# step 3: running ipyrad !
# locale is /n/holyscratch01/hopkins_lab/webster/radseq/ipyRAD_assembly01
# env is ipyrad01

### IMPORTANT NOTICES ###
# 1. Do not include parameters [2] or [3] (you've already demultiplexed, silly!)
# 2. Make sure you control parallelization by hand with ipcluster (otherwise, Odyssey nodes seem to be starting too slowly
#    and ipyrad quits)
# 3. Make sure you've activated conda and are in your environment!

# run steps 1 and 2

#!/bin/bash

#SBATCH -n 20                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 2-00:00              # Runtime in D-HH:MM
#SBATCH -p shared		            # Partition to submit to
#SBATCH --mem=48000             # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ipyrad_s12_%A.out    # File to which STDOUT will be written
#SBATCH -e ipyrad_s12_%A.err    # File to which STDERR will be written
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=sophiewebster@college.harvard.edu      # Email to which notifications will be sent

module load Anaconda3/5.0.1-fasrc02
source activate ipyrad01
ipcluster start --n 20 --daemonize
sleep 60

cd /n/holyscratch01/hopkins_lab/webster/radseq/ipyRAD_assembly01

ipyrad -p params.txt -s 12 --ipcluster


# run step 3 (the Big Cahuna)

#!/bin/bash

#SBATCH -n 20                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 7-00:00              # Runtime in D-HH:MM
#SBATCH -p shared		        # Partition to submit to
#SBATCH --mem=64000             # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ipyrad_s3_%A.out      # File to which STDOUT will be written
#SBATCH -e ipyrad_s3_%A.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=sophiewebster@college.harvard.edu      # Email to which notifications will be sent

source activate ipyrad01
ipcluster start --n 20 --daemonize
sleep 60
cd /n/holyscratch01/hopkins_lab/webster/radseq/ipyRAD_assembly01

ipyrad -p params.txt -s 3 --ipcluster


# errors?
# edit: resolved by manually controlling ipcluster

Warning: technical replicates (same name) present.

-------------------------------------------------------------
ipyrad [v.0.9.50]
Interactive assembly and analysis of RAD-seq data
-------------------------------------------------------------

Encountered an Error.
Message:
ipcluster not found, use 'auto=True' or see docs.
[swebster@boslogin01 ipyrad]$ cd ..

 # AND #

ipyrad.assemble.utils.IPyradError:
            Could not find saved Assembly file (.json) in expected location.
            Checks in: [project_dir]/[assembly_name].json
            Checked: /n/holyscratch01/hopkins_lab/webster/radseq/ipyRAD_assembly01/assembly1.json


            So for example a sample called WatDo_PipPrep_R1_100.fq.gz must
                be referenced in the barcode file as WatDo_PipPrep_100.

# fixing my errors (this time, barcodes being named wrong cuz of L001 and L002)
# resolved because it was unnecessary to specify barcodes!
awk '$3=$3"L001_"'

# run steps 4 and 5

#!/bin/bash

#SBATCH -n 20                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 7-00:00              # Runtime in D-HH:MM
#SBATCH -p shared		        # Partition to submit to
#SBATCH --mem=64000             # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ipyrad_s3_%A.out      # File to which STDOUT will be written
#SBATCH -e ipyrad_s3_%A.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=sophiewebster@college.harvard.edu      # Email to which notifications will be sent

source activate ipyrad01
ipcluster start --n 20 --daemonize
sleep 60
cd /n/holyscratch01/hopkins_lab/webster/radseq/ipyRAD_assembly01

ipyrad -p params.txt -s 45 --ipcluster

## okay, looks like after steps 4 and 5, 54 (55?) samples got tossed.

# run steps 6 and 7

#!/bin/bash

#SBATCH -n 20                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 7-00:00              # Runtime in D-HH:MM
#SBATCH -p shared		        # Partition to submit to
#SBATCH --mem=64000             # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ipyrad_s3_%A.out      # File to which STDOUT will be written
#SBATCH -e ipyrad_s3_%A.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=sophiewebster@college.harvard.edu      # Email to which notifications will be sent

source activate ipyrad01
ipcluster start --n 20 --daemonize
sleep 60
cd /n/holyscratch01/hopkins_lab/webster/radseq/ipyRAD_assembly01

ipyrad -p params.txt -s 67 --ipcluster

# feedback from steps 6 & 7

Step 6: Clustering/Mapping across samples
skipping samples not in state==5:
['1012L001_', '1012L002_', '1013L001_', '1013L002_', '1014L001_', '1014L002_', '1015L001_', '1021L001_', '1021L002_', '1022L001_', '1022L002_', '1023L002_', '1025L001_', '1025L002_', '1029L001_', '1029L002_', '1031L001_', '1031L002_', '1032L001_', '1033L002_', '1036L002_', '1037L001_', '1037L002_', '1038L001_', '1040L001_', '1040L002_', '1041L002_', '1044L001_', '1044L002_', '1045L001_', '1045L002_', '1049L002_', '1050L001_', '1053L001_', '1053L002_', '1054L001_', '1054L002_', '1055L001_', '1055L002_', '1057L001_', '1057L002_', '1058L001_', '1060L002_', '1061L001_', '1061L002_', '1062L001_', '1062L002_', '1068L002_', '1069L002_', '1070L001_', '1071L001_', '1071L002_', '1076L001_', '1077L001_', '1077L002_', 'AA2L001_', 'AA2L002_']
[####################] 100% 0:00:10 | concatenating inputs
[####################] 100% 0:02:27 | clustering across
[####################] 100% 0:00:19 | building clusters
[####################] 100% 0:07:52 | aligning clusters

Step 7: Filtering and formatting output files

Encountered an Error.
Message:
There are samples in this assembly that were not present in step 6. This is
likely due to failed samples retained in the assembly from prior to step 5, or
branching/merging. The following samples are not in the step6 database:
{'1012L001_', '1029L001_', '1061L002_', '1031L001_', '1040L001_', '1033L002_', '1071L002_', '1023L002_', '1049L002_', '1076L001_', '1055L001_', '1022L002_', '1057L001_', '1031L002_', '1036L002_', '1021L002_', '1062L001_', '1040L002_', '1022L001_', '1069L002_', '1061L001_', '1054L002_', '1029L002_', '1053L001_', '1057L002_', '1077L001_', '1013L002_', '1044L002_', 'AA2L001_', '1071L001_', '1014L002_', 'AA2L002_', '1013L001_', '1068L002_', '1037L002_', '1058L001_', '1021L001_', '1054L001_', '1077L002_', '1044L001_', '1038L001_', '1055L002_', '1032L001_', '1050L001_', '1070L001_', '1045L002_', '1015L001_', '1037L001_', '1045L001_', '1012L002_', '1060L002_', '1025L002_', '1062L002_', '1041L002_', '1014L001_', '1025L001_', '1053L002_'}
Simplest solution is to branch and remove these from the assembly.

# yikes. a lot of these files seem quite small!

# simple loop for plucking out # of lines.

ls -lh | awk '{print $5" "$9}' | tail -n 135 > names.txt

for f in *.fastq ; do
  cat $f | wc -l;
done > check.txt

# for egrepping purposes, here are the misbegotten files (see line 592, or step 7 error message)
egrep "1012|1029|1061|1031|1040|1033|1071|1023|1049|1076|1055|1022|1057|1031|1036|1021|1062|1040|1022|1069|1061|1054|1029|1053|1057|1077|1013|1044|AA2|1071|1014|AA2|1013|1068|1037|1058|1021|1054|1077|1044|1038|1055|1032|1050|1070|1045|1015|1037|1045|1012|1060|1025|1062|1041|1014|1025|1053"

# running assembly 2
# steps 1 and 2 generates assembly2_edits, assembly2.json, and assembly2_s1_demultiplex_stats.txt
# step 3 generates assembly2_clust_0.85/
# step 4 and 5 generate assembly2-tmpdir/
