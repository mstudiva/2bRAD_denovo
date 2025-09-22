## Analysis of Next-Generation Sequencing Data (ANGSD) pipeline, version September 17, 2025
# Based on 2bRAD pipeline by Ryan Eckert (reckert2017@fau.edu)
https://ryaneckert.github.io/Stephanocoenia_FKNMS_PopGen/code/
# Modified by Michael Studivan (studivanms@gmail.com)


#------------------------------
## Downloading scripts (if needed)

cd ~/bin/
git clone --recursive https://github.com/samtools/htslib.git
git clone https://github.com/ANGSD/angsd.git
cd htslib;make;cd ../angsd ;make HTSSRC=../htslib
make HTSSRC="systemwide"


#------------------------------
## Fuzzy Genotyping (ANGSD)

copy all *.bam* files to a new directory 'ANGSD'
mkdir ../ANGSD
cd ../ANGSD
mv ../mappedReads/*.bam* .

ls *bam >bamsClones

module load angsd-0.933-gcc-9.2.0-65d64pp

# Initial run ANGSD settings:
# -minMapQ 20: only highly unique mappings (prob of erroneous mapping =< 1%)
# -maxDepth: highest total depth (sum over all samples) to assess; set to 10x number of samples
# -minInd: the minimal number of individuals the site must be genotyped in. Reset to 50% of total N at this stage.
export FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -maxDepth 1030 -minInd 52"
export TODO="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"

echo '#!/bin/bash' >ofavDD.sh
echo angsd -b bamsClones -GL 1 $FILTERS $TODO -P 1 -out dd -nThreads 20 >>ofavDD.sh

sbatch --mem=200GB -o ofavDD.o%j -e ofavDD.e%j --mail-user=studivanms@gmail.com --mail-type=ALL --partition=shortq7 ofavDD.sh

# Summarizing results
module load R/3.6.1
echo '#!/bin/bash' >RQC.sh
echo Rscript ~/bin/plotQC.R prefix=dd >>RQC.sh
echo gzip -9 dd.counts >>RQC.sh
sbatch -e RQC.e%j -o RQC.o%j --dependency=afterok:460550 --mem=200GB --mail-type=ALL --mail-user=studivanms@gmail.com RQC.sh

cat quality.txt # proportion of sites covered at >5X

# scp dd.pdf to laptop to look at distribution of base quality scores, fraction of sites in each sample passing coverage thresholds and fraction of sites passing genotyping rates cutoffs. Use these to guide choices of -minQ, -minIndDepth and -minInd filters in subsequent ANGSD runs
cd local/directory/
scp mstudiva@koko-login.hpc.fau.edu:~/project/directory/ANGSD/\*.pdf .

# Second ANGSD run settings:
module load angsd-0.933-gcc-9.2.0-65d64pp
# change the -minInd flag to 80% of the number of samples you have
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 30 -dosnpstat 1 -doHWE 1 -hwe_pval 1e-5 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 77 -snp_pval 1e-6 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2"

echo '#!/bin/bash' > ofavClones.sh
echo angsd -b bamsClones -GL 1 $FILTERS $TODO -P 1 -out ofavClones >>ofavClones.sh

sbatch --mem=200GB -o radClones.o%j -e radClones.e%j --mail-type=ALL --mail-user=studivanms@gmail.com --partition=shortq7 radClones.sh

# scp ibs matrix and bcf file to local machine
cd local/directory/
scp mstudiva@koko-login.hpc.fau.edu:~/project/directory/ANGSD/radClones.ibsMat .
scp mstudiva@koko-login.hpc.fau.edu:~/project/directory/ANGSD/radClones.bcf .


#------------------------------
# Next, run ANGSD_clones.R script, then return here once finished


#------------------------------
# Removing clones

mkdir clones
mv ofavClones* clones

ls *.bam > bamsNoClones

cat bamsClones | grep -v 'OFAV_urban_359.trim.un.bt2.bam\|OFAV_nmfs_035.trim.un.bt2.bam\|OFAV_nmfs_069.trim.un.bt2.bam\|OFAV_nmfs_041.trim.un.bt2.bam\|OFAV_urban_080_2.trim.un.bt2.bam\|OFAV_urban_080_3.trim.un.bt2.bam\|OFAV_urban_289_2.trim.un.bt2.bam\|OFAV_urban_081.trim.un.bt2.bam\|OFAV_urban_415.trim.un.bt2.bam' >bamsNoClones

module load angsd-0.933-gcc-9.2.0-65d64pp
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 30 -dosnpstat 1 -doHWE 1 -hwe_pval 1e-5 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 70 -snp_pval 1e-6 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2"

echo '#!/bin/bash' > ofavNoClones.sh
echo angsd -b bamsNoClones -GL 1 $FILTERS $TODO -P 1 -out ofavNoClones >> ofavNoClones.sh
sbatch --mem=200GB -o radNoClones.o%j -e radNoClones.e%j -p shortq7 --mail-type=ALL --mail-user=studivanms@gmail.com radNoClones.sh

# How many SNPs are there?
grep "filtering:" ofavNoClones.e*
# Number of sites retained after filtering: 20323

# some housekeeping
zipper.py -a -9 -f bam --launcher -e studivanms@gmail.com
sbatch zip.slurm


#------------------------------
## Population structure (pcangsd)

mkdir ../pcangsd
cd pcangsd

cp ../ANGSD/ofavNoClones.beagle.gz .

source activate pcangsd

echo '#!/bin/bash' > pcangsd.sh
echo 'pcangsd -b ofavNoClones.beagle.gz -o ofavPcangsd --admix --inbreed-samples --pcadapt --selection' >> pcangsd
launcher_creator.py -j pcangsd -n pcangsd -q shortq7 -t 06:00:00 -e $EMAIL -w 1 -N 1
sbatch pcangsd.slurm

# Calculate population structure from genotype likelihoods using NGSadmix (part of ANGSD)
mkdir ../ngsAdmix
cp *beagle* ../ngsAdmix

cd ../ngsAdmix

# Simulating (50x) different values of K (minimum of 1, maximum of the number of your sites + 3)
ngsAdmixLauncher.py -f ofavNoClones.beagle.gz --maxK 16 -r 50 -n ofav --launcher -e studivanms@gmail.com

sbatch --mem=200GB ofavNgsAdmix.slurm

# Now determining the most likely value of K using Clumpak with the Evanno method
>ofavNgsAdmixLogfile
for log in ofav*.log; do
grep -Po 'like=\K[^ ]+' $log >> ofavNgsAdmixLogfile;
done

# Start R in terminal
R

logs <- as.data.frame(read.table("ofavNgsAdmixLogfile"))

# output is organized with 10, 11 preceding 1, 2, 3 etc.
logs$K <- c(rep("10", 50), rep("11", 50), rep("12", 50), rep("13", 50), rep("14", 50), rep("15", 50), rep("16", 50), rep("1", 50), rep("2", 50), rep("3", 50), rep("4", 50), rep("5", 50), rep("6", 50), rep("7", 50), rep("8", 50), rep("9", 50))
write.table(logs[, c(2, 1)], "ofavNgsAdmixLogfile_formatted", row.names = F, col.names = F, quote = F)
quit()
# No need to save workspace image [press 'n']
n

# Check that the formatted log file has the correct number of entries
cat ofavNgsAdmixLogfile_formatted | wc -l
# 800, so yes (K=16 * 50 simulations each)

# Make copies of .qopt files to run structure selector on (.Q files)
for file in ofav*.qopt; do
filename=$(basename -- "$file" .qopt);
cp "$file" "$filename".Q;
done

mkdir ofavQ
mv ofav*Q ofavQ

zip -r ofavQ.zip ofavQ

# scp .zip and formatted logfile to local machine and upload to CLUMPAK (http://clumpak.tau.ac.il/bestK.html) and structure selector (https://lmme.ac.cn/StructureSelector/index.html)

# scp ofavNoClones* to local machine for further analyses with R
