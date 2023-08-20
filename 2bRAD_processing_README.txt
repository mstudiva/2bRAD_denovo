## 2bRAD reads processing pipeline, version August 15, 2023
# Created by Ryan Eckert (reckert2017@fau.edu), modified by Michael Studivan (studivanms@gmail.com) for this project

https://ryaneckert.github.io/Stephanocoenia_FKNMS_PopGen/data/


#------------------------------
## Modules

# Add to your .bashrc
module load angsd-0.933-gcc-9.2.0-65d64pp
module load bayescan-2.1-gcc-8.3.0-7gakqmd
module load qt-5.15.2-gcc-9.2.0-zi7wcem BayeScEnv/1.1
module load bcftools-1.9-gcc-8.3.0-il4d373
module load bowtie2-2.3.5.1-gcc-8.3.0-63cvhw5
module load cdhit-4.8.1-gcc-8.3.0-bcay75d
module load htslib-1.9-gcc-8.3.0-jn7ehrc
module load kraken2-2.1.1-gcc-9.2.0-ocivj3u
module load python-3.7.4-gcc-8.3.0-3tniqr5
module load launcher
module load miniconda3-4.6.14-gcc-8.3.0-eenl5dj
module load ncbi-toolkit-22_0_0-gcc-9.2.0-jjhd2wa
module load ngsadmix-32-gcc-8.3.0-qbnwmpq
module load ngsRelate/v2
module load R/3.6.1
module load samtools-1.10-gcc-8.3.0-khgksad
module load vcftools-0.1.14-gcc-8.3.0-safy5vc


#------------------------------
## Downloading scripts

cd ~/bin
svn checkout https://github.com/RyanEckert/Stephanocoenia_FKNMS_PopGen/trunk/scripts/
mv scripts/* .
rm -rf scripts

wget http://www.cmpg.unibe.ch/software/PGDSpider/PGDSpider_2.0.7.1.zip
unzip PGDSpider_2.0.7.1.zip
rm PGDSpider_2.0.7.1.zip

# Makes all bin scripts executable
chmod +x *.sh *.pl *.py


#------------------------------
## Downloading reads

echo '#!/bin/bash' > downloadReads.sh
echo 'bs download project --concurrency=high -q -n JA23220 -o .' >> downloadReads.sh
# -n is the project name and -o is the output directory

echo "find . -name '*.gz' -exec mv {} . \;" >> downloadReads.sh
echo 'rmdir SA*' >>downloadReads.sh
#echo 'mkdir ../concatReads' >> downloadReads.sh
#echo 'cp *.gz ../concatReads' >> downloadReads.sh
#echo 'cd ../concatReads' >> downloadReads.sh
echo 'mergeReads.sh -o mergeTemp' >> downloadReads.sh
# -o is the directory to put output files in

echo 'rm *L00*' >> downloadReads.sh
echo "find . -name '*.gz' -exec mv {} . \;" >> downloadReads.sh
echo 'gunzip *.gz' >> downloadReads.sh
echo 'rmdir mergeTemp' >> downloadReads.sh

chmod +x downloadReads.sh

launcher_creator.py -b 'srun downloadReads.sh' -n downloadReads -q shortq7 -t 06:00:00 -e studivanms@gmail.com
sbatch downloadReads.slurm

# Count raw reads
echo '#!/bin/bash' >rawReads.sh
echo readCounts.sh -e .fastq -o resistRaw >>rawReads.sh
sbatch -o rawReads.o%j -e rawReads.e%j rawReads.sh --mail-type=ALL --mail-user=studivanms@gmail.com


#------------------------------
## Trimming and filtering

# Deduplicates row pools into separate 3ill-BC's (1-12), using reverse complement as the ID
2bRAD_trim_launch_dedup.pl fastq > trims.sh
launcher_creator.py -j trims.sh -n trims -q shortq7 -t 06:00:00 -e studivanms@gmail.com
sbatch --mem=200GB trims.slurm

# Do we have the correct number of files?
ls -l *.tr0 | wc -l

mkdir trimmedReads
srun mv *.tr0 trimmedReads/ &

# Rezips the raw fastq's for storage
zipper.py -f fastq -a -9 --launcher -e studivanms@gmail.com
sbatch --mem=200GB zip.slurm

cd ../trimmedReads

# Renames files based on two column lookup table (sampleID.csv): filename, then sample ID
srun sampleRename.py -i sampleID -f tr0

# Creating conda environment for cutadapt, since it conflicts with launcher_creator
# module load miniconda3-4.6.14-gcc-8.3.0-eenl5dj
# conda config --add channels defaults
# conda config --add channels bioconda
# conda config --add channels conda-forge
# conda create -n cutadaptenv cutadapt

conda activate cutadaptenv

# For loop to generate a list of commands for each file
echo '#!/bin/bash' > trimse.sh
echo 'module load miniconda3-4.6.14-gcc-8.3.0-eenl5dj' >> trimse.sh
echo 'conda activate cutadaptenv' >> trimse.sh
for file in *.tr0; do
echo "cutadapt -q 15,15 -m 36 -o ${file/.tr0/}.trim $file > ${file/.tr0/}.trimlog.txt" >> trimse.sh;
done

# Old-school way of submitting jobs
sbatch -o trimse.o%j -e trimse.e%j --mem=200GB trimse.sh
sbatch -o trimse2.o%j -e trimse2.e%j --mem=200GB trimse2.sh

conda deactivate

# Do we have the correct number of files?
ls -l *.trim | wc -l

# Counting the trimmed reads
echo '#!/bin/bash' >cleanReads
echo readCounts.sh -e trim -o Filt >>cleanReads
sbatch --mem=200GB cleanReads

mkdir ../filteredReads
mv *.trim ../filteredReads

# Rezips the row pools for storage
zipper.py -f tr0 -a -9 --launcher -e studivanms@gmail.com
sbatch zip.slurm

cat FiltReadCounts


#------------------------------
## Aligning to reference genome

cd ~/db/
GENOME_FASTA=Ofaveolata.fasta

echo '#!/bin/bash' >genomeBuild.sh
echo bowtie2-build $GENOME_FASTA $GENOME_FASTA >>genomeBuild.sh
echo samtools faidx $GENOME_FASTA >>genomeBuild.sh

sbatch -o genomeBuild.o%j -e genomeBuild.e%j --mem=200GB genomeBuild.sh
