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
echo '#!/bin/bash' >rawReads
echo readCounts.sh -e .fastq -o resistRaw >>rawReads
sbatch -o rawReads.o%j -e rawReads.e%j rawReads --mail-type=ALL --mail-user=studivanms@gmail.com


#------------------------------
## Trimming and filtering

2bRAD_trim_launch_dedup.pl fastq > trims.sh
launcher_creator.py -j trims.sh -n trims -q shortq7 -t 06:00:00 -e studivanms@gmail.com
sbatch --mem=200GB trims.slurm
