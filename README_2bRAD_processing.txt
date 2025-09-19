## 2bRAD reads processing pipeline, version September 10, 2025
# Based on 2bRAD pipeline by Ryan Eckert (reckert2017@fau.edu)
https://ryaneckert.github.io/Stephanocoenia_FKNMS_PopGen/code/
# Modified by Michael Studivan (studivanms@gmail.com)

# NOTE: This protocol is only for use with an established coral genome, in this case, Orbicella faveolata (NCBI GCA_042242905.1).
# NOTE: For de novo "fake" genome assembly for 2bRAD, follow Ryan's protocol link above


#------------------------------
## Loading modules

# Add the following to ~/.bashrc using nano .bashrc
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

# Ctrl + x to save/quit


#------------------------------
## Creating conda environments for specialized modules

# uncomment and run below if you don't have conda set up
# module load anaconda3-2021.05-gcc-9.4.0-llhdqho
# conda config --add channels defaults
# conda config --add channels bioconda
# conda config --add channels conda-forge

# these packages don't play well with KoKo modules, or they require specific dependencies
conda create -n cutadapt cutadapt

conda create -n pcangsd bioconda::pcangsd

conda create -n moments bioconda::moments

conda create -n pgdspider bioconda::pgdspider


#------------------------------
## Downloading scripts

cd ~/bin
git clone https://github.com/RyanEckert/Stephanocoenia_FKNMS_PopGen
mv Stephanocoenia_FKNMS_PopGen/scripts/* .

git clone https://github.com/mstudiva/2bRAD_denovo

git clone https://github.com/xiaoming-liu/stairway-plot-v2.git
mv stairway-plot-v2.git/trunk/stairway_plot_v2.1.2.zip .
unzip stairway_plot_v2.1.2.zip
rm -r stairway_plot_v2.1.2.zip stairway-plot-v2.git

# Makes all bin scripts executable
chmod +x *.sh *.pl *.py

# re-launch .bashrc
source .bashrc


#------------------------------
## Downloading raw reads

mkdir rawReads
cd rawReads

# transfer directly from Illumina BaseSpace)
echo '#!/bin/bash' > downloadReads.sh
echo 'bs download project --concurrency=high -q -n ######## -o .' >> downloadReads.sh
# -n is the project name and -o is the output directory

echo "find . -name '*.gz' -exec mv {} . \;" >> downloadReads.sh
echo 'rmdir SA*' >>downloadReads.sh
echo 'mkdir ../concatReads' >> downloadReads.sh
echo 'cp *.gz ../concatReads' >> downloadReads.sh
echo 'cd ../concatReads' >> downloadReads.sh
echo 'mergeReads.sh -o mergeTemp' >> downloadReads.sh
# -o is the directory to put output files in

echo 'rm *L00*' >> downloadReads.sh
echo "find . -name '*.gz' -exec mv {} . \;" >> downloadReads.sh
echo 'gunzip *.gz' >> downloadReads.sh
echo 'rmdir mergeTemp' >> downloadReads.sh

chmod +x downloadReads.sh

launcher_creator.py -b 'srun downloadReads.sh' -n downloadReads -q shortq7 -t 06:00:00 -e studivanms@gmail.com
sbatch --mem=200GB downloadReads.slurm

# Count raw reads
echo '#!/bin/bash' >rawReads.sh
echo readCounts.sh -e gz -o Raw >>rawReads.sh
sbatch -o rawReads.o%j -e rawReads.e%j --mail-type=ALL --mail-user=studivanms@gmail.com rawReads.sh
# scp RawReadCounts to local machine


#------------------------------
## Deduplication and trimming

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

# If any files are not renamed, their barcodes have been mis-sequenced
# move them to a different directory (the hyphen is only found in pre-renamed files)
mkdir unused
mv *-* unused/

# For loop to generate a list of commands for each file
echo '#!/bin/bash' > trimse.sh
echo 'module load miniconda3-4.6.14-gcc-8.3.0-eenl5dj' >> trimse.sh
echo 'conda activate cutadapt' >> trimse.sh
for file in *.tr0; do
echo "cutadapt -q 15,15 -m 36 -o ${file/.tr0/}.trim $file > ${file/.tr0/}.trimlog.txt" >> trimse.sh;
done

# Since this job cannot be run in parallel, split job script up and run each separately
# Creating one script per species
for i in {1..6}; do cp trimse.sh "trimse$i.sh"; done
conda activate cutadapt
sbatch -o trimse.o%j -e trimse.e%j --mem=200GB --mail-type=ALL --mail-user=studivanms@gmail.com trimse.sh
sbatch -o trimse2.o%j -e trimse2.e%j --mem=200GB --mail-type=ALL --mail-user=studivanms@gmail.com trimse2.sh
sbatch -o trimse3.o%j -e trimse3.e%j --mem=200GB --mail-type=ALL --mail-user=studivanms@gmail.com trimse3.sh
sbatch -o trimse4.o%j -e trimse4.e%j --mem=200GB --mail-type=ALL --mail-user=studivanms@gmail.com trimse4.sh
sbatch -o trimse5.o%j -e trimse5.e%j --mem=200GB --mail-type=ALL --mail-user=studivanms@gmail.com trimse5.sh
sbatch -o trimse6.o%j -e trimse6.e%j --mem=200GB --mail-type=ALL --mail-user=studivanms@gmail.com trimse6.sh
sbatch -o trimse7.o%j -e trimse7.e%j --mem=200GB --mail-type=ALL --mail-user=studivanms@gmail.com trimse7.sh
conda deactivate

# Do we have the correct number of files?
ls -l *.trim | wc -l

# Counting the trimmed reads
echo '#!/bin/bash' >cleanReads
echo readCounts.sh -e trim -o Filt >>cleanReads
sbatch --mem=200GB --mail-type=ALL --mail-user=studivanms@gmail.com cleanReads
# scp FiltReadCounts to local machine

mkdir ../filteredReads
mv *.trim ../filteredReads

# Rezips the row pools for storage
zipper.py -f tr0 -a -9 --launcher -e studivanms@gmail.com
sbatch zip.slurm

cat FiltReadCounts


#------------------------------
## Genome formatting

cd ~/db/

module load bowtie2-2.3.5.1-gcc-8.3.0-63cvhw5
module load samtools-1.10-gcc-8.3.0-khgksad

mkdir ~/db/symGenomes

# Concatenated symbiont genomes
# Using concatenated Symbiodiniaceae references from NCBI's most recent genomes:
# Symbiodinium (NCBI GCA_965279495.1), Breviolum (GCA_965643015.1), Cladocopium (GCA_947184155.2), Durusdinium (GCA_963970005.1)

python concatFasta.py -o symbConcatGenome.fasta -s symbConcatGenome_summary.tsv -m symbConcatGenome_contig_to_fakechr.tsv GCA_965279495.1_pySymTrid1.1_genomic.fna GCA_965643015.1_pyBreMinu3.1_genomic.fna GCA_947184155.2_Cgoreaui_SCF055-01_v2.1_genomic.fna GCA_963970005.1_Durusdinium_trenchii_CCMP2556_genomic.fna

# Now building bowtie2 index for concatenated symbiont genomes
echo '#!/bin/bash' >genomeBuild.sh
echo bowtie2-build symbConcatGenome.fasta symbConcatGenome >>genomeBuild.sh
echo samtools faidx symbConcatGenome.fasta >>genomeBuild.sh
sbatch -o genomeBuild.o%j -e genomeBuild.e%j --mail-type=ALL --mail-user=studivanms@gmail.com genomeBuild.sh

# Building bowtie2 index for coral genome
echo '#!/bin/bash' >genomeBuild.sh
echo bowtie2-build Orbicella_faveolata_gen_17.scaffolds.fa OfaveolataGenome >>genomeBuild.sh
echo samtools faidx Orbicella_faveolata_gen_17.scaffolds.fa >>genomeBuild.sh
sbatch -o genomeBuild.o%j -e genomeBuild.e%j --mail-type=ALL --mail-user=studivanms@gmail.com genomeBuild.sh


#------------------------------
## Symbiont Alignment

module load bowtie2-2.3.5.1-gcc-8.3.0-63cvhw5

SYMGENOME=~/db/symGenomes/symbConcatGenome

# aligning reads to concatenated symbiont reference
2bRAD_bowtie2_launcher.py -g $SYMGENOME -f trim -n zooxMaps --split -u un -a zoox --launcher -e studivanms@gmail.com
sbatch zooxMaps.slurm

# some housekeeping
mkdir ../mappedReads
mkdir ../mappedReads/symbionts
mv *.sam ../mappedReads/symbionts
mv *.zoox ../mappedReads/symbionts
cd ../mappedReads/symbionts

# Counting the mapped zoox reads
# calculate mapping efficiency from these values compared to trimmed reads in Excel
echo '#!/bin/bash' >mappedZooxReads
echo readCounts.sh -e zoox -o Zoox >>mappedZooxReads
sbatch --mem=200GB --mail-type=ALL --mail-user=studivanms@gmail.com mappedZooxReads
# scp ZooxReadCounts to local machine

module load samtools-1.10-gcc-8.3.0-khgksad

# making script to generate indexed bam files
>s2b
for file in *.sam; do
echo "samtools sort -O bam -o ${file/.sam/}.bam $file && samtools index ${file/.sam/}.bam">>s2b;
done
launcher_creator.py -j s2b -n s2b -t 6:00:00 -N 5 -e studivanms@gmail.com -q shortq7
sbatch s2b.slurm

# counting the symbiont reads by genera
>ZooxReads
for i in *.bam; do
echo $i >>ZooxReads;
samtools idxstats $i | cut -f 1,3 >>ZooxReads;
done

# some more housekeeping
zipper.py -a -9 -f sam --launcher -e studivanms@gmail.com
sbatch zip.slurm

zipper.py -a -9 -f zoox --launcher -e studivanms@gmail.com
sbatch zip.slurm


#------------------------------
## Host alignment (2bRAD)

cd ~/project/directory/2bRAD/filteredReads
# if your samples are gzipped:
zipper.py -a -9 -f gz --gunzip --launcher -e studivanms@gmail.com
sbatch zip.slurm

HOSTGENOME=~/db/ofavGenome/OfaveolataGenome

mkdir junk

# mapping with --local option, enables clipping of mismatching ends (guards against deletions near ends of RAD tags)
2bRAD_bowtie2_launcher.py -g $HOSTGENOME -f un --split -u junk -a host --undir junk --launcher -e studivanms@gmail.com
sbatch --mem=200GB maps.slurm

ls *.sam | wc -l

echo '#!/bin/bash' >mappedReads
echo readCounts.sh -e al -o Host >>mappedReads
sbatch --mem=200GB --mail-type=ALL --mail-user=studivanms@gmail.com mappedReads
# scp HostReadCounts to local machine

>s2b
for file in *.sam; do
echo "samtools sort -O bam -o ${file/.sam/}.bam $file && samtools index ${file/.sam/}.bam">>s2b;
done

launcher_creator.py -j s2b -n s2b -q shortq7 -t 06:00:00 -e studivanms@gmail.com
sbatch --mem=200GB s2b.slurm

ls *bam | wc -l

# some housekeeping
zipper.py -a -9 -f sam --launcher -e studivanms@gmail.com
sbatch zip.slurm

zipper.py -a -9 -f un --launcher -e studivanms@gmail.com
sbatch zip.slurm

zipper.py -a -9 -f trim --launcher -e studivanms@gmail.com
sbatch zip.slurm

zipper.py -a -9 -f host --launcher -e studivanms@gmail.com
sbatch zip.slurm

mv *.trim.gz ../../trimmedReads
mv *.sam.gz ../mappedReads
mv *.host.gz ../mappedReads
mv *.bam ../mappedReads
mv *.bai ../mappedReads

cd junk
zipper.py -a -9 -f junk --launcher -e studivanms@gmail.com
sbatch zip.slurm


#------------------------------
# Now proceed with README_ANGSD_processing
