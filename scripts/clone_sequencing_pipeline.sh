#-------------------------------#
# clonal sequencing analysis ####
#-------------------------------#

# outline
# sequenced all 48 clones of Alex's experiment
# sequenced all 48 clones of Elze's experiment (think these were exposed to copper treatment too)
# sequeced the ancestor: Pseudomonas fluorescens SBW25

# check conda is installed
conda --version
conda update conda

# install packages 
conda install --channel "bioconda" fastqc
conda install -c bioconda freebayes fastqc bowtie2 bcftools samtools bwa samblaster vcflib
conda install -c bioconda biopet-vcffilter
conda install -c bioconda ncbi-genome-download

# check that the files are there
ls 
cd alex_pfluorescens_clones
ls
cd Trimmed
ls


# downloaded the trimmed files - insert the ENA code here

#The raw Fastq files are trimmed for the presence of Illumina adapter sequences using Cutadapt version 1.2.1. The option -O 3 was used, so the 3' end of any reads which match the adapter sequence for 3 bp. or more are trimmed.
#The reads are further trimmed using Sickle version 1.200 with a minimum window quality score of 20. Reads shorter than 15 bp. after trimming were removed. If only one of a read pair passed this filter, it is included in the R0 file. The output files from Cutadapt and Sickle are available here.

# screen
screen -ls

# set up working directory
wd=~/change_folder_name

# set other directories - make them if they are not present
trimmed=$wd/trimmed_files
fastqc_output=$wd/fastqc_output

mkdir -p $trimmed
mkdir -p $fastqc_output

# check it works
ls $wd

# how many raw files are there
ls -lR $wd/raw/*/*.fastq.gz | wc -l
# 285 - this is correct as there is an R0 file in each

# run fastqc on all files
# make sure directory for output exists
fastqc -t 6 -o $fastqc_output $wd/raw/*/*R1_001.fastq.gz
fastqc -t 6 -o $fastqc_output $wd/raw/*/*R2_001.fastq.gz

#-------------------------#
# map to reference genome #
#-------------------------#

# download reference genome
mkdir -p $wd/ref
cd $wd/ref

# find the accession number for the assembly you want to download
ncbi-genome-download --dry-run -A GCF_000009225.2 bacteria --format fasta
ncbi-genome-download -A GCF_000009225.2 bacteria --format fasta

ref=$wd/ref/refseq/bacteria/GCF_000009225.2/GCF_000009225.2_ASM922v1_genomic.fna.gz
gunzip $ref
ref=$wd/ref/refseq/bacteria/GCF_000009225.2/GCF_000009225.2_ASM922v1_genomic.fna

# prepare the reference for mapping
bwa index $ref
samtools faidx $ref

mkdir -p $wd/processed/bams

ls $wd/processed

# file=/home2/ISAD/dp415/alex_pfluorescens_clones/raw/Sample_1/1_TTACCGAC-CGAATACG_L001_R1_001.fastq.gz

# map reads to reference genome
for file in $wd/raw/*/*R1_001.fastq.gz

    do
        echo $file

        # assign fwd and rev files
        file_fwd=$file
        file_rev=${file%R1_001.fastq.gz}R2_001.fastq.gz

        #stub
        stub=$(basename ${file%_R1_001.fastq.gz})

        # skip if file already exists
        if [ -e $wd/processed/bams/${stub}_bwa.bam ]
        then
                continue
        fi


        # do the mapping for fwd and reverse pairs
        # use samblaster to mark duplicates
        
        # method 2: bwa mem
        bwa mem -t 12 $ref $file_fwd $file_rev | samblaster | samtools view -S -b -o $wd/processed/bams/${stub}_bwa.bam

done

# make directories for mapped and unmapped files
bam_files=$wd/processed/bams

mkdir -p $bam_files/mapped
mkdir -p $bam_files/unmapped
mkdir -p $bam_files/stats

# bam_file=/home2/ISAD/dp415/alex_pfluorescens_clones/processed/bams/16_AGACCTTG-TTACGTGC_L001_bwa.bam

for bam_file in $wd/processed/bams/*.bam

    do
        #echo $bam_file
        stub=${bam_file%.bam}
        stub2=$(basename ${bam_file%.bam})

         # skip if file already exists
        if [ -e $wd/processed/bams/mapped/${stub2}_sorted.bam ]
        then
                continue
        fi

        # 1. split bam into mapped and unmapped reads
        samtools view -b -F4 $bam_file > $bam_files/mapped/${stub2}_mapped.bam
        samtools view -f4 $bam_file > $bam_files/unmapped/${stub2}_unmapped.bam

        # 2. sort mapped file by position in genome and not by order of mapped sequence
        samtools sort -o $bam_files/mapped/${stub2}_sorted.bam $bam_files/mapped/${stub2}_mapped.bam

        # 3. index the sorted bam file
        samtools index $bam_files/mapped/${stub2}_sorted.bam

        # 3. remove intermediate (unsorted mapped bam file)
        rm $bam_files/mapped/${stub2}_mapped.bam

        # 4. extract stats of file
        # sorted and mapped bam file
        samtools flagstat $bam_files/mapped/${stub2}_sorted.bam > $wd/processed/bams/stats/${stub2}_mapped_stats.txt
        
        # raw bam file
        samtools flagstat $bam_file > $wd/processed/bams/stats/${stub2}_raw_stats.txt

done

#------------------------------------------#
# run various snp callers on each bam file #
#------------------------------------------#

# make directory for output
mkdir -p $wd/processed/vcf_output/freebayes

# select a single file for testing
# bam_file=$wd/processed/bams/mapped/ancestors_L001_bwa_pa14_sorted.bam

# select files

cd $wd/processed/bams/mapped
find . -name '*sorted.bam' | wc -l

bam_file=$wd/processed/bams/mapped/2_AGTGACCT-TCTAGGAG_L001_bwa_sorted.bam

for bam_file in $wd/processed/bams/mapped/*sorted.bam

    do
        echo $bam_file
        stub=${bam_file%.bam}
        stub2=$(basename ${bam_file%.bam})
        
        # skip if file already exists
        if [ -e $wd/processed/vcf_output/freebayes/${stub2}_freebayes.vcf ]
        then
                continue
        fi

        # call SNPs using FreeBayes, set ploidy to 1 (individual clone sampling)
        freebayes -f $ref -p 1 $bam_file > $wd/processed/vcf_output/freebayes/${stub2}_freebayes.vcf
        # keep only snps with a QUAL score > 20, based on recommendation of Erik Garrison
        #https://www.biostars.org/p/71299/
        vcffilter -f "QUAL > 20" $wd/processed/vcf_output/freebayes/${stub2}_freebayes.vcf > $wd/processed/vcf_output/freebayes/${stub2}_freebayes_filter.vcf
done

cd $wd/processed/vcf_output/freebayes
find . -name '*freebayes.vcf' | wc -l