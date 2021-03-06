# laevis_GBS_2020

# Background
Received aligned and indexed data from Ben(.bam and .bai) Going to craete structure plots and cal Fst)
Copied those files into 
** all the related scripts can be found in /scratch/premacht/laevis_GBS_2020/saved_scripts_to_be_used **

# Filter and prepare clear bam files (somehow sort at the end had not exactly done the job. So I had to sort and index again in next script)
Run this script in the folder with bam files

```bash
#!/bin/sh
#SBATCH --job-name=bwa_505
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=512gb
#SBATCH --output=bwa505.%J.out
#SBATCH --error=bwa505.%J.err
#SBATCH --account=def-ben

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

module load bwa
module load samtools/1.10
module load nixpkgs/16.09
module load intel/2018.3
module load bcftools/1.10.2
module load vcftools/0.1.16

for i in *.bam; do

# This part selects only chrs and saves an inbetween file

samtools view -b ${i} chr1L chr1S chr2L chr2S chr3L chr3S chr4L chr4S chr5L chr5S chr6L chr6S chr7L chr7S chr8L chr8S chr9_10L chr9_10S > inbetween.bam

#index inbetween bam
samtools index inbetween.bam

# Remove indels
samtools view -h inbetween.bam | awk '$1 ~ "^@" || $6 !~ "I|D"' | samtools view -b > test.bam

# sort for final output
samtools sort test.bam -o ${i}_final_sorted.bam

# Index the final output file

samtools index ${i}_final_sorted.bam;done

# Seperate L subgenome
for i in *_final_sorted.bam;do samtools view -b $i chr1L chr2L chr3L chr4L chr5L chr6L chr7L chr8L chr9_10L > ${i%%cuttrim_sorted.bam_final_sorted.bam}before_final_l_only.bam; done

# Seperate S subgenome
for i in *_final_sorted.bam;do samtools view -b $i chr1S chr2S chr3S chr4S chr5S chr6S chr7S chr8S chr9_10S > ${i%%cuttrim_sorted.bam_final_sorted.bam}before_final_s_only.bam; done

# Seperate whole genome/ without scaffolds
for i in *_final_sorted.bam;do samtools view -b $i chr1S chr2L chr2S chr3L chr3S chr4L chr4S chr5L chr5S chr6L chr6S chr7L chr7S chr8L chr8S chr9_10L chr9_10S > ${i%%cuttrim_sorted.bam_final_sorted.bam}before_final_whole_genome.bam; done

# Make directories for bamfiles for different genomes
mkdir ../filtered_bam_files

mkdir ../filtered_bam_files/whole_genome
mkdir ../filtered_bam_files/s_only
mkdir ../filtered_bam_files/l_only

# Move files to corresponding directories
mv *whole_genome.ba* ../filtered_bam_files/whole_genome
mv *_l_only.ba* ../filtered_bam_files/l_only
mv *_s_only.ba* ../filtered_bam_files/s_only


```

## in one of the directories created by above script for subgenomes(eg-filtered_bam_files/s_only) run the following script.
This will
 

Run these in the filtered_bam_files/l_only folder to download and prepare reference genome
```bash
# Create a directory for reference genome
mkdir ../reference_genome

# Download reference genome, unzip it and move to the reference genome folder
wget http://ftp.xenbase.org/pub/Genomics/JGI/Xenla9.2/XENLA_9.2_genome.fa.gz
gunzip XENLA_9.2_genome.fa.gz
mv XENLA_9.2_genome.fa* ../reference_genome

module load bwa
module load samtools/1.10
samtools faidx ../reference_genome/XENLA_9.2_genome.fa
```
# sort and index bams again, create vcfs for all subgenomes, remove scaffolds and labels from header of the vcfs and store them in seperate folders created for them (These files will be ready to use for analysis)

# ** Processing all files this way might take lot of time. In that case, modify the outer for loop and use this script to submit jobs parallely for seperate genomes(leaving just one genome as j at a time) * Remember to give the whole genome more time and memory when you do this***

Run in filtered_bam_files/l_only

```bash
#!/bin/sh
#SBATCH --job-name=bwa_505
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --mem=512gb
#SBATCH --output=bwa505.%J.out
#SBATCH --error=bwa505.%J.err
#SBATCH --account=def-ben

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

module load bwa
module load samtools/1.10
module load nixpkgs/16.09
module load intel/2018.3
module load bcftools/1.10.2
module load vcftools/0.1.16


# Create a directory for filtered vcf files
mkdir ../../filtered_VCFs

# loop through all subgenomes

for j in  ../l_only ../s_only ../whole_genome; do cd ${j}

# previous sort had not worked correctly. So I had to sort again here before creating VCF
for i in *_before_final*bam; do samtools sort ${i} -o ${i%*_before_final*bam}_final_${j#../}.bam
        samtools index ${i%*_before_final*bam}_final_${j#../}.bam & done

# move older files to a seperate folder
mkdir bams_before_final
mv *before_final* bams_before_final

# Create VCF
bcftools mpileup -a FORMAT/DP -q20 -d8000 -f ../reference_genome/XENLA_9.2_genome.fa *_final_${j#../}.bam | bcftools call -V indels --format-fields GQ -m -O z -O z -o laevis_GBS_2020_${j#../}.vcf.gz

# unzip VCF
gunzip laevis_GBS_2020_${j#../}.vcf.gz

# Get sample list
vcf-query -l laevis_GBS_2020_${j#../}.vcf > all_sample_list

# Select chromosomes and remove scaffolds, remove scaffold from header (Selected chromosomes here moght not affect subgenomes here as bam files used are already filtered for subgenomes)

vcftools --gzvcf laevis_GBS_2020_${j#../}.vcf --keep all_sample_list --chr chr1L --chr chr1S --chr chr2L --chr chr2S --chr chr3L --chr chr3S --chr chr4L --chr chr4S --chr chr5L --chr chr5S --chr chr6L --chr chr6S --chr chr7L --chr chr7S --chr chr8L --chr chr8S --chr chr9_10L --chr chr9_10S --non-ref-ac-any 1 --recode --recode-INFO-all --stdout | grep -v '^##contig=<ID=Scaffold' > laevis_GBS_2020_${j#../}_scaffolds_removed.vcf


# make a folder for the VCF 
mkdir mkdir ../../filtered_VCFs/vcf_${j#../}

#move vcf there
mv laevis_GBS_2020_${j#../}* ../../filtered_VCFs/vcf_${j#../}/ & done

```

## Run following in one of the folders created for filtered bams by the first script to collect finalized bams in a seperate folder(filtered_bam_files/l_only)

 make a directory to collect all finalized bam files
 ```bash
mkdir ../../finalized_bam_files
```

 Collect all finalized bam files looping through all directories and place them in respective folders
```bash
for j in  ../l_only ../s_only ../whole_genome; do cd ${j}
        mkdir ../../finalized_bam_files/finalized_bams_${j#../}
        cp ./*final_${j#../}.bam* ../../finalized_bam_files/finalized_bams_${j#../} ; done
```
# ====>> You can start population structure analysis with these data(finalized bams from the step above) as described after Fst =======>>>>>

# Filtering based on per site coverage

in the reference genome folder,

** YOU MIGHT HAVE TO LOAD GATK BEFORE SUBMITTING JOB TO TYPE YES AND AGREE TERMS**


create GATk dictionary file for reference genome
```bsh
module load nixpkgs/16.09
module load gatk/4.1.2.0
gatk --java-options "-Xmx2G" CreateSequenceDictionary -R   XENLA_9.2_genome.fa
```
 Create depth table using GATK for all genomes(in /filtered_VCFs/vcf_l_only)
 
 ** You will not be able to get enough memory for java(like -Xmx16G) if you use bash to run script(max is 1/4 th of the physical memory you ask for). Therefore to test, submit job like this and keep checking .err file **
 
 In filtered_VCFs/vcf_l_only folder,
 
 ```bash
#!/bin/sh
#SBATCH --job-name=bwa_505
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --mem=512gb
#SBATCH --output=bwa505.%J.out
#SBATCH --error=bwa505.%J.err
#SBATCH --account=def-ben

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

module load nixpkgs/16.09
module load gatk/3.8

for j in ../*; do cd ${j}; java -Xmx16G -jar $EBROOTGATK/GenomeAnalysisTK.jar -T VariantsToTable -R ../../filtered_bam_files/reference_genome/XENLA_9.2_genome.fa -V *_scaffolds_removed.vcf -F CHROM -F POS -GF DP -o laevis_GBS_2020_${j#../}_GVCF_DP.table ; done
```
# Cal moving average depath, cal cutoffs , plot mean_over_cuttoff and output a file with sites to exclude, ready to be used by VCFtools in the next step (to run this for all genomes together, copied this script in saved_scripts_to_be_used and ran the script following this)

Run this where the DP.Table file is from the previous script

```R
# set working directory to current script directory 
# **uncomment this if you use this in local computer. KEEP COMMENTED OUT IF YOU ARE ON COMPUTECANADA) ***
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library (ggplot2)

# get all the files with the site data
files <- list.files(path = ".", pattern = "GVCF_DP.table")
my_data<-c()
temp<-c()

# read in the data and name the df based on the file name
for(f in 1:length(files)) {
  temp <- read.table(files[f], header = T)
  my_data <- rbind(my_data,temp)
  temp<-c()
}  

# my_data now has coverage per site info for each sample for all sites, genomewide
dim(my_data)
colnamez <- colnames(my_data)

# here is a function to calculate moving averages
# https://stackoverflow.com/questions/743812/calculating-moving-average
moving_fun <- function(x, w, FUN, ...) {
  # x: a double vector
  # w: the length of the window, i.e., the section of the vector selected to apply FUN
  # FUN: a function that takes a vector and return a summarize value, e.g., mean, sum, etc.
  # Given a double type vector apply a FUN over a moving window from left to the right, 
  #    when a window boundary is not a legal section, i.e. lower_bound and i (upper bound) 
  #    are not contained in the length of the vector, return a NA_real_
  if (w < 1) {
    stop("The length of the window 'w' must be greater than 0")
  }
  output <- x
  for (i in 1:length(x)) {
    # plus 1 because the index is inclusive with the upper_bound 'i'
    lower_bound <- i - w + 1
    if (lower_bound < 1) {
      output[i] <- NA_real_
    } else {
      output[i] <- FUN(x[lower_bound:i, ...])
    }
  }
  output
}

# example
# v <- seq(1:10)

# compute a MA(2)
# moving_fun(v, 2, mean)


# make a new dataframe that has the moving average for each sample throughout the 
# genome.  No worries if the window goes across chromosomes
mv_ave_df2 <- c()
for (i in 3:length(my_data)){  
  # calculate the moving average for each column
  moving_average <- moving_fun(my_data[,i], 50, mean)
  # add this to a new dataframe
  mv_ave_df2 <- cbind(mv_ave_df2,moving_average)
  # rename the column to match the sample
  colnames(mv_ave_df2)[i-2] = colnamez[i]
}  

# add chromosome and position data to mv_ave_df2
mv_ave_df3 <- data.frame(mv_ave_df2,my_data$CHROM,my_data$POS)


colnames(mv_ave_df3)[ncol(mv_ave_df3)-1] <- "CHROM"
colnames(mv_ave_df3)[ncol(mv_ave_df3)] <- "POS"


# calculate mean depth and sd per sample
mean_depth<- c()
sd_depth<- c()

for (i in 3:length(my_data)){
  x<- mean(my_data[,i],na.rm=T)
  mean_depth<-append(mean_depth,x, after = length(mean_depth))
  y<-sd(my_data[,i],na.rm=T)
  sd_depth<-append(sd_depth,y, after = length(sd_depth))
}  

cutoff_vector <- c()


# identify chr and positions in any sample that are >4 sd above that samples mean coverage
# based on the rolling average
# first make a vector with cutoff values for each sample
for (i in 3:ncol(mv_ave_df3)){
  cutoff <- mean_depth[i-2] + 4*sd_depth[i-2]
  cutoff_vector <- append(cutoff_vector,cutoff, after=length(cutoff_vector))
}

mean_over_cuttoff <- mean_depth/cutoff_vector
cutoff_data<-data.frame(1:length(mean_over_cuttoff),mean_over_cuttoff)
p<-ggplot(cutoff_data)+
  geom_point(aes(x=X1.length.mean_over_cuttoff.,y=mean_over_cuttoff))+
  theme_bw()

ggsave(filename = gsub("GVCF_DP.table","mean_over_cuttoff_plot.pdf",files),plot = p,width = 30,height = 10)
# give the elements in cutoff_vector some names

#manual
#names(cutoff_vector) <- c('bru_PF707.DP','download.DP')


# automated *******************************************
#get sample names from mv_ave_df3
xx<-names(mv_ave_df3)

# select just sample names
yy<-xx[1:(length(xx)-2)]

# take sample names for vector names
names(cutoff_vector)<-yy
# zz<-sub("_cuttrim_sorted_final_l_only.bam.DP","",yy)

# ****************************************************

# now cycle through each column and identify chr and pos of the bad ones
# manual
#subtest <- NULL
#subtest <- subset(mv_ave_df3, 
#              JM_no_label1_Draken_CCACGT_cuttrim_sorted_final_l_only.bam.DP > cutoff_vector["JM_no_label1_Draken_CCACGT_cuttrim_sorted_final_l_only.bam.DP"]  | 
#                JM_no_label2_Draken_TTCAGA_cuttrim_sorted_final_l_only.bam.DP > cutoff_vector["JM_no_label2_Draken_TTCAGA_cuttrim_sorted_final_l_only.bam.DP"] 
# )

# automated - Does the same thing as above**********
sub <- NULL

sub<-subset(mv_ave_df3,
            get(yy)>min(cutoff_vector),
            )


# **************************************************

# automated - only collects sites with a higher depth than that samples own cutoff*****************************************
#sub <- NULL

#for (i in yy) {
  
#  sub <- subset(mv_ave_df3,
#               get(i) >cutoff_vector[i] 
#  )
#}
# **************************************************

dim(sub)
dim(mv_ave_df3)


to_file <- data.frame(sub$CHROM, sub$POS)

# write to file
write.table(to_file, gsub("GVCF_DP.table","positions_to_exclude.txt",files), append = FALSE, sep = "\t", dec = ".",
            row.names = F, col.names = F,quote = FALSE)
```
# Run the above script for all subgenomes at once 

Copied the script in saved_scripts_to_be_used first.
and then,

Load R

```bash
module load nixpkgs/16.09
module load gcc/7.3.0
module load r
```
to run script in all directories, use the following command in filtered_VCFs/vcf_l_only

```bash
for j in ../vcf*;do
       cd ${j}
        cp /scratch/premacht/laevis_GBS_2020/saved_scripts_to_be_used/cal_moving_dp_and_find_excludes.r .
        Rscript cal_moving_dp_and_find_excludes.r ;done
        
```

# Filter above selected sites using vcftools

```bash
module load nixpkgs/16.09
module load intel/2018.3
module load vcftools/0.1.16

for j in ../vcf*;do
       cd ${j}
vcftools --vcf laevis_GBS_2020_${j#../vcf_}_scaffolds_removed.vcf --out laevis_GBS_2020_${j#../}_positions_excluded --exclude-positions *positions_to_exclude.txt --recode ;done

```

#  Run following in the filtered_VCFs/vcf_l_only to collect finalized VCFs in a new directory to be used in the next steps.

```bash
mkdir ../../finalized_vcf_files

for j in  ../*l_only ../*s_only ../*whole_genome; do cd ${j}
mv ./*positions_excluded.recode.vcf ../../finalized_vcf_files ; done

```

Copy parseVCF.py from https://github.com/simonhmartin/genomics_general/blob/master/VCF_processing/parseVCF.py to finalized_vcf_files making a vi file with the same name(parseVCF.py)
** For some reason, you cannot paste python scripts directly to a file named xxxx.py. Therefore you may have to copy it to a vi file with a different name and then rename it as needed for the python script. **

then,

Load python, create vitual env in home directory, activate it, upgrade pip in the environment and install numpy package(ENV  is the name of the empty directory containing your environment)

```bash
module load python/3.8.2
virtualenv --no-download ~/ENV
source ~/ENV/bin/activate
pip install --no-index --upgrade pip
pip install numpy --no-index
```
Then convert VCF into geno format and move them to a newly created directory
```bash
for j in *positions_excluded.recode.vcf ; do python parseVCF.py -i ${j} -o ${j%%recode.vcf}geno ; done
mkdir ../geno_files
mv *.geno ../geno_files/
```


* When you want to exit environment, use
```bash
(ENV) [name@server ~] deactivate
```

# Calculating Fst with geno files

** needed help for these steps can be found in https://github.com/simonhmartin/genomics_general#diversity-and-divergence-analyses-in-sliding-windows **

following files are needed in the directory with geno files for the script to work

1.copy python script from https://github.com/simonhmartin/genomics_general/blob/master/genomics.py to the directory with geno files.
and
2.copy python script from https://github.com/simonhmartin/genomics_general/blob/master/popgenWindows.py 
as well.
3.Create a populatin file (--popsFile) , which has two columns: the first gives sample names and teh second gives population name:(you can use excel and save the output as txt)

eg:
```txt
JM_no_label1_Draken_CCACGT_cuttrim_sorted_final_l_only.bam      popC
JM_no_label2_Draken_TTCAGA_cuttrim_sorted_final_l_only.bam      popD
```
4. Needed geno files

Then 
# Calculate standard population genomic statistics in sliding windows: pi, FST and DXY using following command 
(explanations about the flags used can be found in the help link above)
(you may include command "--writeFailedWindows" if you get an empty output file, it should always write something in the output file, even if all results are NA)

**use this script to take population names and sample names to paste in the next command

```r
#set current path as wd
library("rstudioapi") 
setwd(dirname(getActiveDocumentContext()$path))

pop_list<-read.table("./pop_file.txt")


final_list_to_paste<-"paste following in the terminal:  "

for (pop in unique(pop_list$V2)) {

pop_1<-pop_list$V1[pop_list$V2==pop]

sample_list<-noquote(toString(pop_1))

full_cmd<-paste(" -p ",pop," ",sample_list,sep = "")

final_list_for_this_pop<-gsub(", ",",",full_cmd)

final_list_to_paste<-paste(final_list_to_paste,final_list_for_this_pop,sep = "")
}

print(final_list_to_paste)

```

```bash
python popgenWindows.py -w 50000 -m 50 -g laevis_GBS_2020_vcf_l_only_positions_excluded.geno -o output.csv.gz -f phased --popsFile pops.txt
```

# ********* filter genes ********************

Generate the needed gff files with data needed for filtering like following

```curl
curl http://ftp.xenbase.org/pub/Genomics/JGI/Xenla9.2/XENLA_9.2_Xenbase.gff3 > XLv9.2_xenbase_annotations.gff
grep 'exon' XLv9.2_xenbase_annotations.gff | egrep -v 'mRNA' > XLv9.2_xenbase_annotations_exons_only.gff
``` 
Then filter the genes with the help of that file with the following script by running it as shown below

Prints_Fst_blocks_with_genes.pl
```perl
#!/usr/bin/env perl
use strict;
use warnings;


# This program reads in the output the the general genomics popgenWindows.py script
# and then a gff file that has only the mRNA coordinates.
# It then prints out lines from the general genomics file that have genes.



# first concatenate FST data from all chrs like this:
# ../Makes_inputfile_for_jackknife.pl prefix

# download the xenbase file:
# curl http://ftp.xenbase.org/pub/Genomics/JGI/Xenla9.2/XENLA_9.2_Xenbase.gff3 > XLv9.2_xenbase_annotations.gff

# make a simplified gff file that has only the mRNA coordinates
# grep 'mRNA' XLv9.2_xenbase_annotations.gff | egrep -v 'exon' > XLv9.2_xenbase_annotations_mRNA_only.gff

# make the perl script executable:
# chmod 755 Prints_Fst_blocks_with_genes.pl

# then run this file like this (with four variables following the perl script):
# ./Prints_Fst_blocks_with_genes.pl laevis_GBS_2020_l_only_positions_excluded_final_fst.csv XLv9.2_xenbase_annotations_mRNA_only.gff output_with_only_blocks_containing_genes.out 

# laevis_GBS_2020_l_only_positions_excluded_final_fst.csv is the output file from the Fst analysis with results from all chromosomes concatenated and only one header in the file
# XLv9.2_xenbase_annotations_mRNA_only.gff is the gff file with only mRNA coordinates
# output_with_only_blocks_containing_genes.out is the output file with Fst results from only those blocks containing genes 


my $windowsize=10000; # this is the size of the genomic window 


my $inputfile = $ARGV[0]; # the name of the Fst input file is from the commandline
unless (open DATAINPUT, $inputfile) {
        print "Can not find the Fst input file.\n";
        exit;
}
my @temp;
my @temp1;
my %fst;
my @fst;
my $firstline;

# Read in the Fst file from general genomics and store each block in a 
# hash with the key being the chr and start coordinates separated by an underscore
while ( my $line = <DATAINPUT>) {
        @temp=split(',',$line);
        if($temp[0] ne 'scaffold'){
                $fst{$temp[0]."_".$temp[1]}{"line"} = $line; # this line will be printed to the outfile if the next flag is equal to 1
                $fst{$temp[0]."_".$temp[1]}{"chr"} = $temp[0]; # this is a flag that will be updated if we want to print this line
                $fst{$temp[0]."_".$temp[1]}{"start"} = $temp[1]; # this is a flag that will be updated if we want to print this line
                $fst{$temp[0]."_".$temp[1]}{"stop"} = $temp[2]; # this is a flag that will be updated if we want to print this line
                $fst{$temp[0]."_".$temp[1]}{"print_or_not"} = 0; # this is a flag that will be updated if we want to print this line
        }
        else{
                $firstline=$line;
        }
}
close DATAINPUT;

# OK now check to see if genes are in these windows
my $inputfile1 = $ARGV[1]; # the name of the gff file is from the commandline
unless (open DATAINPUT2, $inputfile1) {
        print "Can not find the gff file.\n";
        exit;
}
my $start;
my $value;
# go through each line in the gff file
while ( my $line = <DATAINPUT2>) {
        @temp=split('\t',$line);
        foreach my $key (keys %fst){
        if (($fst{$key}{"chr"} eq $temp[0])&&
                (
                (($fst{$key}{"start"} <= $temp[3])&&($fst{$key}{"start"} + $windowsize >= $temp[3])) # the start of the gene is in the window
                ||
                (($fst{$key}{"stop"} - $windowsize <= $temp[4])&&($fst{$key}{"stop"} >= $temp[4])) # the end of the gene is in the window
                ))
                {
                        # gene is in this block, so print it
                        $fst{$key}{"print_or_not"} = 1;
        }
        }
}
close DATAINPUT2;

my $outputfile = $ARGV[2]; # the name of the output file is from the commandline
unless (open(OUTFILE, ">$outputfile"))  {
        print "I can\'t write to $outputfile\n";
        exit;
}

print OUTFILE $firstline;
foreach my $key (sort keys %fst){
    if ($fst{$key}{"print_or_not"} == 1) {
             print OUTFILE $fst{$key}{"line"};
    }
}

close OUTFILE;

```
Run the script like following

run_filter_gene_blocks.sh
```bash
#!/bin/sh
#SBATCH --job-name=bwa_505
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --output=bwa505.%J.out
#SBATCH --error=bwa505.%J.err
#SBATCH --account=def-ben

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

./Prints_Fst_blocks_with_genes.pl laevis_GBS_2020_l_only_positions_excluded_5k_final_fst.csv XLv9.2_xenbase_annotations_exons_only.gff window5k_filtered_exons_l_only.csv

./Prints_Fst_blocks_with_genes.pl laevis_GBS_2020_s_only_positions_excluded_5k_final_fst.csv XLv9.2_xenbase_annotations_exons_only.gff window5k_filtered_exons_s_only.csv

```

Then use jacknife and prepare the output for plotting

test_block_jackknife_withFst.pl

```perl
#!/usr/bin/env perl
use strict;
use warnings;


# This program reads in the output of the script called 
# and calculates the standard error of the weighted mean value of fDM with weightings based 
# on the sum of the number of abba and baba sites in each window.

my $inputfile = $ARGV[0];

# print headers first

print "pop1","\t","pop2","\t","Fst","\t","CI_1","\t","CI2","\n";

#loop for selected columns
my @selected_chrs = (84..149);

for my $j (@selected_chrs){

#take column no

my $chr = $j;

unless (open DATAINPUT, $inputfile) {
        print "Can not find the input file.\n";
        exit;
}

my @windowsites;
my @Fst_values;
my $sumsites=0;
my $counter=0;
my @temp;
my $y;
my $x;

while ( my $line = <DATAINPUT>) {
        @temp=split(',',$line);
        if($temp[0] ne 'scaffold'){
                # if fDM is not numeric, it shouldn't contribute to the weighted average
                if($temp[$chr] !~ /nan/){
                        # load the number of sites for each window
                        $windowsites[$counter]=$temp[4];
                        # also load the Fst for each window
                        $Fst_values[$counter]=$temp[$chr];
                        # sum all the sites
                        $sumsites+= $temp[4];
                        # keep track of how many windows are loaded
                        $counter+=1;
                }
        }
}

my $weighted_Fst_average=0;

# calculate the weighted average of Fst
for ($y = 0 ; $y <= $#Fst_values; $y++ ) {
        $weighted_Fst_average+=($windowsites[$y]*$Fst_values[$y]/$sumsites);
}

# Print it
# print "The weighted Fst average is ",sprintf("%.6f",$weighted_Fst_average),"\n";


# now calculate the standard error.
my @jack_sumsites;
my $jack_averagesites;
my @jack_fst_values;
my $jack_weighted_average=0;
my @jack_weighted_fdm_values;
my @jackarray;
my $counter2=0;

for ($y = 0 ; $y < $counter; $y++ ) {
        for ($x = 0 ; $x < $counter; $x++ ) {
                # leave out one row for each jackknfe replicates
                if($y != $x){
                        # load the number of sites for each window
                        $jack_sumsites[$counter2]=$windowsites[$x];
                        # add them up to eventually get the average
                        $jack_averagesites+=$windowsites[$x];
                        # also load the Fst stat for each window
                        $jack_fst_values[$counter2]=$Fst_values[$x];
                        # keep track of the number of windows
                        $counter2+=1;
                }
        }
        # make the $jack_averagesites an average
        $jack_averagesites=$jack_averagesites/($counter2);
        # calculate the weighted average
        for ($x = 0 ; $x < $counter2; $x++ ) {
                $jack_weighted_average+=($jack_sumsites[$x]*$jack_fst_values[$x]/$jack_averagesites);
        }
        push(@jackarray,($jack_weighted_average/$counter2));
        # reset variables
        $jack_weighted_average=0;
        $jack_averagesites=0;
        @jack_fst_values=();
        $counter2=0;
}

# now calculate the variance of the jackknife replicates

# first we need the mean
my $jack_mean=0;
for ($x = 0 ; $x < $counter; $x++ ) {
        $jack_mean+=$jackarray[$x];
}
$jack_mean=$jack_mean/($counter);

my $jack_var=0;
for ($x = 0 ; $x < $counter; $x++ ) {
        $jack_var+=($jack_mean-$jackarray[$x])**2;
}
# for the sample variance, divide by (n-1)
#print "jack_fst_var ",sprintf("%.9f",$jack_var/($counter-1)),"\n";
my $sterr = sqrt($counter*($jack_var/($counter-1)));
#print "The standard error of the weighted fst is ",sprintf("%.5f",$sterr),"\n";
#print "The 95\%CI of the weighted fst is ",
#sprintf("%.6f",($weighted_Fst_average-1.96*$sterr))," - ",sprintf("%.6f",($weighted_Fst_average+1.96*$sterr)),"\n";


#get col names
my $header = `head -n 1 $inputfile`;
chomp ($header);
my @colnames = split (/,/,$header);

my $pops = @colnames[$chr];

chomp ($pops);
my @popnames = split (/_/,$pops);


# print "\n","\n","\n","summary","\n","\n","\n";
# trying tab seperated fst outputs

#print "pop1","\t","pop2","\t","Fst","\t","CI_1","\t","CI2","\n";

print $popnames[1],"\t",$popnames[2],"\t",sprintf("%.6f",$weighted_Fst_average),"\t",sprintf("%.6f",($weighted_Fst_average-1.96*$sterr)),"\t",sprintf("%.6f",($weighted_Fst_average+1.96*$sterr)),"\n";
}

close DATAINPUT;

```
run it like follwing

```bash
#!/bin/sh
#SBATCH --job-name=bwa_505
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --output=bwa505.%J.out
#SBATCH --error=bwa505.%J.err
#SBATCH --account=def-ben

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

for i in *_only.csv ; do perl test_block_jackknife_withFst.pl ${i} > ${i}_fst_summary.tsv ; done

```

Then download the file and check whether the number of rows to read match the following r code 
Put the output files from above jacknife script in a folder named 'fst_outputs'.
save this r script in in the same location where the 'fst_outputs' folder is.
change the file names to read as needed and run the r script.

All the plots wil be there with the R script


plot_fst_with_CI.R
``` r
#set current path as wd
library("rstudioapi") 
setwd(dirname(getActiveDocumentContext()$path))

library(ggplot2)

# import data tables ( an extra row had been added at the end of the table because of a mistake and removed it here-the 67th was just empty)
l_only_data<-read.table("./fst_outputs/window5k_filtered_exons_l_only_fst_summary.tsv",nrows = 65,header = TRUE)
s_only_data<-read.table("./fst_outputs/window5k_filtered_exons_s_only_fst_summary.tsv",nrows = 65,header = TRUE)
#whole_genome_data<-read.table("./fst_outputs/laevis_GBS_2020_whole_genome_positions_excluded_final_fst.csv_fst_summary.tsv",nrows = 65,header = TRUE)


# make a list for 3 plot cluster
three_plot_list<-list()

# loop through all population combinations

for (i in 1:65) {
  

#prepare l only
l_only_pop_name<-c(1)
l_only_fst<-c(as.numeric(l_only_data[[i,3]]))
l_only_CI_1<-c(as.numeric(l_only_data[[i,4]]))
l_only_CI_2<-c(as.numeric(l_only_data[[i,5]]))

l_only_full_df<-data.frame(l_only_pop_name,l_only_fst,l_only_CI_1,l_only_CI_2)



colnames(l_only_full_df)<-c(as.character(l_only_data[[i,1]]),as.character(l_only_data[[i,2]]),"l_only_CI_1","l_only_CI_2")

#prepare s only

s_only_pop_name<-c(2)
s_only_fst<-c(as.numeric(s_only_data[[i,3]]))
s_only_CI_1<-c(as.numeric(s_only_data[[i,4]]))
s_only_CI_2<-c(as.numeric(s_only_data[[i,5]]))

s_only_full_df<-data.frame(s_only_pop_name,s_only_fst,s_only_CI_1,s_only_CI_2)


colnames(s_only_full_df)<-c(as.character(s_only_data[[i,1]]),as.character(s_only_data[[i,2]]),"s_only_CI_1","s_only_CI_2")

#prepare whole genome( removed this part below in the plot)

#whole_genome_pop_name<-c(2)
#whole_genome_fst<-c(as.numeric(whole_genome_data[[i,3]]))
#whole_genome_CI_1<-c(as.numeric(whole_genome_data[[i,4]]))
#whole_genome_CI_2<-c(as.numeric(whole_genome_data[[i,5]]))

#whole_genome_full_df<-data.frame(whole_genome_pop_name,whole_genome_fst,whole_genome_CI_1,whole_genome_CI_2)


#colnames(whole_genome_full_df)<-c(as.character(whole_genome_data[[i,1]]),as.character(whole_genome_data[[i,2]]),"whole_genome_CI_1","whole_genome_CI_2")

# get column name for geom_point
current_pop_l<-as.character(l_only_data[[i,2]])
current_pop_s<-as.character(s_only_data[[i,2]])

#plot

plot_cluster<-ggplot(data = l_only_full_df,aes(x = l_only_full_df[,1],y = l_only_full_df[,2]))+
  ylim(0,0.5)+
  scale_x_continuous(breaks = c(1,2),labels = c("L","S"))+
  geom_errorbar(data = l_only_full_df, mapping = aes(x = l_only_full_df[,1],y = l_only_full_df[,2], ymin =l_only_CI_1, ymax = l_only_CI_2), size=2, width=.3,color="blue")+
  geom_point(data = l_only_full_df, mapping = aes_string(x = l_only_full_df[,1],y = current_pop_l),shape=15, size=3,color="black")+
  geom_errorbar(data = s_only_full_df, mapping = aes(x = s_only_full_df[,1],y = s_only_full_df[,2], ymin =s_only_CI_1, ymax = s_only_CI_2), size=2, width=.3,color="red")+
  geom_point(data = s_only_full_df, mapping = aes_string(x = s_only_full_df[,1],y = current_pop_s),shape=15, size=3,color="black")+
 #geom_errorbar(data = whole_genome_full_df, mapping = aes(x = whole_genome_full_df[,1],y = whole_genome_full_df[,2], ymin =whole_genome_CI_1, ymax = whole_genome_CI_2), size=1, width=.1,color="black")+
 # geom_point(data = whole_genome_full_df, mapping = aes(x = whole_genome_full_df[,1],y = whole_genome_full_df[,2]),shape=15, size=1,color="grey")+
  theme_bw()+
  ylab(as.character(l_only_data[[i,1]]))+
  xlab(as.character(l_only_data[[i,2]]))+
  theme(text = element_text(size=20,face = "bold"))+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# save plot

three_plot_list[[i]]<-plot_cluster
}

#create final plot adding all plots together

library(gridExtra)
Final_plot_grid<-grid.arrange(
  grobs = three_plot_list,
  ncol=8
  
  #widths = c(2, 2),
  #heights=c(2,1)
)

ggsave(filename = "full_pop_list_plots.pdf",plot = Final_plot_grid,height = 45,width =25 )

# ************ selected plots **************

# select only the populations with at lease 5 samples

selected_plot_list<-list()
i <-1
exclude <- c(5,10,11:20,24,29,32,37,39,44,45,50:55,59,62,64,65)
for(selected_pops in (1:65)[-exclude]){
  selected_plot_list[[i]]<-three_plot_list[[selected_pops]]
  i<-i+1
}

selected_plot_grid<-grid.arrange(
  grobs = selected_plot_list,
  ncol=10
  
  #widths = c(2, 2),
  #heights=c(2,1)
)

ggsave(filename = "selected_pop_list_plots.pdf",plot = selected_plot_grid,height = 45,width =25 )

# create an identical - non visible plot to fill the gaps( just using whatever data I have in here)

plot_blank<-ggplot(data = l_only_full_df,aes(x = l_only_full_df[,1],y = l_only_full_df[,2]))+
  ylim(0,0.5)+
  scale_x_continuous(breaks = c(1,2),labels = c("L only","S only"))+
  geom_errorbar(data = l_only_full_df, mapping = aes(x = l_only_full_df[,1],y = l_only_full_df[,2], ymin =l_only_CI_1, ymax = l_only_CI_2), size=2, width=.3,color="white")+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  #uncomment the following line to see just gridlines
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank())

#preparing plot lists seperately for each of the populations combined with blank plots to create triangle shaped plot

final_plot_list_with_blanks<-list()

# flip all thes axes so sommon axes comes to y ( so I can reduce labelling used)

# **********************

# for Inhaca
# first plot with Y lab
final_plot_list_with_blanks[[1]]<-selected_plot_list[[1]]

#second to end without labs
for (plot_no in 2:8) {
  final_plot_list_with_blanks[[plot_no]]<-selected_plot_list[[plot_no]]+
    theme(axis.title.y = element_blank(),
          axis.text.y=element_blank())
}

# **********************

# for DeDorn

# first plot as a blank

final_plot_list_with_blanks[[9]]<-plot_blank

# plot with Y lab
final_plot_list_with_blanks[[10]]<-selected_plot_list[[9]]

#without labs
for (plot_no in 11:16) {
  final_plot_list_with_blanks[[plot_no]]<-selected_plot_list[[plot_no-1]]+
    theme(axis.title.y = element_blank(),
          axis.text.y=element_blank())
}



# **********************

# for Laigns

# first plots as a blank

for (blank_plot_no in 17:18) {
  final_plot_list_with_blanks[[blank_plot_no]]<-plot_blank
}

# plot with Y lab
final_plot_list_with_blanks[[19]]<-selected_plot_list[[16]]

#without labs
for (plot_no in 20:24) {
  final_plot_list_with_blanks[[plot_no]]<-selected_plot_list[[plot_no-3]]+
    theme(axis.title.y = element_blank(),
          axis.text.y=element_blank())
}



# **********************

# for BW

# first plots as a blank

for (blank_plot_no in 25:27) {
  final_plot_list_with_blanks[[blank_plot_no]]<-plot_blank
}

# plot with Y lab
final_plot_list_with_blanks[[28]]<-selected_plot_list[[22]]

#without labs
for (plot_no in 29:32) {
  final_plot_list_with_blanks[[plot_no]]<-selected_plot_list[[plot_no-6]]+
    theme(axis.title.y = element_blank(),
          axis.text.y=element_blank())
}



# **********************

# for GRNP

# first plots as a blank

for (blank_plot_no in 33:36) {
  final_plot_list_with_blanks[[blank_plot_no]]<-plot_blank
}

#  plot with Y lab
final_plot_list_with_blanks[[37]]<-selected_plot_list[[27]]

# without labs
for (plot_no in 38:40) {
  final_plot_list_with_blanks[[plot_no]]<-selected_plot_list[[plot_no-10]]+
    theme(axis.title.y = element_blank(),
          axis.text.y=element_blank())
}



# **********************

# for Kimber

# first plots as a blank

for (blank_plot_no in 41:45) {
  final_plot_list_with_blanks[[blank_plot_no]]<-plot_blank
}

# plot with Y lab
final_plot_list_with_blanks[[46]]<-selected_plot_list[[31]]

#without labs
for (plot_no in 47:48) {
  final_plot_list_with_blanks[[plot_no]]<-selected_plot_list[[plot_no-15]]+
    theme(axis.title.y = element_blank(),
          axis.text.y=element_blank())
}



# **********************

# for newou

# first plots as a blank

for (blank_plot_no in 49:54) {
  final_plot_list_with_blanks[[blank_plot_no]]<-plot_blank
}

#  plot with Y lab
final_plot_list_with_blanks[[55]]<-selected_plot_list[[34]]

# without labs
for (plot_no in 56) {
  final_plot_list_with_blanks[[plot_no]]<-selected_plot_list[[plot_no-21]]+
    theme(axis.title.y = element_blank(),
          axis.text.y=element_blank())
}



# **********************

# for threesis

# last plots as a blank

for (blank_plot_no in 57:63) {
  final_plot_list_with_blanks[[blank_plot_no]]<-plot_blank
}

# first plot with Y lab
final_plot_list_with_blanks[[64]]<-selected_plot_list[[36]]

#second to end without labs - no plots without labs here


# preparing layout matrix for plot layout
lay<-rbind(c(1,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8),
      c(9,9,9,10,10,10,10,11,11,11,12,12,12,13,13,13,14,14,14,15,15,15,16,16,16),
      c(17,17,17,18,18,18,19,19,19,19,20,20,20,21,21,21,22,22,22,23,23,23,24,24,24),
      c(25,25,25,26,26,26,27,27,27,28,28,28,28,29,29,29,30,30,30,31,31,31,32,32,32),
      c(33,33,33,34,34,34,35,35,35,36,36,36,37,37,37,37,38,38,38,39,39,39,40,40,40),
      c(41,41,41,42,42,42,43,43,43,44,44,44,45,45,45,46,46,46,46,47,47,47,48,48,48),
      c(49,49,49,50,50,50,51,51,51,52,52,52,53,53,53,54,54,54,55,55,55,55,56,56,56),
      c(57,57,57,58,58,58,59,59,59,60,60,60,61,61,61,62,62,62,63,63,63,64,64,64,64))

# print plot list with blanks
selected_plot_grid_with_blanks<-grid.arrange(
  grobs = final_plot_list_with_blanks,
  
  #widths = c(2.5,2,2,2,2,2,2,2),
  layout_matrix =lay 
  
  #heights=c(2,1)
)

ggsave(filename = "selected_pop_list_with_blanks.pdf",plot = selected_plot_grid_with_blanks,height = 75,width =15 ,limitsize = FALSE)


```




## ===========>>>>> POPULATION STRUCTURE ANALYSIS ===========>>>>>>>

This starts with the finalized bam files from above

First copy only the bam files to a new directory
in the directory with the finalized_bam_files folder from above,

```bash
#!/bin/sh
#SBATCH --job-name=bwa_505
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:00:00
#SBATCH --mem=100gb
#SBATCH --output=bwa505.%J.out
#SBATCH --error=bwa505.%J.err
#SBATCH --account=def-ben

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

find ./finalized_bam_files/ -name '*.bam' | cpio -pdm  pop_structure/

```

Then *in the pop_structure folder* run the following.
This will,
load ANGSD, Set paths to the programs and the data(use the corresponding directories) and get angsd outputs for all genomes and collect them in corresponding folders
```
#!/bin/sh
#SBATCH --job-name=bwa_505
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --mem=32gb
#SBATCH --output=bwa505.%J.out
#SBATCH --error=bwa505.%J.err
#SBATCH --account=def-ben

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

module load nixpkgs/16.09
module load intel/2018.3
module load angsd/0.929

AA=/scratch/premacht/xlaevis_and_xgilli/ANGSD

for j in finalized_bam_files/*; do find ${j} | grep bam$ > ${j}_file_list
        angsd -bam ${j}_file_list -GL 2 -doMajorMinor 1 -doMaf 1 -SNP_pval 2e-6 -minMapQ 20 -minQ 20 -doCounts 1 -doDepth 1 -setMinDepth 2 -setMaxDepth 100  -minInd 1 -minMaf 0.05 -doGlf 2 -out ${j}_angsd_output -p 1; done

        mkdir angsd_outputs
        for i in l_only s_only whole_genome; do mkdir angsd_outputs/${i}; mv finalized_bam_files/*${i}_angsd_output* angsd_outputs/${i} ; done
```
# Download and install NGSadmix 

First, save this python script on saved_scripts_to_be_used (whereveer you keep all the re-using scripts and change address accordingly in the future) as clumpp_input_maker_new.py

```python
"""

This will take files from multiple runs of a particular K (from, say, Structure, or NGSAdmix), create a clumpp input file, and recommend the clumpp algorithm to choose (by calculating the D stat from the Clumpp manual).


Ex.

python3.7 clumpp_input_maker.py -in test/run*/2*qopt -type ngsadmix -out test/forClumpp

## produces test/forClumpp_clumpp.param & test/forClumpp_clumpp.in

"""


def get_cmdline_args():
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument("-in", dest='input_files', nargs = '+',
						help = "pattern for input file name. Something unique"
						"to grab all names. should be able to search through"
						"all subdirectories holding the replicates (run_1/,"
						"run_2/, etc.): E.g., pattern: run*/2_*qopt ",
						required = True
						)
	parser.add_argument("-type", "--cluster_type",
						help = "was this and 'ngsadmix' or 'structure' run?",
						required = True
						)
	# parser.add_argument("-param", "--write_clumpp_param", help = "write a parameter file of this name for clumpp based on recommendations here.")
	parser.add_argument("-out", "--base_output_name",
						help = "base name to add to the output files for"
						 "clumpp"
						 )
	args = parser.parse_args()
	return(args)

def read_input(name):
	"""
	take a file name, read it in.
	return - list of lines
	"""
	input = list()
	IN = open(name)
	for line in IN.readlines():
		input.append(line)
	return(input)


def collect_input_data(files, filetype):
	"""
	read in the input files and collect the data
	"""
	data = list()
	for n in files:
		data.append(read_input(n))

	return(data, len(files))

def get_num_individuals(data):
	"""
	determines the number of individuals in the analysis
	"""
	return(len(data[0]))

def get_k(data):
	"""
	determines the number of clusters
	"""
	indv1 = data[0][0]
	indv1 = indv1.strip().split()# "\t"
	return(len(indv1)) ## THIS ASSUMES IT IS ONLY CLUSTER ASSIGNMENT DATA HERE

def compute_D(C, K, R, N = 1000):
	"""
	C: # individuals,
	K: # clusters,
	R: # replicates,
	N: number of replicates to test
	"""
	from math import factorial

	try: # b/c some values, the factorial will be too large for python to handle
		T_fullsearch = (factorial(K) ** (R-1)) * (R * (R-1)/2) * K * C # make sure this math works
		D_fullsearch = T_fullsearch * C * 1

		T_greedy = (factorial(K)) * (R * (R-1)/2) * K * C
		D_greedy = T_greedy * C * N
	except:
		return(3)

	# T_largeK = (R * (R − 1)/2) * (K ** 2) * C

	if D_fullsearch <= 10 ** 13:
		print("\n use fullsearch, M = 1")
		return(1)
	elif D_greedy <= 10 ** 13:
		print("\n use Greedy, M = 2, ideally a greedy_option = 1 would be set, but that may be too slow (all possible orders tested), greedy_option = 2 could be used with high R (e.g., 1000). But if K gets large, may need greedy_option = 3 and high R (e.g., 1000).")
		return(2)
	else:
		print("\n use largeKGreedy, M = 3, again, ideally greedy_option 1, but most likely it will be 2 or 3. Use large R for 2 or 3.")
		return(3)

def print_clumpp_input(data, basename = 'forClumpp', type = 'ngsadmix'):
	"""
	make input data file for clumpp
	"""
	import re
	outfile = basename + "_clumpp.in"

	## need to add the starting 5 columns to this including an individual identifier, which will not exist for NGSADMIX.
	with open(outfile, 'w') as out_fh:
		for set in data:
			# for indv in set:
			for i in range(1, len(set)+1):
				if re.search("ngsadmix", type, re.I):
					# create set of indvidual IDs and other expected columns
					line = [str(i), i, "(0)", 99, ":"]
					for l in line:
						out_fh.write("{}\t".format(l))
					out_fh.write("{}".format(set[i-1]))
					# print(set[i-1], end = '')
				else:
					print(indv, end = '')
					out_fh.write("{}".format(set[i-1]))
			out_fh.write("\n")


def print_param_file(C, K, R, M, basename):
	"""
	prints the parameter file for clumpp
	"""
	print("\nSetting GREEDY_OPTION 2, and testing 1000 repeats of input order. Unless M = 1, in which case these will be ignored.\n")

	param_outf = basename + "_clumpp.param"
	header = ["DATATYPE\t0","INDFILE\t" + basename + "_clumpp.in",
				"POPFILE\tnotneeded.popfile","OUTFILE\t" + basename + "_clumpp.out",
				"MISCFILE\tnotneeded.miscfile", "K\t" + str(K), "C\t" + str(C) , "R\t" + str(R),
				"M\t" + str(M), "W\t0", "S\t2", "GREEDY_OPTION\t2", "REPEATS\t1000","PERMUTATIONFILE\tNOTNEEDED.permutationfile",
				"PRINT_PERMUTED_DATA\t0",
				"PERMUTED_DATAFILE\tk4.perm_datafile",
				"PRINT_EVERY_PERM\t0\t",
				"EVERY_PERMFILE\tnotneeded.every_permfile",
				"PRINT_RANDOM_INPUTORDER\t0",
				"RANDOM_INPUTORDERFILE\tk4.random_inputorderfile",
				"OVERRIDE_WARNINGS\t0",
				"ORDER_BY_RUN\t0"]
	with open(param_outf, 'w') as out_fh:
		for l in header:
			out_fh.write("{}\n".format(l))

# DATATYPE 0
# INDFILE k2_forClumpp
# POPFILE notneede.popfile
# OUTFILE k2.outfile
# MISCFILE NOTNEEDED.miscfile
# K 2
# C 69
# R 15
# M 2
# W 0
# S 2
# GREEDY_OPTION 2
# REPEATS 1000
# "PERMUTATIONFILE\tNOTNEEDED.permutationfile",
# "PRINT_PERMUTED_DATA\t0",
# "PERMUTED_DATAFILE\tk4.perm_datafile",
# "PRINT_EVERY_PERM\t0\t",
# "EVERY_PERMFILE\tnotneeded.every_permfile",
# "PRINT_RANDOM_INPUTORDER\t0",
# "RANDOM_INPUTORDERFILE\tk4.random_inputorderfile",
# "OVERRIDE_WARNINGS\t0",
# "ORDER_BY_RUN\t0"





def main():
	import re


	args = get_cmdline_args()

	if not re.search("ngsadmix", args.cluster_type, re.I):
		exit("\n\nonly for ngsadmix data right now. Assumes the input cluster files only have cluster assignment on each line, no other information (Q-files).\n\n")

	cluster_data, R = collect_input_data(args.input_files, args.cluster_type)

	C = get_num_individuals(cluster_data)

	K = get_k(cluster_data)

	if re.search("ngsadmix", args.cluster_type, re.I):
		print_clumpp_input(
							cluster_data,
							args.base_output_name,
							args.cluster_type
							)

	M = compute_D(C, K, R) ## not actually doing anything about this. add that.

	print_param_file(C, K, R, M, args.base_output_name)

if __name__ == '__main__':
	main()
 ```
Then,

In pop_structure folder
Install CLUMPP and NGSadmix
```bash
mkdir CLUMPP
cd CLUMPP
wget https://rosenberglab.stanford.edu/software/CLUMPP_Linux64.1.1.2.tar.gz

gunzip CLUMPP_Linux64.1.1.2.tar.gz

tar xvf CLUMPP_Linux64.1.1.2.tar

cd ../

mkdir ngsadmix
cd ngsadmix
wget https://raw.githubusercontent.com/ANGSD/angsd/master/misc/ngsadmix32.cpp
g++ ngsadmix32.cpp -O3 -lpthread -lz -o NGSadmix
cd ../
ls
```
Run this in /angsd_outputs/l_only to get multiple runs of ngsadmix outputs , create clumpp inputs, 

```bash


for i in ../*; do 
cd ${i}
for run in `seq 10`;   do   mkdir run_$run ; cd run_$run/;   for K in `seq 5`;     do ../../../ngsadmix/NGSadmix -likes ../*.beagle.gz -K $K -P 10 -o $K\_outfiles -minMaf 0.05;   done;   cd ../ ; done
cp /scratch/premacht/laevis_GBS_2020/saved_scripts_to_be_used/clumpp_input_maker_new.py .
for k in `seq 5` ; do python3 clumpp_input_maker_new.py -in run_*/${k}_*qopt -type ngsadmix -out k$k ; done
for f in k*param ; do ../../CLUMPP/CLUMPP_Linux64.1.1.2/CLUMPP $f; done & done

```

Run this in the pop_structure folder to collect all the needed files to be downloaded

```bash
find angsd_outputs/*/run_* -name '*' | cpio -pdm all_final_outputs_for_R/
find angsd_outputs/*/k* -name '*' | cpio -pdm all_final_outputs_for_R/
mkdir all_final_outputs_for_R/angsd_outputs/file_lists
cp finalized_bam_files/*file_list all_final_outputs_for_R/angsd_outputs/file_lists/
```

Run this in all_final_outputs_for_R/l_only to arrange clumpp outputs, runs and file lists
```bash
for j in ../l_only  ../s_only  ../whole_genome ; 
do cd ${j}  
mkdir runs
mv -f run* runs/
mkdir clumpp_files
mv k* clumpp_files
mv ../file_lists/*${j##../}* . ; done
rmdir ../file_lists/
```
Then you can download all the needed files for R with a single command from desktop

```bash
scp -r premacht@graham.computecanada.ca:/scratch/premacht/laevis_GBS_2020/testing_ground/pop_structure/all_final_outputs_for_R/angsd_outputs ./laevis_GBS_2020/
```

# ========>>>>> CONTINUES IN THE DESKTOP MACHINE ===========>>>>>>>

after the files are downloaded,

copy mobile pack folder into the directory where we have
```
1. clumpp_files folder
2. runs folder
3. downloaded file list
```
mobile pack folder has 3 Rscripts

# 1. get_a_file_for_location_info.r

Run this and it will create a directory called locality info and a file named full_summary inside it

```r
# This outputs a file to enter locations for samples

# set wd to the place where file is
library("rstudioapi") 
setwd(dirname(getActiveDocumentContext()$path)) 

# get the name of the file list file name
list_of_files<-list.files(path = './..',pattern = 'file_list' )
# read list_of_files
sample_address_list<-readLines(paste("./../",list_of_files,sep =""))

# extract file names
sample_list<-basename(sample_address_list)


# make an empty dataframe for sample info and add file list

sample_info<-setNames(data.frame(matrix(ncol = 5, nrow = length(sample_list))), c("IND","ID","LOCATION","SPECIES","FILE"))
sample_info[,5]<-sample_list

# Fill IND automatically
sample_info[,1]<-paste("IND",seq.int(nrow(sample_info)),sep = "")

# Give default IDs(just bumbers to be edited later)
sample_info[,2]<-paste(seq.int(nrow(sample_info)),"enter_real_ID_here",sep = "_")

# Default value for locations
sample_info[,3]<-paste(seq.int(nrow(sample_info)),"enter_location_here",sep = "_")

# Default value for species
sample_info[,4]<-"enter_species_here"

# make a directory for location info
dir.create(path = './locality_info')
write.table(sample_info,"./locality_info/full_summary.txt",sep="\t",row.names=FALSE,quote = FALSE)
```

*Open the created file and edit the essential data in it.(by default this script assumes all the samples are from the same species but different populations)
*Then save it without changing file name or format

then run the second file inside mobile pack

# 2. make lnlk

```r

# set wd to the place where file is
library("rstudioapi") 
setwd(dirname(getActiveDocumentContext()$path)) 

# get the name of the file list file name
list_of_files<-list.files(path = './..',pattern = 'file_list' )
# read list_of_files
sample_address_list<-readLines(paste("./../",list_of_files,sep =""))

# extract file names
sample_list<-basename(sample_address_list)


# make an empty dataframe for sample info and add file list

sample_info<-setNames(data.frame(matrix(ncol = 5, nrow = length(sample_list))), c("IND","ID","LOCATION","SPECIES","FILE"))
sample_info[,5]<-sample_list

# Fill IND automatically
sample_info[,1]<-paste("IND",seq.int(nrow(sample_info)),sep = "")

# Give default IDs(just bumbers to be edited later)
sample_info[,2]<-paste(seq.int(nrow(sample_info)),"enter_real_ID_here",sep = "_")

# Default value for locations
sample_info[,3]<-paste(seq.int(nrow(sample_info)),"enter_location_here",sep = "_")

# Default value for species
sample_info[,4]<-"enter_species_here"

# make a directory for location info
dir.create(path = './locality_info')
write.table(sample_info,"./locality_info/full_summary.txt",sep="\t",row.names=FALSE,quote = FALSE)
```
This creates the files bam_names_list and lnlks_allRuns

Then we can run the third script to plot all data

# 3. plot_nsgadmix

```r
######
# plot Ks and lnlks
######

#set current path as wd
library("rstudioapi") 
setwd(dirname(getActiveDocumentContext()$path))
#import all locations to sort
#****************************#
# you can observe sort list from here and might have to paste output for levels below

library(tidyverse)
require("readxl")
library("readxl")
library(dplyr)
library(cowplot)
library(ggplot2)

full_table<-read_tsv("./locality_info/full_summary.txt")
location_list<-(distinct(full_table[3]))

#************************#
#used desc to take vic to top
sort_list<-arrange(.data = location_list,desc(LOCATION)) %>%toString %>% noquote
cat(sort_list)
#************************#

## ---- Functions -----
get_k <- function(file, k_pos = 1) {
  basename(file) %>% str_split('[_*]|K') %>%
    unlist() %>% pluck(k_pos)
}

import_names_file <- function(file){
  ids <- read_tsv(file, col_names = T)
  names(ids) <- c("indv","location","ID","species")
  ids
}

import_clumpp <- function(file, ids) {
  # need file of names in the name order as the bams, which should be same
  # as the clumpp order.
  clust <- read_table(file, col_names = F)
  clust <- read_table(file, col_names = F) %>%
    mutate(k = paste('K', ncol(clust[,c(6:ncol(clust))]), sep = ''))
  
  clust <- bind_cols(clust, ids)
  
  
  ##======>>>YOU CAN COLLECT LOCATION LIST HERE AND ARRANGE IN THE FOLLOWING LINE =====>>>
  
  automated_list<-c(location_list)
 
  
  clust %>%
    
    gather(key = 'clust', value = 'prop', -c(X1, X2,X3,X4,X5, k, indv, location, ID,species)) %>%
    mutate(location = factor(location,
                             
                             # use default order for location order like this
                             
                             #levels = unlist(automated_list, use.names=FALSE),
                             
                             #==>Or extract the location list by running 'levels' above and change oreder here
                             levels = c('GRNP', 'DeDorn', 'Vred','Citrus', 'BW', 'Laigns', 'Niewou','Kimber', 'Draken', 'no_lbl','Inhaca','VicW', 'ThreeS'),
                             
                             #====and==>>
                             
                             #use numbers instead location names
                             #labels = 1:nrow(location_list)
                             
                             # or use real names
                             #labels = unlist(automated_list, use.names=FALSE)
                             
                             # ==> and labels for pre assigned levels if you changed it above
                             
                             labels = c('GRNP', 'DeDorn', 'Vred','Citrus', 'BW', 'Laigns', 'Niewou','Kimber', 'Draken', 'no_lbl','Inhaca','VicW', 'ThreeS')
    )
    )
  
}

##========================================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#without labs
plot_k_wtout_labs <- function(clumpp) {
  ggplot(clumpp, aes(x = indv, y = prop, fill = clust)) +
    geom_col(width = 1) +
    facet_grid( ~ location+species, switch = "x"
                , scales = "free_x"
                , space = 'free') +
    theme_bw() +
    theme(
      legend.position = 'none',
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      strip.background.x = element_blank(),
      
      #uncomment following if you want file names
      #axis.text.x = element_text(angle = 90),
      
      
      axis.text.x  = element_blank(),
      
      #Change label size and angle here/ set this to 0 to remove labels from first 3 plots to reduce conflicts
      strip.text.x = element_text(size = 0),
      
      axis.ticks.x = element_blank(),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 16),
      # remove individual x labelling for all plots
      # strip.text.x = element_blank()
    ) +
    # scale_fill_brewer(type = "qual",palette = "Paired") +
    scale_fill_manual(
      values =get(paste("pal",i,sep = "",collapse = NULL)),
      limits = names(pal)
    )+
    xlab("") +
    scale_y_continuous(breaks = c(0, 0.5, 1.0))
  
}

## --- end of functions




## -- Handle lnlk
lnlks <- read_tsv('lnlks_allRuns.txt', col_names = T)
names(lnlks) <- c("rep","k",'lnlk')

quants <-lnlks %>%
  group_by(k) %>%
  summarise(med = median(lnlk),
            upper = quantile(lnlk, 0.975),
            lower = quantile(lnlk, 0.025)
  )



#final plot with labs
plot_k_wt_labs <- function(clumpp) {
  ggplot(clumpp, aes(x = indv, y = prop, fill = clust)) +
    geom_col(width = 1) +
    facet_grid( ~ location+species, switch = "x"
                , scales = "free_x"
                , space = 'free') +
    theme_bw() +
    theme(
      legend.position = 'none',
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      strip.background.x = element_blank(),
      
      #uncomment following if you want file names
      #axis.text.x = element_text(angle = 90),
      
      
      axis.text.x  = element_blank(),
      
      # **** Change label size and angle here ****
      strip.text.x = element_text(size = 20, face = "italic",angle = 60),
      
      #axis.ticks.x = element_blank(),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 16),
      #remove individual x labelling for all plots
      #strip.text.x = element_blank()
    ) +
    # scale_fill_brewer(type = "qual",palette = "Paired") +
    scale_fill_manual(
      values = pal,
      limits = names(pal)
    )+
    xlab("") +
    scale_y_continuous(breaks = c(0, 0.5, 1.0))
  
}

## --- end of functions


## -- Handle lnlk
lnlks <- read_tsv('lnlks_allRuns.txt', col_names = T)
names(lnlks) <- c("rep","k",'lnlk')

quants <-lnlks %>%
  group_by(k) %>%
  summarise(med = median(lnlk),
            upper = quantile(lnlk, 0.975),
            lower = quantile(lnlk, 0.025)
  )


lnlk_multi_plot <-
  ggplot(filter(quants, k < 6), aes(x = k, y = med, group = k)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, size = 1) +
  geom_point(size = 5)  +
  theme_classic() +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16)
  ) +
  xlab(expression(paste("Number of Clusters (", italic("K"),")"))) +
  ylab(expression(italic(log-likelihood)))


## -- handle cluster plots

sample_names <- import_names_file("bam_names.txt")
## In the same order as samples were given to NgsAdmix (list of bams order)
## E.g., (tab seperated)
# AMNH17273_Xt_SierL_male	sierra_leone	AMNH17273
# XEN091_Xt_SierL_female	sierra_leone	XEN091
# XEN092_Xt_SierL_female	sierra_leone	XEN092
# XEN094_Xt_SierL_male	sierra_leone	XEN094


c_files <- list.files("./../clumpp_files/", "out", full.names = T) %>%
  set_names(., map(., get_k))
# A directory that contains files named: k4_clumpp.out, k5_clumpp.out, etc

clumpp_dats <- map(c_files, import_clumpp, ids = sample_names)

#set colors manually
#assign colour set
left_corner<-"darkblue"
right_corner<-"lightblue"
left_middle<-"forestgreen"
right_middle<-"pink"
upper_little<-"purple"

pal1 <- c(
  "X6" = "lightblue",
  "X7" = "darkblue", 
  "X8" = "pink", 
  "X9" = "forestgreen",
  "X10"= "purple"
)

pal2 <- c(
  "X6" = "pink",
  "X7" = "lightblue", 
  "X8" = "darkblue", 
  "X9" = "forestgreen",
  "X10"= "purple"
)

pal3 <- c(
  "X6" = "lightblue",
  "X7" = "darkblue", 
  "X8" = "pink", 
  "X9" = "forestgreen",
  "X10"= "purple"
)

pal4 <- c(
  "X6" = "purple",
  "X7" = "darkblue", 
  "X8" = "pink", 
  "X9" = "forestgreen",
  "X10"= "lightblue"
)

pal <- c(
  "X6" = "pink",
  "X7" = "lightblue", 
  "X8" = "purple", 
  "X9" = "forestgreen",
  "X10"= "darkblue"
)


# ****change cluster order here*********#

# first cluster change

#change k here
k<-"k4"

#sort data in the end
attach(clumpp_dats[[k]])

#change the priority order of clusters here ** the cluster at on the right end gets plotted at the bottom compared to one in the left
#( list the X value related to the color that should be in the top in the left most end of the following list)
clumpp_dats[[k]][["clust"]] <- factor(clumpp_dats[[k]][["clust"]], levels = c("X10","X9", "X8","X6","X7"))

clumpp_dats[[k]] <- clumpp_dats[[k]][order(clust),]

# second cluster change

#change k here
k<-"k5"

#sort data in the end
attach(clumpp_dats[[k]])

#change the priority order of clusters here ** the cluster at on the right end gets plotted at the bottom compared to one in the left
#( list the X value related to the color that should be in the top in the left most end of the following list)
clumpp_dats[[k]][["clust"]] <- factor(clumpp_dats[[k]][["clust"]], levels = c("X8","X9", "X6","X7","X10"))

clumpp_dats[[k]] <- clumpp_dats[[k]][order(clust),]
# ******************************************
plot_name_list<-1:4
for (i in 1:2) {
  plot_name_list[i] <- map(clumpp_dats[i],plot_k_wtout_labs)
  i<-i+1
}




#plots with manual orientation
i<-3
plot_name_list[i] <- map(clumpp_dats[i],plot_k_wtout_labs)

#k_plots <- map(clumpp_dats[1:3], plot_k_wtout_labs)

#last plot with labels
plot_name_list[4] <- map(clumpp_dats[4], plot_k_wt_labs)

library(gridExtra)
multi_k_plots<-grid.arrange(
  grobs = plot_name_list,
  ncol=1
  
  #widths = c(2, 2),
  #heights=c(2,1)
)



# Create a folder for outputs
dir.create(path = "./../final_plots")

ggsave("./../final_plots/K2-5_ngsadmix.pdf", multi_k_plots,width = 30, height = 20)


## -- combine lnlk and cluster plots
join_plots <-
  cowplot::plot_grid(multi_k_plots, lnlk_multi_plot,
                     labels = "AUTO", ncol = 2
  )
xx<-ggdraw(add_sub(join_plots, "Label", vpadding=grid::unit(0,"lines"),y=6, x=0.5, vjust=4.5))


ggsave("./../final_plots/k_lnlk_and_admix_plots.pdf", join_plots,
       useDingbats=FALSE, width = 30, height = 20)





```
# ===>>> This will create all the finalized plots and they can be found in a newly created directory called final_plots in the parent directory

# Do same with all the directories you have for different genomes


# ************************************************************************************************************************************************************* END OF POPULATION STRUCTURE ANALYSIS ************************************************************************************
























# =========>>>>>  FOR ABBABABA ============>>>>

**To be used for ABBABABA** it is necessary to swap any astrisks with Ns: run following in geno_files folder
then move older files to a new directory
```bash
for j in *.geno; do sed -i 's/\*/N/g' ${j}
gzip -c ${j} > ${j%%.geno}_astrisks_swapped.geno ; done
mkdir old_genos
mv *excluded.geno old_genos/
```







 ## ************************* DID NOT USE BELOW THIS POINT IN FINAL ATTEMPT  *******************************

# Cal depth and move depth files to the folder(If needed. Not going to use here now)
********** (always check file size after calculation to make sure the used region was present in all samples. If not, change the region)********
bash script
```bash
#!/bin/sh
module load bwa
module load samtools/1.10
for i in ./../../bam_files/*.bam ; do samtools depth -r chr1L:1-10000000 $i > $i"_depth" ; done
mv ./../../bam_files/*_depth ./../depth_files
```
run by
```bash
bash cal_depth.sh
```


# plot depth
```rscript
library (ggplot2)
#set working directory here
working_directory<-"./../depth_files"

setwd(working_directory)
#remove scientific notation
options(scipen=999)

#*******create concat file****

#set path with raw files here
raw_files_are_in<-working_directory

#create output directory
dir.create(paste(working_directory,"/out",sep = ""))

#save concat file in
save_path<-paste(working_directory,"/out",sep = "")

require(dplyr)
require(data.table)


# read file path raw
all_paths_raw_pre <-
          list.files(path = raw_files_are_in,
                                  full.names = TRUE,recursive = FALSE)
all_paths_raw<-all_paths_raw_pre[!file.info(all_paths_raw_pre)$isdir]

#create first data frame
x=1

all_data<-read.table(file=all_paths_raw[1])
colnames(all_data)<-c("chr","pos","dep")
all_data["sample_name"]=basename(all_paths_raw[1])
all_data$sample_no<-x

#add data
for (i in all_paths_raw[2:length(all_paths_raw)]) {
          x=x+1
  
  current_table<-read.table(i)
    setnames(current_table,c("chr","pos","dep"))
    current_table["sample_name"]=basename(i)
      current_table$sample_no<-x
      all_data<-rbind(all_data,current_table)
        
}
write.table(all_data,file =paste(save_path,"/concat",sep=''))

#*****end of concat******

#set new working directory
setwd(save_path)

all_data<-read.table(file = "concat")

p <- ggplot(data = all_data, aes(x=sample_no, y=dep))  +
          ylim(0,60)+
            ylab("Depth")+
              xlab("sample_number")+
                geom_boxplot(aes(fill=sample_name))+
                 theme(legend.position = "none")

q <- ggplot(data = all_data, aes(x=sample_no, y=dep))  +
                  ylim(0,60)+
                              ylab("Depth")+
                                            xlab("sample_number")+
                                                            geom_boxplot(aes(fill=sample_name))

        ggsave(filename = "all_samples_without_legend.pdf",plot = p,height = 10,width = 20)
        ggsave(filename = "all_samples_with_legend.pdf",plot = q,height = 10,width = 20)
```
# run plot_depth
```bash
#!/bin/sh
#SBATCH --job-name=bwa_505
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=64gb
#SBATCH --output=bwa505.%J.out
#SBATCH --error=bwa505.%J.err
#SBATCH --account=def-ben

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

module load nixpkgs/16.09
module load gcc/7.3.0
module load r

Rscript plot_depth.r
```
Removed a file with a low depth before the next steps.
# ********************************************************

# copied bam files into
```bash
/scratch/premacht/laevis_GBS_2020/bams_for_seperate_genomes/bams_whole_genome
```
# Seperate subgenomes

After creating seperate directories for subgenomes

```bash
module load nixpkgs/16.09  intel/2016.4
module load samtools/0.1.17

for i in *.bam; do samtools index ${i}; done

for i in *.bam;do samtools view -b $i chr1L chr2L chr3L chr4L chr5L chr6L chr7L chr8L chr9_10L > ../bams_L_only/L_only_${i} ;done


for i in *.bam;do samtools view -b $i chr1S chr2S chr3S chr4S chr5S chr6S chr7S chr8S chr9_10S > ../bams_S_only/S_only_${i} ;done
```
# Cal GL

In each of the folders for subgenomes..
```bash
#!/bin/sh
#SBATCH --job-name=bwa_505
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=64gb
#SBATCH --output=bwa505.%J.out
#SBATCH --error=bwa505.%J.err
#SBATCH --account=def-ben

#SBATCH --mail-user=tharindutaemail@gmail.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL



module load nixpkgs/16.09
module load intel/2018.3
module load angsd/0.929

AA=/scratch/premacht/xlaevis_and_xgilli/ANGSD
BAMFOLDER=../../bams_for_seperate_genomes/bams_L_only
find $BAMFOLDER |  grep bam$ > all.files

angsd -bam all.files -GL 2 -doMajorMinor 1 -doMaf 1 -SNP_pval 2e-6 -minMapQ 20 -minQ 20 -doCounts 1 -doDepth 1 -setMinDepth 2 -setMaxDepth 100  -minInd 15 -minMaf 0.05 -doGlf 2 -out all -P 1
```

# Then in each seperate folders

```bash
~/../../scratch/premacht/xlaevis_and_xgilli/NDSadmix/NGSadmix -likes all.beagle.gz -K 3 -minMaf 0.05 -seed 1 -o all

for run in `seq 10`;   do   mkdir run_$run ; cd run_$run/;   for K in `seq 5`;     do ~/../../scratch/premacht/xlaevis_and_xgilli/NDSadmix/NGSadmix -likes ../all.beagle.gz -K $K -P 10 -o $K\_outfiles -minMaf 0.05;   done;   cd ../; done

module load python/3.5.4

cp ~/../../scratch/premacht/xlaevis_and_xgilli/clumpp_input_maker_new.py .

for k in `seq 5` ; do python3 clumpp_input_maker_new.py -in run_*/${k}_*qopt -type ngsadmix -out k$k ; done

for f in k*param ; do ../CLUMPP_Linux64.1.1.2/CLUMPP $f ; done

```

















# ****************************** Cal Fst ***************************

# Background
Selected the bam files with a good depth value from above depth check step and used those bam files from here onwards

# Create a VCF using selected bam files

Downloaded reference genome 
```bash
wget http://ftp.xenbase.org/pub/Genomics/JGI/Xenla9.2/XENLA_9.2_genome.fa.gz
```
Removed scaffolds from reference genome

```bash
module load bwa
module load samtools/1.10

samtools faidx XENLA_9.2_genome.fa chr1S chr2L chr2S chr3L chr3S chr4L chr4S chr5L chr5S chr6L chr6S chr7L chr7S chr8L chr8S chr9_10L chr9_10S > XENLA_9.2_genome_chrs_only.fs
# In bam folder ,
Observed addresses for all bam files
```bash
ls *.bam | tr "\n" " "
```

Then used following script to create VCF using observed file list,Remove non varient sites, remove scaffolds , remove scaffold labels from header and seperate into subgenomes


```bash
#!/bin/sh
#SBATCH --job-name=bwa_505
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=512gb
#SBATCH --output=bwa505.%J.out
#SBATCH --error=bwa505.%J.err
#SBATCH --account=def-ben

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

module load bwa
module load samtools/1.10
module load nixpkgs/16.09
module load intel/2018.3
module load bcftools/1.10.2
module load vcftools/0.1.16

samtools mpileup -q20 -d8000 -ugf ../fst/reference_genome/XENLA_9.2_genome.fa 2014_Inhaca_10_Inhaca_ATATGT_cuttrim_sorted.bam 2014_Inhaca_150_Inhaca_ATCGTA_cuttrim_sorted.bam 2014_Inhaca_152_Inhaca_CATCGT_cuttrim_sorted.bam 2014_Inhaca_24_Inhaca_CGCGGT_cuttrim_sorted.bam 2014_Inhaca_38_Inhaca_CTATTA_cuttrim_sorted.bam 2014_Inhaca_52_Inhaca_GCCAGT_cuttrim_sorted.bam 2014_Inhaca_65_Inhaca_GGAAGA_cuttrim_sorted.bam 946_Draken_TCGTT_cuttrim_sorted.bam 993_Draken_GGTTGT_cuttrim_sorted.bam BJE3508_DeDorn_ATTGA_cuttrim_sorted.bam BJE3509_DeDorn_CATCT_cuttrim_sorted.bam BJE3510_DeDorn_CCTAG_cuttrim_sorted.bam BJE3511_DeDorn_GAGGA_cuttrim_sorted.bam BJE3512_DeDorn_GGAAG_cuttrim_sorted.bam BJE3513_DeDorn_GTCAA_cuttrim_sorted.bam BJE3514_DeDorn_TAATA_cuttrim_sorted.bam BJE3515_DeDorn_TACAT_cuttrim_sorted.bam BJE3525_Laigns_GAATTCA_cuttrim_sorted.bam BJE3526_Laigns_GAACTTG_cuttrim_sorted.bam BJE3527_Laigns_GGACCTA_cuttrim_sorted.bam BJE3528_Laigns_GTCGATT_cuttrim_sorted.bam BJE3529_Laigns_AACGCCT_cuttrim_sorted.bam BJE3530_Laigns_AATATGG_cuttrim_sorted.bam BJE3531_Laigns_ACGTGTT_cuttrim_sorted.bam BJE3532_Laigns_ATTAATT_cuttrim_sorted.bam BJE3533_Laigns_ATTGGAT_cuttrim_sorted.bam BJE3534_BW_CTCG_cuttrim_sorted.bam BJE3535_BW_TGCA_cuttrim_sorted.bam BJE3536_BW_ACTA_cuttrim_sorted.bam BJE3537_BW_CAGA_cuttrim_sorted.bam BJE3538_BW_AACT_cuttrim_sorted.bam BJE3539_BW_GCGT_cuttrim_sorted.bam BJE3541_BW_CGAT_cuttrim_sorted.bam BJE3542_BW_GTAA_cuttrim_sorted.bam BJE3543_BW_AGCG_cuttrim_sorted.bam BJE3544_BW_GATG_cuttrim_sorted.bam BJE3545_BW_TCAG_cuttrim_sorted.bam BJE3546_BW_TGCGA_cuttrim_sorted.bam BJE3547_GRNP_TAGGAA_cuttrim_sorted.bam BJE3548_GRNP_GCTCTA_cuttrim_sorted.bam BJE3549_GRNP_CCACAA_cuttrim_sorted.bam BJE3550_GRNP_CTTCCA_cuttrim_sorted.bam BJE3551_GRNP_GAGATA_cuttrim_sorted.bam BJE3552_GRNP_ATGCCT_cuttrim_sorted.bam BJE3553_GRNP_AGTGGA_cuttrim_sorted.bam BJE3554_GRNP_ACCTAA_cuttrim_sorted.bam BJE3573_VicW_CGCGGAGA_cuttrim_sorted.bam BJE3574_VicW_CGTGTGGT_cuttrim_sorted.bam BJE3575_Kimber_GTACTT_cuttrim_sorted.bam BJE3576_Kimber_GTTGAA_cuttrim_sorted.bam BJE3578_Kimber_TGGCTA_cuttrim_sorted.bam BJE3579_Kimber_TATTTTT_cuttrim_sorted.bam BJE3580_Kimber_CTTGCTT_cuttrim_sorted.bam BJE3581_Kimber_ATGAAAG_cuttrim_sorted.bam BJE3582_Kimber_AAAAGTT_cuttrim_sorted.bam BJE3632_Niewou_CATAAGT_cuttrim_sorted.bam BJE3633_Niewou_CGCTGAT_cuttrim_sorted.bam BJE3640_Niewou_CGGTAGA_cuttrim_sorted.bam BJE3641_Niewou_CTACGGA_cuttrim_sorted.bam BJE3642_Niewou_GCGGAAT_cuttrim_sorted.bam BJE3644_Niewou_TAGCGGA_cuttrim_sorted.bam BJE3645_Niewou_TCGAAGA_cuttrim_sorted.bam BJE3647_Niewou_TCTGTGA_cuttrim_sorted.bam BJE3654_ThreeSis_TGCTGGA_cuttrim_sorted.bam BJE3655_ThreeSis_ACGACTAG_cuttrim_sorted.bam BJE3656_ThreeSis_TAGCATGG_cuttrim_sorted.bam BJE3657_ThreeSis_TAGGCCAT_cuttrim_sorted.bam BJE3658_ThreeSis_TGCAAGGA_cuttrim_sorted.bam BJE3659_ThreeSis_TGGTACGT_cuttrim_sorted.bam BJE3660_ThreeSis_TCTCAGTG_cuttrim_sorted.bam BJE3661_ThreeSis_CGCGATAT_cuttrim_sorted.bam BJE3662_ThreeSis_CGCCTTAT_cuttrim_sorted.bam BJE3663_ThreeSis_AACCGAGA_cuttrim_sorted.bam BJE3664_ThreeSis_ACAGGGA_cuttrim_sorted.bam BJE3665_ThreeSis_ACGTGGTA_cuttrim_sorted.bam BJE3666_ThreeSis_CCATGGGT_cuttrim_sorted.bam BJE3667_Citrus_CGCTT_cuttrim_sorted.bam BJE3668_Citrus_TCACG_cuttrim_sorted.bam BJE3669_Citrus_CTAGG_cuttrim_sorted.bam BJE3670_Citrus_ACAAA_cuttrim_sorted.bam BJE3671_Citrus_TTCTG_cuttrim_sorted.bam BJE3672_Citrus_AGCCG_cuttrim_sorted.bam BJE3673_Citrus_GTATT_cuttrim_sorted.bam BJE3674_Citrus_CTGTA_cuttrim_sorted.bam BJE3675_Citrus_ACCGT_cuttrim_sorted.bam BJE3676_Citrus_GCTTA_cuttrim_sorted.bam BJE3677_Citrus_GGTGT_cuttrim_sorted.bam BJE3678_Citrus_AGGAT_cuttrim_sorted.bam JM_no_label1_Draken_CCACGT_cuttrim_sorted.bam JM_no_label2_Draken_TTCAGA_cuttrim_sorted.bam Vred_8_Vred_GCTGTGGA_cuttrim_sorted.bam | bcftools call -V indels --format-fields GQ -m -O z -O z -o laevis_GBS_2020.vcf.gz

gunzip laevis_GBS_2020.vcf.gz

vcf-query -l laevis_GBS_2020.vcf > all_sample_list

vcftools --gzvcf laevis_GBS_2020.vcf --keep all_sample_list --chr chr1L --chr chr1S --chr chr2L --chr chr2S --chr chr3L --chr chr3S --chr chr4L --chr chr4S --chr chr5L --chr chr5S --chr chr6L --chr chr6S --chr chr7L --chr chr7S --chr chr8L --chr chr8S --chr chr9_10L --chr chr9_10S --non-ref-ac-any 1 --recode --recode-INFO-all --stdout | grep -v '^##contig=<ID=Scaffold' > all_samples_varient_chrs_only_header_scaffolds_cleared.vcf

vcftools --vcf ./all_samples_varient_chrs_only_header_scaffolds_cleared.vcf --keep all_sample_list --chr chr1L  --chr chr2L --chr chr3L --chr chr4L --chr chr5L  --chr chr6L  --chr chr7L  --chr chr8L  --chr chr9_10L  --recode --recode-INFO-all --stdout > L_only_varient_chrs_only_header_scaffolds_cleared.vcf

vcftools --vcf ./all_samples_varient_chrs_only_header_scaffolds_cleared.vcf --keep all_sample_list --chr chr1S  --chr chr2S --chr chr3S --chr chr4S --chr chr5S  --chr chr6S  --chr chr7S  --chr chr8S  --chr chr9_10S  --recode --recode-INFO-all --stdout > S_only_varient_chrs_only_header_scaffolds_cleared.vcf

```

Then seperated these VCF files into seperate directories for subgenomes and whole genome



