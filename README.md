# laevis_GBS_2020

# Background
Received aligned and indexed data from Ben(.bam and .bai) Going to craete structure plots and cal Fst)
Copied those files into 

# Filter and prepare clear bam files (somehow sort at the end had not exactly done the job. So I had to sort and index again in next script)

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

samtools view -b ${i} chr1S chr2L chr2S chr3L chr3S chr4L chr4S chr5L chr5S chr6L chr6S chr7L chr7S chr8L chr8S chr9_10L chr9_10S > inbetween.bam

#index inbetween bam
samtools index inbetween.bam

# Remove indels
samtools view -h inbetween.bam | awk '$1 ~ "^@" || $6 !~ "I|D"' | samtools view -b > test.bam

# sort for final output
samtools sort test.bam -o ${i}_final_sorted.bam

# Index the final output file

samtools index ${i}_final_sorted.bam;done

# Seperate L subgenome
for i in *_final_sorted.bam;do samtools view -b $i chr1L chr2L chr3L chr4L chr5L chr6L chr7L chr8L chr9_10L > ${i}_L_only; done

# Seperate S subgenome
for i in *_final_sorted.bam;do samtools view -b $i chr1S chr2S chr3S chr4S chr5S chr6S chr7S chr8S chr9_10S > ${i}_S_only; done

# Make directories for bamfiles for different genomes
mkdir ../filtered_bam_files

mkdir ../filtered_bam_files/whole_genome
mkdir ../filtered_bam_files/s_only
mkdir ../filtered_bam_files/l_only

# Move files to corresponding directories
mv *_final_sorted.bam ../filtered_bam_files/whole_genome
mv *_L_only ../filtered_bam_files/l_only
mv *_S_only ../filtered_bam_files/s_only


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

** Processing all files this way might take lot of time. In that case, modify the outer for loop and use this script to submit jobs parallely for seperate genomes(leaving just one genome as j at a time)

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
for i in *bam_final_sorted.bam*; do samtools sort ${i} -o ${i%.bam_final_sorted.bam*}_final_${j#../}.bam
        samtools index ${i%.bam_final_sorted.bam*}_final_${j#../}.bam;done


# Create VCF
bcftools mpileup -a FORMAT/DP -q20 -d8000 -f ../reference_genome/XENLA_9.2_genome.fa *_final_${j#../}.bam | bcftools call -V indels --format-fields GQ -m -O z -O z -o laevis_GBS_2020_${j#../}.vcf.gz

# unzip VCF
gunzip laevis_GBS_2020_${j#../}.vcf.gz

# Get sample list
vcf-query -l laevis_GBS_2020_${j#../}.vcf > all_sample_list

# Select chromosomes and remove scaffolds, remove scaffold from header (Selected chromosomes here moght not affect subgenomes here as bam files used are already filtered for subgenomes)

vcftools --gzvcf laevis_GBS_2020_${j#../}.vcf --keep all_sample_list --chr chr1L --chr chr1S --chr chr2L --chr chr2S --chr chr3L --chr chr3S --chr chr4L --chr chr4S --chr chr5L --chr chr5S --chr chr6L --chr chr6S --chr chr7L --chr chr7S --chr chr8L --chr chr8S --chr chr9_10L --chr chr9_10S --non-ref-ac-any 1 --recode --recode-INFO-all --stdout | grep -v '^##contig=<ID=Scaffold' > laevis_GBS_2020_${j#../}_scaffolds_removed.vcf


# create a folder for older files
mkdir older_files

#move older files there 

mv *.bam_final_sorted.bam* ./older_files

# make a folder for the VCF 
mkdir mkdir ../../filtered_VCFs/vcf_${j#../}

#move vcf there
mv laevis_GBS_2020_${j#../}* ../../filtered_VCFs/vcf_${j#../}/ ;done

```
## Run following in one of the folders created for filtered bams by the first script to collect finalized bams in a seperate folder(filtered_bam_files/l_only, PNLY IF NEEDED)

 make a directory to collect all finalized bam files
 ```bash
mkdir ../../finalized_bam_files
```

 Collect all finalized bam files looping through all directories and place them in respective folders
```bash
for j in  ../l_only ../s_only ../whole_genome; do cd ${j}
        mkdir ../../finalized_bam_files/finalized_bams_${j#../}
        mv ./*final_${j#../}.bam* ../../finalized_bam_files/finalized_bams_${j#../} ; done
```
# Filtering based on per site coverage

```bash
module load nixpkgs/16.09
module load gatk/4.1.2.0
```
you may have to type “yes” here to access GATK
Index reference genome (Use this if it has not been done already
```bash
module load bwa
module load samtools/1.10
samtools faidx XENLA_9.2_genome.fa
```
create GATk dictionary file for reference genome
```bsh
module load nixpkgs/16.09
module load gatk/4.1.2.0
gatk --java-options "-Xmx2G" CreateSequenceDictionary -R   XENLA_9.2_genome.fa
```
 Create depth table using GATK for all genomes(in /filtered_VCFs/vcf_l_only)
 
 ** You will not be able to get enough memory for java(like -Xmx16G) if you use bash to run script(max is 1/4 th of the physical memory you ask for). Therefore to test, submit job like this and keep checking .err file **
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

for j in ../*; do cd ${j}; java -Xmx16G -jar $EBROOTGATK/GenomeAnalysisTK.jar -T VariantsToTable -R ../../../reference_genome/XENLA_9.2_genome.fa -V *_scaffolds_removed.vcf -F CHROM -F POS -GF DP -o laevis_GBS_2020_${j#../}_GVCF_DP.table ; done
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

then,

```bash





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



