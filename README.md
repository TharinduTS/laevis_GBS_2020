# laevis_GBS_2020

# Background
Received aligned and indexed data from Ben(.bam and .bai) Going to craete structure plots and cal Fst)
Copied those files into 

# Cal depth and move depth files to the folder
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

















# ****************************** Cal Fst ***************************

# Background
Selected the bam files with a good depth value from above depth check step and used those bam files from here onwards

# Create a VCF using selected bam files

Downloaded reference genome 
```bash
wget http://ftp.xenbase.org/pub/Genomics/JGI/Xenla9.2/XENLA_9.2_genome.fa.gz
```
Observed addresses for all bam files
```bash
find ../all_bam_files/ | tr "\n" " "
```

Then used following script to create VCF using observed file list

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

samtools mpileup -q20 -d8000 -ugf ../reference_genome/XENLA_9.2_genome.fa ../../bam_files/2014_Inhaca_10_Inhaca_ATATGT_cuttrim_sorted.bam ../../bam_files/2014_Inhaca_150_Inhaca_ATCGTA_cuttrim_sorted.bam ../../bam_files/2014_Inhaca_152_Inhaca_CATCGT_cuttrim_sorted.bam ../../bam_files/2014_Inhaca_24_Inhaca_CGCGGT_cuttrim_sorted.bam ../../bam_files/2014_Inhaca_38_Inhaca_CTATTA_cuttrim_sorted.bam ../../bam_files/2014_Inhaca_52_Inhaca_GCCAGT_cuttrim_sorted.bam ../../bam_files/2014_Inhaca_65_Inhaca_GGAAGA_cuttrim_sorted.bam ../../bam_files/946_Draken_TCGTT_cuttrim_sorted.bam ../../bam_files/993_Draken_GGTTGT_cuttrim_sorted.bam ../../bam_files/amnh17260_Nigeria_GTGAGGGT_cuttrim_sorted.bam ../../bam_files/BJE2897_Lendu_TTCCTGGA_cuttrim_sorted.bam ../../bam_files/BJE3252_Cameroon_TATCGGGA_cuttrim_sorted.bam ../../bam_files/BJE3508_DeDorn_ATTGA_cuttrim_sorted.bam ../../bam_files/BJE3509_DeDorn_CATCT_cuttrim_sorted.bam ../../bam_files/BJE3510_DeDorn_CCTAG_cuttrim_sorted.bam ../../bam_files/BJE3511_DeDorn_GAGGA_cuttrim_sorted.bam ../../bam_files/BJE3512_DeDorn_GGAAG_cuttrim_sorted.bam ../../bam_files/BJE3513_DeDorn_GTCAA_cuttrim_sorted.bam ../../bam_files/BJE3514_DeDorn_TAATA_cuttrim_sorted.bam ../../bam_files/BJE3515_DeDorn_TACAT_cuttrim_sorted.bam ../../bam_files/BJE3525_Laigns_GAATTCA_cuttrim_sorted.bam ../../bam_files/BJE3526_Laigns_GAACTTG_cuttrim_sorted.bam ../../bam_files/BJE3527_Laigns_GGACCTA_cuttrim_sorted.bam ../../bam_files/BJE3528_Laigns_GTCGATT_cuttrim_sorted.bam ../../bam_files/BJE3529_Laigns_AACGCCT_cuttrim_sorted.bam ../../bam_files/BJE3530_Laigns_AATATGG_cuttrim_sorted.bam ../../bam_files/BJE3531_Laigns_ACGTGTT_cuttrim_sorted.bam ../../bam_files/BJE3532_Laigns_ATTAATT_cuttrim_sorted.bam ../../bam_files/BJE3533_Laigns_ATTGGAT_cuttrim_sorted.bam ../../bam_files/BJE3534_BW_CTCG_cuttrim_sorted.bam ../../bam_files/BJE3535_BW_TGCA_cuttrim_sorted.bam ../../bam_files/BJE3536_BW_ACTA_cuttrim_sorted.bam ../../bam_files/BJE3537_BW_CAGA_cuttrim_sorted.bam ../../bam_files/BJE3538_BW_AACT_cuttrim_sorted.bam ../../bam_files/BJE3539_BW_GCGT_cuttrim_sorted.bam ../../bam_files/BJE3541_BW_CGAT_cuttrim_sorted.bam ../../bam_files/BJE3542_BW_GTAA_cuttrim_sorted.bam ../../bam_files/BJE3543_BW_AGCG_cuttrim_sorted.bam ../../bam_files/BJE3544_BW_GATG_cuttrim_sorted.bam ../../bam_files/BJE3545_BW_TCAG_cuttrim_sorted.bam ../../bam_files/BJE3546_BW_TGCGA_cuttrim_sorted.bam ../../bam_files/BJE3547_GRNP_TAGGAA_cuttrim_sorted.bam ../../bam_files/BJE3548_GRNP_GCTCTA_cuttrim_sorted.bam ../../bam_files/BJE3549_GRNP_CCACAA_cuttrim_sorted.bam ../../bam_files/BJE3550_GRNP_CTTCCA_cuttrim_sorted.bam ../../bam_files/BJE3551_GRNP_GAGATA_cuttrim_sorted.bam ../../bam_files/BJE3552_GRNP_ATGCCT_cuttrim_sorted.bam ../../bam_files/BJE3553_GRNP_AGTGGA_cuttrim_sorted.bam ../../bam_files/BJE3554_GRNP_ACCTAA_cuttrim_sorted.bam ../../bam_files/BJE3573_VicW_CGCGGAGA_cuttrim_sorted.bam ../../bam_files/BJE3574_VicW_CGTGTGGT_cuttrim_sorted.bam ../../bam_files/BJE3575_Kimber_GTACTT_cuttrim_sorted.bam ../../bam_files/BJE3576_Kimber_GTTGAA_cuttrim_sorted.bam ../../bam_files/BJE3578_Kimber_TGGCTA_cuttrim_sorted.bam ../../bam_files/BJE3579_Kimber_TATTTTT_cuttrim_sorted.bam ../../bam_files/BJE3580_Kimber_CTTGCTT_cuttrim_sorted.bam ../../bam_files/BJE3581_Kimber_ATGAAAG_cuttrim_sorted.bam ../../bam_files/BJE3582_Kimber_AAAAGTT_cuttrim_sorted.bam ../../bam_files/BJE3632_Niewou_CATAAGT_cuttrim_sorted.bam ../../bam_files/BJE3633_Niewou_CGCTGAT_cuttrim_sorted.bam ../../bam_files/BJE3640_Niewou_CGGTAGA_cuttrim_sorted.bam ../../bam_files/BJE3641_Niewou_CTACGGA_cuttrim_sorted.bam ../../bam_files/BJE3642_Niewou_GCGGAAT_cuttrim_sorted.bam ../../bam_files/BJE3644_Niewou_TAGCGGA_cuttrim_sorted.bam ../../bam_files/BJE3645_Niewou_TCGAAGA_cuttrim_sorted.bam ../../bam_files/BJE3647_Niewou_TCTGTGA_cuttrim_sorted.bam ../../bam_files/BJE3654_ThreeSis_TGCTGGA_cuttrim_sorted.bam ../../bam_files/BJE3655_ThreeSis_ACGACTAG_cuttrim_sorted.bam ../../bam_files/BJE3656_ThreeSis_TAGCATGG_cuttrim_sorted.bam ../../bam_files/BJE3657_ThreeSis_TAGGCCAT_cuttrim_sorted.bam ../../bam_files/BJE3658_ThreeSis_TGCAAGGA_cuttrim_sorted.bam ../../bam_files/BJE3659_ThreeSis_TGGTACGT_cuttrim_sorted.bam ../../bam_files/BJE3660_ThreeSis_TCTCAGTG_cuttrim_sorted.bam ../../bam_files/BJE3661_ThreeSis_CGCGATAT_cuttrim_sorted.bam ../../bam_files/BJE3662_ThreeSis_CGCCTTAT_cuttrim_sorted.bam ../../bam_files/BJE3663_ThreeSis_AACCGAGA_cuttrim_sorted.bam ../../bam_files/BJE3664_ThreeSis_ACAGGGA_cuttrim_sorted.bam ../../bam_files/BJE3665_ThreeSis_ACGTGGTA_cuttrim_sorted.bam ../../bam_files/BJE3666_ThreeSis_CCATGGGT_cuttrim_sorted.bam ../../bam_files/BJE3667_Citrus_CGCTT_cuttrim_sorted.bam ../../bam_files/BJE3668_Citrus_TCACG_cuttrim_sorted.bam ../../bam_files/BJE3669_Citrus_CTAGG_cuttrim_sorted.bam ../../bam_files/BJE3670_Citrus_ACAAA_cuttrim_sorted.bam ../../bam_files/BJE3671_Citrus_TTCTG_cuttrim_sorted.bam ../../bam_files/BJE3672_Citrus_AGCCG_cuttrim_sorted.bam ../../bam_files/BJE3673_Citrus_GTATT_cuttrim_sorted.bam ../../bam_files/BJE3674_Citrus_CTGTA_cuttrim_sorted.bam ../../bam_files/BJE3675_Citrus_ACCGT_cuttrim_sorted.bam ../../bam_files/BJE3676_Citrus_GCTTA_cuttrim_sorted.bam ../../bam_files/BJE3677_Citrus_GGTGT_cuttrim_sorted.bam ../../bam_files/BJE3678_Citrus_AGGAT_cuttrim_sorted.bam ../../bam_files/JM_no_label1_Draken_CCACGT_cuttrim_sorted.bam ../../bam_files/JM_no_label2_Draken_TTCAGA_cuttrim_sorted.bam ../../bam_files/RT5_Botsw_GGATTGGT_cuttrim_sorted.bam ../../bam_files/Vred_8_Vred_GCTGTGGA_cuttrim_sorted.bam | bcftools call -V indels --format-fields GQ -m -O z -O z -o laevis_GBS_2020.vcf
```



