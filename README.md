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
Removed sample with a lower depth value

BJE3577_Kimber_TAACGA

# seperate subgenomes 

In the directory with all bam files
``bash
module load nixpkgs/16.09  intel/2016.4
module load samtools/0.1.17

for i in *.bam; do samtools index ${i}; done

for i in *.bam;do samtools view -b $i chr1L chr2L chr3L chr4L chr5L chr6L chr7L chr8L chr9_10L > ../bams_L_only/L_only_${i} ;done


for i in *.bam;do samtools view -b $i chr1S chr2S chr3S chr4S chr5S chr6S chr7S chr8S chr9_10S > ../bams_S_only/S_only_${i} ;done
```
# plot depth
```rscript



