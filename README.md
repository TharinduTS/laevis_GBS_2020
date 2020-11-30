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



