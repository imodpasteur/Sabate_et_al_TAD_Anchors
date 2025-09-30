# Micro-C analysis

Code to analyze raw sequencing reads from Capture Micro-C and to detect and score TADs and loops in the GSE104333 dataset from Rao et al, Cell, 2017.<br>


## Requirements:

1. Juicer 1.19.02
2. HiCPro 3.1.0
3. FastQC v0.12.1
4. Cutadapt v0.6.10
5. Bowtie2 v2.4.4
6. Cooler v0.9.1


## Micro Capture-C sequencing processing

Raw sequencing reads from Capture Micro-C were analyzed by Dr. Michael Szalay

### Raw sequencing data quality check by FastQC command:
```
fastqc seqfile1.fq.gz
```
### Trimming to 50 bp using TrimGalore command:
```
trimgalore --hardtrim5 50
```
### HiC-Pro pipeline, Bowtie2 settings:
```
HiC-Pro -i /micro-c/samples/ -o /micro-c/hic_out_samples/ -c HiC-Pro_3.1.0/config-hicpro.txt -s mapping -s proc_hic -s quality_checks -s merge_persample
#BOWTIE2_GLOBAL_OPTIONS = --very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder
#BOWTIE2_LOCAL_OPTIONS =  --very-sensitive -L 20 --score-min L,-0.6,-0.2 --end-to-end --reorder
```
### Cooler command
```
cooler cload pairs -c1 2 -p1 3 -c2 5 -p2 6 ./scripts/chrom_sizes.txt:100
```
### .mcool files command
```
cooler zoomify -r 100, 400, 1000, 2000, 4000, 10000, 20000, 50000, 100000, 1000000 file.cool -o file.mcool 
```

## Hi-C analysis of Rao et al, Cell, 2017
### Loop detection and scoring using Hiccups from Juicer
```
java -Xms8000m -Xmx27000m -Djava.library.path=/PathTo/JCuda-All-10.0.0 -jar juicer_tools_1.19.02.jar hiccups -m 2048 -r 5000,10000 -k KR -f 0.1 -p 4,2 -i 7,5 -t 0.02,1.5,1.75,2 -d 20000,20000 --threads 18 GSE104333_Rao-2017-untreated_combined.hic hiccups_untreated_5kb
```
### TAD detection and scoring using Arrowhead from Juicer
```
java -jar juicer_tools_1.19.02.jar arrowhead -m 2000 -r 5000 -k KR --threads 10 GSE104333_Rao-2017-untreated_combined.hic TADs_5000_untreated 
```
