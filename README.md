# LinkedSV



### Prerequisites
GCC that supports c99/gnu99. Succesfully tested on GCC 4.8.5. Higher vesion should be good. 

Python packages: psutil, subprocess, sklearn, scipy, numpy 

### Installation of LinkedSV
```
git clone https://github.com/WGLab/LinkedSV.git 
cd LinkedSV/htslib-1.3/
make
make check
cp libhts.a ../lib/
cd .. 
sh build.sh 
```

### Running LinkedSV

```
usage: linkedsv.py [-h] -i input.sorted.bam -d out_directory -r ref.fasta
                   [-v reference_version] [--gap_region_bed gap_region.bed]
                   [--black_region_bed black_region.bed] [-t num_thread]
                   [-q min_map_qual]
                   [--min_fragment_length min_fragment_length]
                   [--min_supp_barcodes min_supporting_barcodes]
                   [--samtools path/to/samtools] [--bedtools path/to/bedtools]
                   [--wgs] [--targeted] [--target_region target_region.bed]
                   [--gap_distance_cut_off gap_distance_cut_off]
```
Examples: 
`python linkedsv.py -i input.bam -d path/to/output/dir/ -r ref.fasta -v ref_version -t num_threads`

The ref.fasta file should be the same fasta file that was used for alignment. 
`ref_version` is used to tell LinkedSV which black_list file and gap_region file should be used. Currently LinkedSV only supports hg19 and b37, but we will soon support hg38. If `ref_version` is not specifed or you are using a different reference file, please generate these files by yourself and specify the `--gap_region_bed` and `--black_region_bed` parameters. 

