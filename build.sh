#!/bin/bash
cd src

mkdir -p lib/
cd htslib-1.3/
make
make check
cp libhts.a ../lib/
cd .. 

gcc -g -O0 -std=c99 -I ./include -L ./lib/ -o extract_barcode_info extract_barcode_info.c              -l hts -l z -l m -l pthread 
gcc -g -O0 -std=c99 -I ./include -L ./lib/ -o output_bam_coreinfo  output_bam_coreinfo.c               -l hts -l z -l m -l pthread 
gcc -g -O0 -std=c99                        -o remove_sparse_nodes  remove_sparse_nodes.c          tk.c
gcc -g -O0 -std=c99                        -o grid_overlap         grid_overlap.c                 tk.c
gcc -g -O0 -std=c99                        -o cal_read_depth_from_bcd21 cal_read_depth_from_bcd21.c tk.c       -l z
gcc -g -O0 -std=c99                        -o cal_barcode_depth_from_bcd21 cal_barcode_depth_from_bcd21.c tk.c -l z
gcc -g -O0 -std=c99                        -o cal_twin_win_bcd_cnt cal_twin_win_bcd_cnt.c         tk.c -l z

cp extract_barcode_info output_bam_coreinfo remove_sparse_nodes grid_overlap cal_read_depth_from_bcd21 cal_barcode_depth_from_bcd21 cal_twin_win_bcd_cnt ../bin/
