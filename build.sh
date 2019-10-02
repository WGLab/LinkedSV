#!/bin/bash
mkdir -p bin/

## pigz ##
cd pigz
make
mv pigz ../bin/
mv unpigz ../bin/
cd ../

## fermikit ##
cd fermikit
make
cd ..


## cpp ##
cd src

mkdir -p lib/
cd htslib-1.3/
./configure
make
make check
cp libhts.a ../lib/
cd .. 

g++ -g -O1 -std=c++11 -I ./include/ -L ./lib/ -o cluster_reads cluster_reads.cpp tk.cpp cgranges.cpp -l z -pthread
g++ -g -O1 -std=c++11 -I ./include/ -L ./lib/ -o extract_barcode_info extract_barcode_info.cpp  -l hts -l z -l m -l pthread
g++ -g -O1 -std=c++11 -I ./include/ -L ./lib/ -o output_bam_coreinfo output_bam_coreinfo.cpp -l hts -l z -l m -l pthread
g++ -g -O1 -std=c++11 -I ./include/ -L ./lib/ -o remove_sparse_nodes  remove_sparse_nodes.cpp  tk.cpp cgranges.cpp  -lz

g++ -g -O1 -std=c++11 -I ./include/ -L ./lib/ -o cal_hap_read_depth_from_bcd21 cal_hap_read_depth_from_bcd21.cpp  tk.cpp cgranges.cpp  -l z
g++ -g -O1 -std=c++11 -I ./include/ -L ./lib/ -o grid_overlap grid_overlap.cpp tk.cpp cgranges.cpp  -lz
g++ -g -O1 -std=c++11 -I ./include/ -L ./lib/ -o cal_read_depth_from_bcd21 cal_read_depth_from_bcd21.cpp tk.cpp cgranges.cpp  -lz
g++ -g -O1 -std=c++11 -I ./include/ -L ./lib/ -o cal_barcode_depth_from_bcd21 cal_barcode_depth_from_bcd21.cpp tk.cpp cgranges.cpp  -lz

g++ -g -O1 -std=c++11 -I ./include/ -L ./lib/ -o cal_twin_win_bcd_cnt cal_twin_win_bcd_cnt.cpp tk.cpp cgranges.cpp  -lz
g++ -g -O1 -std=c++11 -I ./include/ -L ./lib/ -o cal_centroid_from_read_depth cal_centroid_from_read_depth.cpp tk.cpp cgranges.cpp  -lz
g++ -g -O1 -std=c++11 -I ./include/ -L ./lib/ -o cal_2d_overlapping_barcodes cal_2d_overlapping_barcodes.cpp tk.cpp cgranges.cpp  -lz
g++ -g -O1 -std=c++11 -I ./include/ -L ./lib/ -o call_small_deletions call_small_deletions.cpp tk.cpp cgranges.cpp  -l z

mv cluster_reads extract_barcode_info output_bam_coreinfo remove_sparse_nodes ../bin/
mv cal_hap_read_depth_from_bcd21 grid_overlap cal_read_depth_from_bcd21 cal_barcode_depth_from_bcd21 ../bin/
mv cal_twin_win_bcd_cnt cal_centroid_from_read_depth cal_2d_overlapping_barcodes call_small_deletions ../bin/

