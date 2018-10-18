cd htslib-1.3/
mkdir -p lib/
make
make check
cp libhts.a ../lib/
cd .. 
gcc  -g -O2 -I ./include  -L ./lib/ -o extract_barcode_info extract_barcode_info.c  -l hts -l z -l m -l pthread 
gcc  -g -O2 -I ./include  -L ./lib/ -o output_bam_coreinfo  output_bam_coreinfo.c   -l hts -l z -l m -l pthread 
gcc  -g -O2 -std=gnu99              -o remove_sparse_nodes  remove_sparse_nodes.c  tk.c

