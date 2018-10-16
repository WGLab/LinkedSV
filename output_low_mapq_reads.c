#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include "utils.h"

#include <htslib/sam.h>

#define UNKNOWN_HAP_TYPE 0 
#define LINE_MAX 4096
#define N_SPLIT_FILE 1


int output_low_mapq_reads(const char * in_bam, const char * out_bam, int map_qual_cutoff)
{

    samFile * in, *out;
    bam_hdr_t * hdr;
    bam1_t * b;
    bam1_core_t * c;
    int ret;

    b = bam_init1();
    c = &b->core;

    in = sam_open(in_bam, "rb");
    if (NULL == in) {
        fprintf (stderr, "ERROR: failed to open file: %s\n", in_bam);
        abort();
    }
    out = sam_open(out_bam, "w");
    if (NULL == out) {
        fprintf (stderr, "ERROR: failed to open file: %s\n", out_bam);
        abort();
    }

    hdr = sam_hdr_read(in);
    sam_hdr_write(out, hdr);

    while (ret = sam_read1(in, hdr, b) >= 0)
    {
        if (c->flag & BAM_FUNMAP) {  // ignore unmapped reads
            continue;
        }
        if (c->flag & BAM_FSUPPLEMENTARY) { // ignore supplementary alignments
            continue;
        }   

        if (c->flag & BAM_FSECONDARY) {   // ingnore secondary alignments
            continue;
        }   

        if (c->flag & BAM_FDUP) {  // ignore duplications
            continue;
        }

        if (c->qual < map_qual_cutoff){
            sam_write1(out, hdr, b);
        }

    }
    sam_close (in);
    sam_close (out);
    return 0;
}

int usage(FILE * fp)
{
    fprintf (fp, "Usage: output_low_mapq_reads <input_bam> <output_sam> <map_qual_cut_off>\n");
    return 0;
}

int main (int argc, char * argv[])
{
    char * in_bam;
    char * out_bam;
    int map_qual_cutoff;

    int c;

    in_bam        = NULL;
    out_bam     = NULL;

    if (argc < 4){
        usage(stderr);
        return 1;
    }

    in_bam        = argv[1];
    out_bam     = argv[2];
    map_qual_cutoff = atoi(argv[3]);

    if (NULL == in_bam){
        fprintf (stderr, "ERROR! input_bam is required!\n");
        usage(stderr);
        return 1;
    }
    if (NULL == out_bam){
        fprintf (stderr, "ERROR! output_file is required!\n");
        usage(stderr);
        return 1;
    }
    if (map_qual_cutoff < 0 || map_qual_cutoff > 30){
        fprintf (stderr, "ERROR! map_qual_cutoff should be in the range of [0, 30].\n");
        usage(stderr);
    }
    output_low_mapq_reads(in_bam, out_bam, map_qual_cutoff);

    return 0;    
}
