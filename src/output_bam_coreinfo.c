#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>

#include <htslib/sam.h>

#define UNKNOWN_HAP_TYPE 0 
#define LINE_MAX 4096
#define N_SPLIT_FILE 1


int output_coreinfo_bam (const char * in_bam)
{

    samFile * in, *out;
    bam_hdr_t * hdr;
    bam1_t * b;
    bam1_core_t * c;
    int ret;
    uint8_t * seq, *aux,*p1, *p2;


    b = bam_init1();
    c = &b->core;

    in = sam_open(in_bam, "rb");
    if (NULL == in) {
        fprintf (stderr, "ERROR: failed to open file: %s\n", in_bam);
        abort();
    }

    hdr = sam_hdr_read(in);
    out = sam_open("-", "w");
    if (NULL == out) {
        fprintf (stderr, "ERROR: failed to write to stdout\n");
        abort();
    }
    sam_hdr_write(out, hdr);
    while (ret = sam_read1(in, hdr, b) >= 0)
    {
        if (c->flag & BAM_FUNMAP) {
            continue;
        }
        seq = bam_get_seq(b);
        aux = bam_get_aux(b);
        p1 = seq;
        p2 = aux;
        while (p2 < b->data + b->l_data)
        {
            *p1 = *p2;
            p1 ++; 
            p2 ++;
        }

        b->l_data -= aux - seq;
        c->l_qseq = 0;
        sam_write1(out, hdr, b);
    }

    sam_close (in);
    sam_close (out);
    return 0;
}

int usage(FILE * fp)
{
    fprintf (fp, "Usage: extract_barcode <input_bam>\n");
    return 0;
}

int main (int argc, char * argv[])
{
    char * in_bam;
	
    in_bam      = NULL;
    
    if (argc < 2){
        usage(stderr);
        return 1;
    }

    in_bam      = argv[1];

    if (NULL == in_bam){
        fprintf (stderr, "ERROR! input_bam is required!\n");
        usage(stderr);
        return 1;
    }

    output_coreinfo_bam(in_bam);

    return 0;    
}
