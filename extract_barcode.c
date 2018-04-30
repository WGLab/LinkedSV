#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <utils.h>

#include <htslib/sam.h>

#define UNKNOWN_HAP_TYPE 0 
#define LINE_MAX 4096

int extract_barcode_info (const char * inbam, const char * out_file, char * stat_file, int min_map_qual)
{
    samFile * in;
    bam_hdr_t * hdr;
    bam1_t * b;
    bam1_core_t * c;
    int ret;
    int flag_filter;
    FILE * out_fp;
    FILE * stat_fp;
    int no_bx_cnt;
    int no_hp_cnt;
    int tot_rd_cnt;
    int output_rd_cnt;
    int mapq10_cnt;
    int mapq20_cnt;
    int mapq30_cnt;
    int unmapped_rd_cnt;
    int suppl_align_cnt;
    int secdn_align_cnt;
    int duplication_cnt;
    int32_t pair_tid;
    int32_t pair_pos;
    int64_t pair_abs_pos;

    char * barcode = NULL;
    int32_t hap_type;
    uint8_t * s;
    uint32_t bcd_int;


    b = bam_init1();
    c = &b->core;
    out_fp = NULL;

    in = sam_open(inbam, "rb");
    if (NULL == in) {
        fprintf (stderr, "ERROR: failed to open file: %s\n", inbam);
        abort();
    }
    // to be added: check sort bam //


    stat_fp = fopen(stat_file, "w");
    if (NULL == stat_fp){
        fprintf (stderr, "ERROR: failed to open file for writing: %s\n", stat_file);
        abort();
    }


    hdr = sam_hdr_read(in);
    tot_rd_cnt = output_rd_cnt = 0;
    no_bx_cnt = no_hp_cnt = 0;
    unmapped_rd_cnt = duplication_cnt = secdn_align_cnt = suppl_align_cnt = 0;
    mapq10_cnt = mapq20_cnt = mapq30_cnt = 0;

    out_fp  = fopen (out_file, "w");
    if (NULL == out_fp){
        fprintf (stderr, "ERROR: failed to open file for writing: %s\n", out_file);
        abort();
    }
    while (ret = sam_read1(in, hdr, b) >= 0)
    {
        tot_rd_cnt++;
        if (c->flag & BAM_FUNMAP) {  // ignore unmapped reads
            unmapped_rd_cnt ++;
            continue;
        }
        if (c->flag & BAM_FSUPPLEMENTARY) { // ignore supplementary alignments
            suppl_align_cnt ++;
            continue;
        }

        if (c->flag & BAM_FSECONDARY) {   // ingnore secondary alignments
            secdn_align_cnt++;
            continue;
        }

        if (c->flag & BAM_FDUP) {  // ignore duplications
            duplication_cnt++;
            continue;
        }

        s = bam_aux_get(b, "BX");
        if (NULL == s){
            no_bx_cnt ++; // reads without barcode 
            continue;
        }

        barcode = bam_aux2Z(s);

        s = bam_aux_get(b, "HP");

        if (NULL == s){
            no_hp_cnt++; // reads without haplotype 
            hap_type = UNKNOWN_HAP_TYPE;
        }else{
            hap_type = bam_aux2i(s);
        }
        
        if (c->qual >= 30){
            mapq30_cnt++;
        }else if (c->qual >= 20){
            mapq20_cnt++;
        }else if (c->qual >= 10){
            mapq10_cnt++;
        }

        if (c->qual < min_map_qual){
            continue;
        }
        
        fprintf(out_fp,"%d\t%d\t%d\t%d\t%s\t%d\n", c->tid, c->pos, bam_endpos(b), c->qual, barcode, hap_type);
        output_rd_cnt++;

    }

    fprintf (stat_fp, "####statistics####\n\n");
    fprintf (stat_fp, "bam file:                     %s\n", inbam);
    fprintf (stat_fp, "total alignment count:        %d\n", tot_rd_cnt);
    fprintf (stat_fp, "unmapped read count:          %d (%.2f%%)\n", unmapped_rd_cnt, (double)unmapped_rd_cnt/(double)tot_rd_cnt*100.0);
    fprintf (stat_fp, "supplementary alignments:     %d (%.2f%%)\n", suppl_align_cnt, (double)suppl_align_cnt/(double)tot_rd_cnt*100.0);
    fprintf (stat_fp, "secondary alignments:         %d (%.2f%%)\n", secdn_align_cnt, (double)secdn_align_cnt/(double)tot_rd_cnt*100.0);
    fprintf (stat_fp, "duplications:                 %d (%.2f%%)\n", duplication_cnt, (double)duplication_cnt/(double)tot_rd_cnt*100.0);
    fprintf (stat_fp, "\n");
    fprintf (stat_fp, "output read count:            %d\n", output_rd_cnt);
    fprintf (stat_fp, "reads with MapQ >=30:         %d (%.2f%%)\n", mapq30_cnt, (double)mapq30_cnt/(double)output_rd_cnt * 100.0);
    fprintf (stat_fp, "reads with MapQ in [20,30):   %d (%.2f%%)\n", mapq20_cnt, (double)mapq20_cnt/(double)output_rd_cnt * 100.0);
    fprintf (stat_fp, "reads with MapQ in [10,20):   %d (%.2f%%)\n", mapq10_cnt, (double)mapq10_cnt/(double)output_rd_cnt * 100.0);
    fprintf (stat_fp, "reads without barcode info:   %d (%.2f%%)\n", no_bx_cnt, (double)no_bx_cnt/(double)output_rd_cnt * 100.0);
    fprintf (stat_fp, "reads without phasing info:   %d (%.2f%%)\n", no_hp_cnt, (double)no_hp_cnt/(double)output_rd_cnt * 100.0);

    sam_close (in);
    fclose(out_fp);
    return 0;
}

int usage(FILE * fp)
{
    fprintf (fp, "Usage: extract_barcode <input_bam> <output_file> <statistics_file> <min_map_qual>\n");
    return 0;
}

int main (int argc, char * argv[])
{
    char * inbam;
    char * out_file;
    char * stat_file;
	int min_map_qual;
	
    int c;

    inbam        = NULL;
    out_file     = NULL;
    stat_file    = NULL;
    
    if (argc < 5){
        usage(stderr);
        return 1;
    }

    inbam        = argv[1];
    out_file      = argv[2];
	stat_file    = argv[3];
	min_map_qual = atoi(argv[4]);

    if (NULL == inbam){
        fprintf (stderr, "ERROR! input_bam is required!\n");
        usage(stderr);
        return 1;
    }
    if (NULL == out_file){
        fprintf (stderr, "ERROR! output_file is required!\n");
        usage(stderr);
        return 1;
    }
    if (NULL == stat_file){
        fprintf (stderr, "ERROR! statistics_file is required!\n");
        usage(stderr);
        return 1;
	}
	if (min_map_qual < 0 || min_map_qual > 30){
        fprintf (stderr, "ERROR! min_map_qual should be in the range of [0, 30].\n");
		usage(stderr);
	}
    extract_barcode_info(inbam, out_file, stat_file, min_map_qual);

    return 0;    
}
