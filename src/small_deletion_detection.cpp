/* C++ header files */
#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <string>
#include <algorithm>

/* C header files */
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include <zlib.h>

/* custom header files */
#include "tk.h"
#include "cnv.h"


int output_cnvcall_vector(const CHR_INFO * chr_info, const std::vector <CNVCall> & cnv_call_vector, const std::string & out_file)
{
    FILE * out_fp;
    out_fp = fopen(out_file.c_str(), "a");
    if (NULL == out_fp)
    {
        fprintf(stderr, "ERROR! Failed to open file: %s\n", out_file.c_str());
        exit(1);
    }
    for (size_t i = 0; i < cnv_call_vector.size(); i++)
    {
        if (cnv_call_vector[i].cnv_type != 0 && cnv_call_vector[i].score[cnv_call_vector[i].cnv_hap_type] >= min_pass_score)
        {
            cnv_call_vector[i].output_bedpe_line(out_fp, chr_info);
            fprintf(out_fp, "\n");
        }
    }
    fclose(out_fp);
    return 0;
}

int get_depth_cutoff_for1hap(const QuantileNumbers & depth_quantiles, float & upper_depth_cutoff, float & lower_depth_cutoff, float &upper_depth_pvalue, float & lower_depth_pvalue)
{
    float mean_depth;

    mean_depth = depth_quantiles.q[500];
    lower_depth_cutoff = mean_depth / 10.0;
    upper_depth_cutoff = mean_depth * 1.5;

    upper_depth_pvalue = 1.0;
    lower_depth_pvalue = 1.0;

    for (size_t i = 0; i < depth_quantiles.q.size(); i++)
    {
        if (depth_quantiles.q[i] > lower_depth_cutoff)
        {
            lower_depth_pvalue = (float) depth_quantiles.b[i];
            break;
        }
    }

    for (size_t i = 0; i < depth_quantiles.q.size(); i++)
    {
        if (depth_quantiles.q[i] > upper_depth_cutoff)
        {
            upper_depth_pvalue = (float)1.0 - (float)depth_quantiles.b[i];
            break;
        }
    }

    return 0;
}

int get_depth_rate_cutoff_for1hap (const QuantileNumbers & depth_ratio_quantiles, float & lower_depth_ratio_cutoff, float & upper_depth_ratio_cutoff, float & lower_depth_ratio_pvalue, float & upper_depth_ratio_pvalue)
{
    lower_depth_ratio_pvalue = 0.01;
    upper_depth_ratio_pvalue = 0.01;
    lower_depth_ratio_cutoff = depth_ratio_quantiles.q[010];
    upper_depth_ratio_cutoff = depth_ratio_quantiles.q[990];
    return 0;
}

int calculate_depth_rate_distributions(const std::vector <std::vector <CNVInterval> > & wg_2d_interval_vector, std::vector <QuantileNumbers> & out_depth_rate_quantiles)
{
    int frm_depth;
    float read_depth;
    float depth_rate;
    double vector_size;
    size_t idx; 

    std::vector <float> depth_rate_vector;

    for (int hap_type = 0; hap_type < 4; hap_type ++)
    {
        depth_rate_vector.clear();
        for (int tid = 0; tid < wg_2d_interval_vector.size(); tid++)
        {
            for (int idx = 0; idx < wg_2d_interval_vector[tid].size(); idx ++)
            {
                read_depth = wg_2d_interval_vector[tid][idx].depth[hap_type];
                frm_depth = wg_2d_interval_vector[tid][idx].frm_depth[hap_type];
                if (frm_depth == 0)  { continue;  }
                if (read_depth < 0.001) { continue; }
        
                depth_rate = read_depth / frm_depth;
                depth_rate_vector.push_back(depth_rate);
            }
        }
        std::sort (depth_rate_vector.begin(), depth_rate_vector.end());
        vector_size = depth_rate_vector.size();

        for (int i = 0; i < 1001; i++)
        {
            idx = (int)(vector_size * out_depth_rate_quantiles[hap_type].b[i]+0.5);
            out_depth_rate_quantiles[hap_type].q[i] = depth_rate_vector[idx];
        }
        fprintf(stderr, "\n\n\nhap_type =%d\n", hap_type);
        fprintf(stderr, "depth ratio quantiles:\n");
        out_depth_rate_quantiles[hap_type].print(stderr, 0);
    }

    return 0;
}


int determine_frm_hap_type(int hap0_cnt, int hap1_cnt, int hap2_cnt)
{
    int sum; 
    sum = hap0_cnt + hap1_cnt + hap2_cnt; 
    if (hap0_cnt * 2 > sum) { return 0; }
    if (hap1_cnt * 2 > sum) { return 1; }
    if (hap2_cnt * 2 > sum) { return 2; }
    return 0;
}

int calculate_high_qual_read_count(char * map_qual, int & frm_read_count, int & frm_high_qual_read_count)
{
    int i = 0; 
    frm_read_count = 0;
    frm_high_qual_read_count = 0; 

    while (1)
    {
        if (map_qual[i] == '0')
        {
            frm_read_count ++; 
        }else if (map_qual[i] == '1')
        {
            frm_read_count ++; 
            frm_high_qual_read_count++;
        }else{
            break;
        }
        i++;
    }

    return 0; 
}

int get_frm_depth_for_each_interval(std::string bcd22_file, std::vector <std::vector <CNVInterval> > & wg_2d_interval_vector, const CHR_INFO * chr_info, int bin_size)
{
    int frm_tid, frm_start, frm_end; 
    int frm_hap_type;
    int frm_read_count, frm_high_qual_read_count; 
    int hap0_cnt, hap1_cnt, hap2_cnt; 
    char * map_pos;
    char * map_qual;
    int start_idx, end_idx; 
    const int min_frm_length = 2000;
    int min_frm_high_qual_reads = 6;

    char * chrname;
    gzFile bcd22_fp;

    CNVInterval frm_itv; 
    int index;
    char * line;
    
    
    line = new char [LOCAL_LINE_MAX];
    map_pos = new char [LOCAL_LINE_MAX];
    map_qual = new char [LOCAL_LINE_MAX];

    int total_frm_cnt = 0;
    int short_frm_cnt = 0;
    int low_qual_frm_cnt = 0;

    bcd22_fp = gzopen(bcd22_file.c_str(), "r");
    if (Z_NULL == bcd22_fp) {
        fprintf(stderr, "ERROR! Failed to open file for reading: %s\n", bcd22_file.c_str());
        exit(1);
    }

    int line_cnt = 0; 

    fprintf(stderr, "processing file: %s\n", bcd22_file.c_str());
    while (gzgets(bcd22_fp, line, LOCAL_LINE_MAX))
    {
        line_cnt ++;
        if (error_level > 1 && line_cnt % 10000000 == 0)
        {
            fprintf(stderr, "read %d lines in file: %s\n", line_cnt, bcd22_file.c_str());
        }
        if (line[0] == '#') { continue;}

        sscanf(line, "%d\t%d\t%d\t%*s\t%*s\t%*s\t%*s\t%d\t%d\t%d\t%s\t%s\t%*s\n", &frm_tid, &frm_start, &frm_end, &hap0_cnt, &hap1_cnt, &hap2_cnt, map_pos, map_qual);
        total_frm_cnt++;

        if (frm_end - frm_start < min_frm_length) { short_frm_cnt ++; continue; }

        calculate_high_qual_read_count(map_qual, frm_read_count, frm_high_qual_read_count); 

        if (frm_high_qual_read_count < min_frm_high_qual_reads || frm_high_qual_read_count * 2 < frm_read_count )  { low_qual_frm_cnt++; continue; }

        frm_hap_type = determine_frm_hap_type(hap0_cnt, hap1_cnt, hap2_cnt);
        start_idx = frm_start / bin_size;
        end_idx = frm_end / bin_size;
        for (int idx = start_idx; idx < end_idx; idx++)
        {
            if (wg_2d_interval_vector[frm_tid][idx].frm_depth[frm_hap_type] < MAX_FRM_DEPTH)
            {
                wg_2d_interval_vector[frm_tid][idx].frm_depth[frm_hap_type]++;
            }
        }
    }

 
    CNVInterval * ptr_itv;
    for (int tid = 0; tid < wg_2d_interval_vector.size(); tid++)
    {
        for (int idx = 0; idx < wg_2d_interval_vector[tid].size(); idx ++)
        {
            ptr_itv = &wg_2d_interval_vector[tid][idx];
            ptr_itv->frm_depth[3] = ptr_itv->frm_depth[0] + ptr_itv->frm_depth[1] + ptr_itv->frm_depth[2];
        }
    }

    delete [] line;
    delete [] map_pos;
    delete [] map_qual;
    gzclose(bcd22_fp);

    return 0; 
}


int read_weird_reads_cluster_file(const CHR_INFO * chr_info, const std::string & weird_reads_cluster_file, std::vector <CNVCall> & weird_reads_del_vector)
{
    FILE * weird_reads_cluster_fp;
    char * line;
    const size_t line_max = 1048576;
    std::vector <std::string > col_str_vector; 
    

    weird_reads_cluster_fp = fopen(weird_reads_cluster_file.c_str(), "r");
    if (NULL == weird_reads_cluster_fp) {
        fprintf(stderr, "Failed to open file: %s\n", weird_reads_cluster_file.c_str());
        exit(1);
    }

    line = new char[line_max];

    while (fgets(line, line_max, weird_reads_cluster_fp))
    {
        if (line[0] == '#') { continue; }
        CNVCall temp_cnv_call;
        col_str_vector = split_cstring_with_delimiter(line, "\t");
        if (col_str_vector.size() < 14)
        {
            fprintf(stderr, "ERROR! %s should have 14 columns!\n", weird_reads_cluster_file.c_str());
            exit(1);
        }
        temp_cnv_call.tid = chrname2tid(col_str_vector[0].c_str(), chr_info);
        temp_cnv_call.start_pos = std::stoi(col_str_vector[1]);
        temp_cnv_call.end_pos   = std::stoi(col_str_vector[4]);
        temp_cnv_call.num_pe_supp[3] = std::stoi(col_str_vector[9]);
        temp_cnv_call.num_pe_supp[0] = std::stoi(col_str_vector[10]);
        temp_cnv_call.num_pe_supp[1] = std::stoi(col_str_vector[11]);
        temp_cnv_call.num_pe_supp[2] = std::stoi(col_str_vector[12]);
        temp_cnv_call.aux_info = col_str_vector[13];
        temp_cnv_call.cnv_type = -1;
        weird_reads_del_vector.push_back(temp_cnv_call);
    }

    fclose(weird_reads_cluster_fp);
    delete [] line;

    return 0;
}

int validate_discordant_reads(CNVSettings & global_settings, const std::string & weird_reads_cluster_file, const std::vector <std::vector <CNVInterval> > & wg_2d_interval_vector, std::vector <CNVCall> & weird_reads_del_vector, const std::string & out_file)
{

    int bin_size;
    int start_idx, end_idx;
    float depth;
    float depth_rate;
    float frm_depth;

    float wg_mean_depth[4];
    float wg_mean_frm_depth[4];
    float wg_mean_depth_rate[4];
    float bin_cnt;

    for (int hap_type = 0; hap_type < 4; hap_type++)
    {
        wg_mean_depth[hap_type] = global_settings.depth_quantiles_vector[hap_type].q[500];
        wg_mean_depth_rate[hap_type] = global_settings.depth_rate_quantiles_vector[hap_type].q[500];
        if (wg_mean_depth_rate[hap_type] > 0)
        {
            wg_mean_frm_depth[hap_type] =  wg_mean_depth[hap_type] / wg_mean_depth_rate[hap_type];
            fprintf(stderr, "wg_mean_frm_depth[%d] = %.2f\n", hap_type, wg_mean_frm_depth[hap_type]);
        }else{
            fprintf(stderr, "ERROR! wg_mean_depth_rate[%d] = %f\n", hap_type, wg_mean_depth_rate[hap_type]);
            exit(1);
        }
    }
    
    fprintf(stderr, "reading file: %s\n", weird_reads_cluster_file.c_str());
    read_weird_reads_cluster_file(global_settings.chr_info, weird_reads_cluster_file, weird_reads_del_vector);

    if (wg_2d_interval_vector.size() == 0 || wg_2d_interval_vector[0].size() == 0)
    {
        fprintf(stderr, "ERROR! wg_2d_interval_vector is empty!\n");
        exit(1);
    }

    bin_size = wg_2d_interval_vector[0][0].end_pos - wg_2d_interval_vector[0][0].start_pos;
    for (auto & del_call : weird_reads_del_vector)
    {
        start_idx = del_call.start_pos / bin_size;
        end_idx = del_call.end_pos / bin_size;
        del_call.cnv_type = 0;      // means cnv not called
        del_call.cnv_hap_type = -1; // invalid value
        for (int hap_type = 0; hap_type < 4; hap_type++)
        {
            del_call.num_cnv_bin[hap_type] = 0;
            del_call.depth[hap_type] = 0.0;
            del_call.high_qual_depth[hap_type] = 0.0;
            del_call.frm_depth[hap_type] = 0;
            del_call.score[hap_type] = 0.0;
        }
        
        for (int hap_type = 0; hap_type < 4; hap_type++)
        {
            for (int idx = start_idx; idx < end_idx; idx++)
            {
                depth = wg_2d_interval_vector[del_call.tid][idx].depth[hap_type];
                frm_depth = wg_2d_interval_vector[del_call.tid][idx].frm_depth[hap_type];
                del_call.depth[hap_type] += depth;
                del_call.frm_depth[hap_type] += frm_depth;
                del_call.high_qual_depth[hap_type] += wg_2d_interval_vector[del_call.tid][idx].high_qual_depth[hap_type];
                
                if (depth < global_settings.lower_depth_cutoff[hap_type] && frm_depth > wg_mean_frm_depth[hap_type] / 4.0)
                {
                    del_call.num_cnv_bin[hap_type] ++;
                }
            }
            bin_cnt = end_idx - start_idx;
            del_call.depth[hap_type] /= bin_cnt;
            del_call.high_qual_depth[hap_type] /= bin_cnt;
            del_call.frm_depth[hap_type] /= bin_cnt;

            if ( (double) del_call.num_cnv_bin[hap_type]  > 0.667 * (double) (end_idx - start_idx) && del_call.num_pe_supp[hap_type] >= 2)
            {
                del_call.score[hap_type] = del_call.num_pe_supp[hap_type] * 2 + del_call.num_cnv_bin[hap_type] * 1;
            }else{
                del_call.score[hap_type] = 0.0;
            }
            if (del_call.score[hap_type] > 20)
            {
                del_call.cnv_type = -1;
                if (del_call.cnv_hap_type == -1 || del_call.score[hap_type] > del_call.score[del_call.cnv_hap_type])
                {
                    del_call.cnv_hap_type = hap_type;
                }
            }
        }
    }

    output_cnvcall_vector(global_settings.chr_info, weird_reads_del_vector, out_file);

    return 0;
}



int call_small_deletions(const std::string hap_type_read_depth_file, const std::string weird_reads_cluster_file, const std::string bcd22_file, const std::string faidx_file, const std::string gap_region_bed_file, const std::string out_file)
{
    CNVSettings global_settings;
    std::vector <std::vector <CNVInterval> > wg_2d_interval_vector; 
    std::vector <CNVCall> wg_candidate_cnv_region_vector;
    std::vector <CNVCall> weird_reads_del_vector;

    global_settings.bin_size = 100;

    global_settings.chr_info = get_chr_info(faidx_file.c_str());
    
    generate_whole_genome_interval_vector(wg_2d_interval_vector, global_settings.chr_info, global_settings.bin_size);
  
    get_depth_for_each_interval(hap_type_read_depth_file, wg_2d_interval_vector, global_settings.chr_info, global_settings.bin_size);
   
    calculate_depth_distributions (wg_2d_interval_vector, global_settings.depth_quantiles_vector, global_settings.hiqh_qual_depth_quantiles_vector);

    global_settings.wg_mean_depth = global_settings.depth_quantiles_vector[3].q[500];
    global_settings.wg_sd_depth = global_settings.wg_mean_depth;

    if (error_level > 1)
    {
        fprintf(stderr, "whole genome average depth = %.2f\n", global_settings.wg_mean_depth);
    }

    get_frm_depth_for_each_interval(bcd22_file, wg_2d_interval_vector, global_settings.chr_info, global_settings.bin_size);

    calculate_depth_rate_distributions(wg_2d_interval_vector, global_settings.depth_rate_quantiles_vector);

    fprintf(stderr, "\n");
    for (int hap_type = 1; hap_type < 4; hap_type++)
    {
        get_depth_cutoff_for1hap (global_settings.depth_quantiles_vector[hap_type], global_settings.upper_depth_cutoff[hap_type], global_settings.lower_depth_cutoff[hap_type], global_settings.upper_depth_pvalue[hap_type], global_settings.lower_depth_pvalue[hap_type]);

        get_depth_rate_cutoff_for1hap (global_settings.depth_rate_quantiles_vector[hap_type], global_settings.lower_depth_ratio_cutoff[hap_type], global_settings.upper_depth_ratio_cutoff[hap_type], global_settings.lower_depth_ratio_pvalue[hap_type], global_settings.upper_depth_ratio_pvalue[hap_type]);

        if (error_level > 1)
        {
            fprintf(stderr, "lower_depth_cutoff[%d]=%f\n", hap_type, global_settings.lower_depth_cutoff[hap_type] );
            fprintf(stderr, "lower_depth_pvalue[%d]=%f\n\n", hap_type, global_settings.lower_depth_pvalue[hap_type] );

            fprintf(stderr, "upper_depth_cutoff[%d]=%f\n", hap_type, global_settings.upper_depth_cutoff[hap_type] );
            fprintf(stderr, "upper_depth_pvalue[%d]=%f\n\n", hap_type, global_settings.upper_depth_pvalue[hap_type] );

            fprintf(stderr, "upper_depth_ratio_cutoff[%d]=%f\n", hap_type, global_settings.lower_depth_ratio_cutoff[hap_type] );
            fprintf(stderr, "upper_depth_ratio_pvalue[%d]=%f\n\n", hap_type, global_settings.upper_depth_ratio_cutoff[hap_type] );
        }

    }

    validate_discordant_reads(global_settings, weird_reads_cluster_file, wg_2d_interval_vector, weird_reads_del_vector, out_file);

    return 0;    
}

int main (int argc, char * argv[])
{
    std::string usage; 
    std::string hap_type_read_depth_file, weird_reads_cluster_file, bcd22_file, faidx_file, gap_region_bed_file, out_file;
    double alpha;
    double min_LRR;
    int32_t min_size;

    usage = "Usage: call_small_deletions <in.hap_type_read_depth_file> <in_weird_reads_cluster_file> <in.bcd22> <faidx_file> <gap_region.bed> <out_file>\n";

    int min_arg_num = 0; 
    for (int i = 0; i < usage.size(); i++)
    {
        if (usage[i] == '<')
        {
            min_arg_num++; 
        }
    }
    
    if (argc < min_arg_num + 1){
        std::cerr << usage << std::endl;
        return 1;
    }

    hap_type_read_depth_file = argv[1];
    weird_reads_cluster_file = argv[2];
    bcd22_file               = argv[3];
    faidx_file               = argv[4];
    gap_region_bed_file      = argv[5];
    out_file                 = argv[6];

    FILE * fp;
    fp = fopen(out_file.c_str(), "w"); 
    fclose(fp);


    call_small_deletions(hap_type_read_depth_file, weird_reads_cluster_file, bcd22_file, faidx_file, gap_region_bed_file, out_file);

    return 0;
}
