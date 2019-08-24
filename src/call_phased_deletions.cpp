/* C++ header files */
#include <iostream>
#include <vector>
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


const int64_t FIX_LENGTH = (int64_t)(1e10);
const int LOCAL_LINE_MAX = 1<<20;
const int MAX_DEPTH      = 10000 * 100;

const int HAP_UNK     = 1 << 0;
const int HAP1_DEL    = 1 << 1;
const int HAP2_DEL    = 1 << 2;

const int HAP1_DUP    = 1 << (1+16);
const int HAP2_DUP    = 1 << (2+16);

const int flag_dup    = HAP1_DUP | HAP2_DUP;
const int flag_del    = HAP1_DEL | HAP2_DEL;

bool debug = true;

class Settings
{
public:
    std::vector <QuantileNumbers> depth_quantiles_vector;
    std::vector <QuantileNumbers> hiqh_qual_depth_quantiles_vector;

    std::vector <float> upper_depth_cutoff;
    std::vector <float> lower_depth_cutoff;
    std::vector <float> upper_depth_pvalue;
    std::vector <float> lower_depth_pvalue;
    CHR_INFO * chr_info;
    int bin_size;

    cgranges_t * black_list_cr_interval_tree;

    Settings()
    {
        for (size_t i = 0; i < 3; i++)
        {
            upper_depth_cutoff.push_back(0.0);
            lower_depth_cutoff.push_back(0.0);
            upper_depth_pvalue.push_back(0.0);
            lower_depth_pvalue.push_back(0.0);

            depth_quantiles_vector.push_back(QuantileNumbers());
            hiqh_qual_depth_quantiles_vector.push_back(QuantileNumbers());
        }
    }
};


class CNVInterval
{
public:
    int32_t tid;
    int32_t start_pos;
    int32_t end_pos;
    float depth[3] = {0.0}; // depth each hap_type (0, 1, 2). 0 for unknown
    float high_qual_depth[3] = {0.0}; // depth of reads with high map qual

    CNVInterval()
    {
        tid = -1;
        start_pos = -1;
        end_pos = -1;
        depth[0] = depth[1] = depth[2] = 0.0;
        high_qual_depth[0] = high_qual_depth[1] = high_qual_depth[2] = 0.0;
    };


    int length() const
    {
        return end_pos - start_pos;
    }
    Interval to_interval(CHR_INFO * chr_info) const
    {
        Interval itv;
        itv.tid = tid;
        itv.start_pos = start_pos;
        itv.end_pos = end_pos;
        itv.ctg = chr_info->chrname_list->data_list[tid];
        return itv;
    }

    int output_line(FILE * out_fp, const CHR_INFO * const chr_info) const
    {
        char * chrname = chr_info->chrname_list->data_list[tid];
        fprintf(out_fp, "%s\t%d\t%d\t%d\t", chrname, start_pos, end_pos, length()); 
        for (int hap_type = 0; hap_type < 3; hap_type++ )
        {
            fprintf(out_fp, "\t%.2f\t%.2f", depth[hap_type], high_qual_depth[hap_type]);
        }
        return 0; 
    }   
};


class CNVCall : public CNVInterval 
{
public:
    double  score[3] = {0}; // variant score of the CNV call. 0, 1, and 2 for hap_type 0, 1 and 2. 
    uint32_t flag = 0;
    int cnv_hap_type; 
    int cnv_type; // 1 for dup, -1 for del
    CNVCall()
    {
        tid = -1;
        start_pos = -1;
        end_pos = -1;
        depth[0] = depth[1] = depth[2] = 0.0;
        high_qual_depth[0] = high_qual_depth[1] = high_qual_depth[2] = 0.0;

        flag = 0;
        score[0] = score[1] = score[2] = 0.0;
    }
    
    int output_line(FILE * out_fp, const CHR_INFO * const chr_info) const
    {
        char * chrname = chr_info->chrname_list->data_list[tid];
        fprintf(out_fp, "%s\t%d\t%d\t%d\t", chrname, start_pos, end_pos, length()); 
        fprintf(out_fp, "%u\t%d\t%d\t%.2f\t%.2f\t%f", flag, cnv_type, cnv_hap_type, depth[cnv_hap_type], high_qual_depth[cnv_hap_type], score[cnv_hap_type]);
        for (int hap_type = 0; hap_type < 3; hap_type++ )
        {
            fprintf(out_fp, "\t%.2f\t%.2f\t%f", depth[hap_type], high_qual_depth[hap_type], score[hap_type]);
        }
        return 0; 
    }   
    

};

int calculate_interval_avg_depth(const std::vector <std::vector <CNVInterval> > wg_2d_interval_vector, int tid, int start_pos, int end_pos, int hap_type, float & out_depth, float & out_high_qual_depth)
{
    // average depth of a interval
    float total_depth;
    float total_high_qual_depth;

    int bin_size;
    int start_idx, end_idx;

    if (wg_2d_interval_vector.size() == 0 || wg_2d_interval_vector[0].size() == 0)
    {
        fprintf(stderr, "ERROR! wg_2d_interval_vector is empty!\n");
        exit(1);
    }

    if (end_pos <= start_pos)
    {
        fprintf(stderr, "ERROR! end_pos is smaller than start_pos!\n");
        exit(1);
    }

    if (hap_type < 0 || hap_type > 2)
    {
        fprintf(stderr, "ERROR! hap_type should be 0, 1, 2 !\n");
        exit(1);
    }


    bin_size = wg_2d_interval_vector[0][0].end_pos - wg_2d_interval_vector[0][0].start_pos;
    start_idx = start_pos / bin_size;
    end_idx = end_pos / bin_size;

    total_depth = 0;
    for (int idx = start_idx; idx < end_idx; idx++)
    {
        total_depth += wg_2d_interval_vector[tid][idx].depth[hap_type];
        total_high_qual_depth += wg_2d_interval_vector[tid][idx].high_qual_depth[hap_type];
    }

    out_depth = total_depth / (end_pos - start_pos);
    out_high_qual_depth = total_high_qual_depth / (end_pos - start_pos);

    if (debug)
    {
        fprintf(stderr, "%d\t%d\t%d\thap_type=%d\tlength = %d\t%.2f\n", tid, start_pos, end_pos, hap_type, end_pos-start_pos, out_depth);
    }

    return 0;
}

int round2_merge_candidate_regions(const Settings & global_settings, const std::vector <std::vector <CNVInterval> > wg_2d_interval_vector, std::vector <CNVCall> & wg_r1_merged_cnv_region_vector, std::vector <CNVCall> & wg_r2_merged_cnv_region_vector, bool recursive)
{
    CNVCall temp_cnv_call; 
    CNVCall * ptr_next_cnvcall, * ptr_cnvcall; 

    int merged_start_pos;
    int merged_end_pos;
    int gap_start_pos;
    int gap_end_pos;
    float merged_avg_depth;
    float merged_high_qual_depth;
    float gap_avg_depth;
    float gap_high_qual_depth;
    static int num_loops = 0;
    int gap_length;
 
    num_loops++;
    fprintf(stderr, "merging candidatge cnv regions (round 2, loop %d)\n", num_loops);
    wg_r2_merged_cnv_region_vector.clear();

    if (wg_r1_merged_cnv_region_vector.size() == 0) { return 0; }

    temp_cnv_call = wg_r1_merged_cnv_region_vector[0];
    for (int idx = 1; idx < wg_r1_merged_cnv_region_vector.size(); idx++)
    {
        ptr_next_cnvcall = & wg_r1_merged_cnv_region_vector[idx];
        merged_start_pos = temp_cnv_call.start_pos;
        merged_end_pos = ptr_next_cnvcall->end_pos;
        gap_start_pos = temp_cnv_call.end_pos;
        gap_end_pos = ptr_next_cnvcall->start_pos;
        gap_length = gap_end_pos - gap_start_pos;

        if (debug && idx % 100 == 0)
        {
            fprintf(stderr, "processed %d CNV regions\n", idx); 
        }

        if (ptr_next_cnvcall->tid != temp_cnv_call.tid || ptr_next_cnvcall->cnv_type != temp_cnv_call.cnv_type || ptr_next_cnvcall->cnv_hap_type != temp_cnv_call.cnv_hap_type || gap_length * 2 > (merged_end_pos - merged_start_pos))  
        { 
            wg_r2_merged_cnv_region_vector.push_back(temp_cnv_call);  
            temp_cnv_call = *ptr_next_cnvcall; 
            continue;
        }

        calculate_interval_avg_depth(wg_2d_interval_vector, temp_cnv_call.tid, gap_start_pos, gap_end_pos, temp_cnv_call.cnv_hap_type, gap_avg_depth, gap_high_qual_depth);

        merged_avg_depth = temp_cnv_call.depth[temp_cnv_call.cnv_hap_type] * temp_cnv_call.length() + gap_avg_depth * gap_length + ptr_next_cnvcall->depth[temp_cnv_call.cnv_hap_type] * ptr_next_cnvcall->length();
        merged_avg_depth /= (merged_end_pos - merged_start_pos);

        if (temp_cnv_call.cnv_type == -1 && merged_avg_depth < global_settings.lower_depth_cutoff[temp_cnv_call.cnv_hap_type] ) 
        {
            temp_cnv_call.end_pos = merged_end_pos; 
            for (int hap_type = 0; hap_type < 3; hap_type++)
            {
                calculate_interval_avg_depth(wg_2d_interval_vector, temp_cnv_call.tid, gap_start_pos, gap_end_pos, hap_type, gap_avg_depth, gap_high_qual_depth);

                temp_cnv_call.depth[hap_type] = temp_cnv_call.depth[temp_cnv_call.cnv_hap_type] * temp_cnv_call.length() + gap_avg_depth * gap_length + ptr_next_cnvcall->depth[temp_cnv_call.cnv_hap_type] * ptr_next_cnvcall->length();
                temp_cnv_call.depth[hap_type] /= (merged_end_pos - merged_start_pos);

                temp_cnv_call.high_qual_depth[hap_type] = temp_cnv_call.high_qual_depth[temp_cnv_call.cnv_hap_type] * temp_cnv_call.length() + gap_high_qual_depth * gap_length + ptr_next_cnvcall->high_qual_depth[temp_cnv_call.cnv_hap_type] * ptr_next_cnvcall->length();
                temp_cnv_call.high_qual_depth[hap_type] /= (merged_end_pos - merged_start_pos);
            }
            continue;
        }else{
            
            wg_r2_merged_cnv_region_vector.push_back(temp_cnv_call);  
            temp_cnv_call = *ptr_next_cnvcall; 
            continue;
        }
    }
    wg_r2_merged_cnv_region_vector.push_back(temp_cnv_call);

    

    if (recursive && wg_r1_merged_cnv_region_vector.size() != wg_r2_merged_cnv_region_vector.size() )
    {
        wg_r1_merged_cnv_region_vector.clear();
        wg_r1_merged_cnv_region_vector = wg_r2_merged_cnv_region_vector;
        wg_r2_merged_cnv_region_vector.clear();
        round2_merge_candidate_regions(global_settings, wg_2d_interval_vector, wg_r1_merged_cnv_region_vector, wg_r2_merged_cnv_region_vector, recursive);
    }else{
        if (debug)
        {
            fprintf(stderr, "number of loops in round 2 = %d\n", num_loops);
            std::string out_file = "wg_r2_merged_cnv_region_vector.txt";
            FILE * out_fp;
            out_fp = fopen(out_file.c_str(), "w");
            if (out_fp == NULL)
            {
                fprintf(stderr, "ERROR! Failed to open file: %s\n", out_file.c_str());
            }
            for (auto & cnv_call : wg_r1_merged_cnv_region_vector)
            {
                cnv_call.output_line(out_fp, global_settings.chr_info);
                fprintf(out_fp, "\n");
            }
            fclose(out_fp);
        }
    }

    return 0;
}


int round1_merge_candidate_regions(const Settings & global_settings, std::vector <CNVCall> & wg_candidate_cnv_region_vector, std::vector <CNVCall> & wg_r1_merged_cnv_region_vector)
{
    // naive merge 
    
    CNVCall temp_cnv_call; 
    CNVCall * ptr_next_cnvcall, * ptr_cnvcall; 
    std::vector <CNVCall> each_haptype_cnv_region_vector;

    int merged_start_pos;
    int merged_end_pos; 
    int cnv_type = -1;
    int hap_type_flag; 

    wg_r1_merged_cnv_region_vector.clear();

    if (wg_candidate_cnv_region_vector.size() == 0) { return 0; }


    for (int hap_type = 1; hap_type < 3; hap_type++)
    {

        hap_type_flag = 1<< hap_type;

        each_haptype_cnv_region_vector.clear();
        each_haptype_cnv_region_vector.shrink_to_fit();
        for (int idx= 0; idx < wg_candidate_cnv_region_vector.size(); idx++)
        {
            ptr_cnvcall = & wg_candidate_cnv_region_vector[idx];
            if (ptr_cnvcall->flag & hap_type_flag)
            {
                ptr_cnvcall->cnv_hap_type = hap_type;
                ptr_cnvcall->cnv_type = cnv_type;
                each_haptype_cnv_region_vector.push_back(wg_candidate_cnv_region_vector[idx]);
            }
        }
        each_haptype_cnv_region_vector.shrink_to_fit();
        temp_cnv_call = each_haptype_cnv_region_vector[0];
        for (int idx = 1; idx < each_haptype_cnv_region_vector.size(); idx++)
        {
            ptr_next_cnvcall = & each_haptype_cnv_region_vector[idx];
            if (ptr_next_cnvcall->tid == temp_cnv_call.tid && ptr_next_cnvcall->start_pos <= temp_cnv_call.end_pos && ptr_next_cnvcall->cnv_type == temp_cnv_call.cnv_type && ptr_next_cnvcall->cnv_hap_type == temp_cnv_call.cnv_hap_type)
            {
                merged_start_pos = temp_cnv_call.start_pos;
                merged_end_pos = ptr_next_cnvcall->end_pos;
                
                for (int hap_type = 0; hap_type < 3; hap_type++)
                {
                    temp_cnv_call.depth[hap_type] = temp_cnv_call.depth[hap_type] * temp_cnv_call.length() + ptr_next_cnvcall->depth[hap_type] * ptr_next_cnvcall->length();
                    temp_cnv_call.depth[hap_type] /= merged_end_pos - merged_start_pos;

                    temp_cnv_call.high_qual_depth[hap_type] = temp_cnv_call.high_qual_depth[hap_type] * temp_cnv_call.length() + ptr_next_cnvcall->high_qual_depth[hap_type] * ptr_next_cnvcall->length();
                    temp_cnv_call.high_qual_depth[hap_type] /= merged_end_pos - merged_start_pos;
                }
                temp_cnv_call.end_pos = merged_end_pos;
                continue;
            }else{
                wg_r1_merged_cnv_region_vector.push_back(temp_cnv_call);  
                temp_cnv_call = *ptr_next_cnvcall; 
                continue;
            }
        }
        wg_r1_merged_cnv_region_vector.push_back(temp_cnv_call);
    }

    
    if (debug)
    {
        std::string out_file = "wg_r1_merged_cnv_region_vector.txt";
        FILE * out_fp;
        out_fp = fopen(out_file.c_str(), "w");
        if (out_fp == NULL)
        {
            fprintf(stderr, "ERROR! Failed to open file: %s\n", out_file.c_str());
        }
        for (auto & cnv_call : wg_r1_merged_cnv_region_vector)
        {
            cnv_call.output_line(out_fp, global_settings.chr_info);
            fprintf(out_fp, "\n");
        }
        fclose(out_fp);
    }
    
    return 0; 
}


int get_depth_cutoff_for1hap(int hap_type, const std::vector <QuantileNumbers> & depth_quantiles_vector, const std::vector <QuantileNumbers> & hiqh_qual_depth_quantiles_vector, float & upper_depth_cutoff, float & lower_depth_cutoff, float &upper_depth_pvalue, float & lower_depth_pvalue)
{
    float mean_depth;

    mean_depth = depth_quantiles_vector[hap_type].q[500];
    lower_depth_cutoff = mean_depth / 10.0;
    upper_depth_cutoff = mean_depth * 1.5;

    upper_depth_pvalue = 1.0;
    lower_depth_pvalue = 1.0;
    for (size_t i = 0; i < depth_quantiles_vector[hap_type].q.size(); i++)
    {
        if (depth_quantiles_vector[hap_type].q[i] > lower_depth_cutoff)
        {
            lower_depth_pvalue = (float) depth_quantiles_vector[hap_type].b[i];
            break;
        }
    }

    for (size_t i = 0; i < depth_quantiles_vector[hap_type].q.size(); i++)
    {
        if (depth_quantiles_vector[hap_type].q[i] > upper_depth_cutoff)
        {
            upper_depth_pvalue = (float)1.0 - (float)depth_quantiles_vector[hap_type].b[i];
            break;
        }
    }

    return 0;
}

int calculate_depth_distributions(const std::vector <std::vector <CNVInterval> > & wg_2d_interval_vector, std::vector <QuantileNumbers> & out_depth_quantiles, std::vector <QuantileNumbers> & out_high_qual_depth_quantiles)
{
    int depth;
   
    std::vector <int> depth_count_vector (MAX_DEPTH+1, 0);
    std::vector <int> high_qual_depth_count_vector (MAX_DEPTH+1, 0); 

    std::vector <std::vector <int> > depth_count_2d_vector;
    std::vector <std::vector <int> > hiqh_qual_depth_count_2d_vector;


    for (int hap_type = 0; hap_type < 3; hap_type ++)
    {
        std::vector <int> v;
        depth_count_2d_vector.push_back(v);
        hiqh_qual_depth_count_2d_vector.push_back(v);

        for (int i = 0; i < MAX_DEPTH+1; i++)
        {
            depth_count_2d_vector[hap_type].push_back(0);
            hiqh_qual_depth_count_2d_vector[hap_type].push_back(0);
        }

        for (int tid = 0; tid < wg_2d_interval_vector.size(); tid++)
        {
            for (int idx = 0; idx < wg_2d_interval_vector[tid].size(); idx ++)
            {
                depth = wg_2d_interval_vector[tid][idx].depth[hap_type] * 100;
                if (wg_2d_interval_vector[tid][idx].depth[0] + wg_2d_interval_vector[tid][idx].depth[1] + wg_2d_interval_vector[tid][idx].depth[2]< 0.001) //skip regions where total depth is 0
                {
                    continue;
                }

                if (depth > MAX_DEPTH) { depth = MAX_DEPTH;}
                if (depth >= 0)
                {
                    depth_count_2d_vector[hap_type][depth]++;
                }
                
                depth = wg_2d_interval_vector[tid][idx].high_qual_depth[hap_type] * 100;
                if (depth > MAX_DEPTH) { depth = MAX_DEPTH;}
                if (depth >= 0)
                {
                    hiqh_qual_depth_count_2d_vector[hap_type][depth]++;
                }
            }
        }

        calculate_distribution_from_count_vector(depth_count_2d_vector[hap_type], out_depth_quantiles[hap_type]);
        calculate_distribution_from_count_vector(hiqh_qual_depth_count_2d_vector[hap_type], out_high_qual_depth_quantiles[hap_type]);

        for (int i = 0; i < out_depth_quantiles[hap_type].q.size(); i++)
        {
            out_depth_quantiles[hap_type].q[i] /= 100.0;
            out_high_qual_depth_quantiles[hap_type].q[i] /= 100.0;
        }

        fprintf(stderr, "\n\n\nhap_type =%d\n", hap_type);
        fprintf(stderr, "total depth quantiles:\n");
        out_depth_quantiles[hap_type].print(stderr, 0);
        fprintf(stderr, "\nhigh qual depth quantiles:\n");
        out_high_qual_depth_quantiles[hap_type].print(stderr, 0);
        fprintf(stderr, "\n\n");

    }
    return 0;   
}

bool has_overlap_with_cr_interval_tree(const CNVInterval cnv_interval, cgranges_t * cr_interval_tree, CHR_INFO * chr_info)
{
    int ret;
    std::vector <size_t> output_index_vector;

    ret = search_overlap_from_cr_interval_tree(cr_interval_tree, cnv_interval.to_interval(chr_info), output_index_vector);

    if (ret > 0)
    {
        return true;
    }else{
        return false;
    }
}

int get_candidate_regions(const std::vector <std::vector <CNVInterval> > & wg_2d_interval_vector, std::vector <CNVCall> & wg_candidate_cnv_region_vector, Settings & global_settings)
{ 
    
    CNVCall candidate_cnv; 
    const CNVInterval * ptr_itv; 
    uint32_t flag;

    for (int hap_type = 1; hap_type < 3; hap_type++)
    {
        get_depth_cutoff_for1hap(hap_type, global_settings.depth_quantiles_vector, global_settings.hiqh_qual_depth_quantiles_vector, global_settings.upper_depth_cutoff[hap_type], global_settings.lower_depth_cutoff[hap_type], global_settings.upper_depth_pvalue[hap_type], global_settings.lower_depth_pvalue[hap_type]);
        fprintf(stderr, "lower_depth_cutoff[%d]=%f\n", hap_type, global_settings.lower_depth_cutoff[hap_type] );
        fprintf(stderr, "lower_depth_pvalue[%d]=%f\n\n", hap_type, global_settings.lower_depth_pvalue[hap_type] );

        fprintf(stderr, "upper_depth_cutoff[%d]=%f\n", hap_type, global_settings.upper_depth_cutoff[hap_type] );
        fprintf(stderr, "upper_depth_pvalue[%d]=%f\n\n\n\n", hap_type, global_settings.upper_depth_pvalue[hap_type] );
    }

    for (int tid = 0; tid < wg_2d_interval_vector.size(); tid++)
    {
        candidate_cnv.tid = tid;
        for (int idx = 0; idx < wg_2d_interval_vector[tid].size(); idx++)
        {
            flag = 0; 
            ptr_itv = &(wg_2d_interval_vector[tid][idx]); 
            if (has_overlap_with_cr_interval_tree(wg_2d_interval_vector[tid][idx], global_settings.black_list_cr_interval_tree, global_settings.chr_info)) { continue; }
            for (int hap_type = 1; hap_type < 3; hap_type ++)
            {
                if (ptr_itv->depth[0] > (ptr_itv->depth[0] + ptr_itv->depth[1] + ptr_itv->depth[2]) * 0.33) { continue; }
                if (ptr_itv->depth[hap_type] < global_settings.lower_depth_cutoff[hap_type])
                {
                    flag |= (1<<hap_type); 
                }
            }
            if (flag == 0){ continue; }
            
            candidate_cnv.start_pos = ptr_itv->start_pos; 
            candidate_cnv.end_pos = ptr_itv->end_pos;
            candidate_cnv.flag = flag;

            for (int hap_type = 0; hap_type < 3; hap_type ++)
            {
                candidate_cnv.depth[hap_type] = ptr_itv->depth[hap_type];
                candidate_cnv.high_qual_depth[hap_type] = ptr_itv->high_qual_depth[hap_type];
            }
            wg_candidate_cnv_region_vector.push_back(candidate_cnv);
        }
    }

    fprintf(stderr, "number of candidate regions: %lu\n", wg_candidate_cnv_region_vector.size());

    return 0; 
}



int generate_whole_genome_interval_vector(std::vector <std::vector <CNVInterval>> & wg_2d_interval_vector, CHR_INFO * chr_info, int bin_size)
{
    int n_chr;
    int tid;
    char * chrname;
    int chr_length;
    
    CNVInterval itv;
    
    n_chr = chr_info->chr_length_list->size;
    fprintf(stderr, "generating whole genome intervals\n");

    for (tid = 0; tid < n_chr; tid++)
    {
        std::vector <CNVInterval> one_chr_interval_vector;
        wg_2d_interval_vector.push_back(one_chr_interval_vector);
    }

    for (tid = 0; tid < n_chr; tid++)
    {
        chr_length = chr_info->chr_length_list->data_list[tid]; 
        chrname = chr_info->chrname_list->data_list[tid];
        itv.tid = tid; 
        for (int start_pos = 0; start_pos < chr_length + bin_size; start_pos += bin_size)
        {
            itv.start_pos = start_pos;
            itv.end_pos = start_pos + bin_size; 
            wg_2d_interval_vector[tid].push_back(itv);
        }
    }

    for (tid = 0; tid < n_chr; tid++)
    {
        if(debug)
        {
            fprintf(stderr, "before fit: wg_2d_interval_vector[%d].capacity = %d\n", tid, wg_2d_interval_vector[tid].capacity());
        }
        
        wg_2d_interval_vector[tid].shrink_to_fit();

        if (debug)
        {
            fprintf(stderr, "after fit: wg_2d_interval_vector[%d].capacity = %d\n", tid, wg_2d_interval_vector[tid].capacity());
        }
        
    }
    
    return 0;

}


int output_cnv_candidate_regions(const std::vector <CNVCall> & wg_candidate_cnv_region_vector, std::string out_file, const CHR_INFO * chr_info)
{
    FILE * out_fp;
    char * chrname;
    int cnv_type; // 1 for duplication, -1 for deletion

    CNVCall temp_cnvcall;
    
    fprintf(stderr, "writing to the output file: %s\n", out_file.c_str());

    fprintf(stderr, "wg_candidate_cnv_region_vector.size()=%d\n", wg_candidate_cnv_region_vector.size());

    out_fp = fopen(out_file.c_str(), "w");

    for (int idx = 0; idx < wg_candidate_cnv_region_vector.size(); idx++)
    {
        temp_cnvcall = wg_candidate_cnv_region_vector[idx];
        temp_cnvcall.cnv_hap_type = 0;
        temp_cnvcall.cnv_type = -1;
        if (temp_cnvcall.flag & HAP1_DEL)
        {
            temp_cnvcall.cnv_hap_type = 1;
        }
        if (temp_cnvcall.flag & HAP2_DEL)
        {
            temp_cnvcall.cnv_hap_type = 2; 
        }
        if (temp_cnvcall.flag & HAP1_DEL && temp_cnvcall.flag & HAP2_DEL)
        {
            temp_cnvcall.cnv_hap_type = 3;
        }
        temp_cnvcall.output_line(out_fp, chr_info);  
        fprintf(out_fp, "\n");
    }

    fclose(out_fp);


    return 0;
}


int calculate_depth_for_each_interval(std::string hap_type_read_depth_file, std::vector <std::vector <CNVInterval> > & wg_2d_interval_vector, const CHR_INFO * chr_info, int bin_size)
{
    gzFile hap_type_read_depth_fp;
    char * line;
    int tid, start_pos, end_pos, start_idx; 
    float all_hap_high_qual_read_depth, all_hap_total_read_depth;
    float hap1_high_qual_read_depth, hap1_total_read_depth;
    float hap2_high_qual_read_depth, hap2_total_read_depth;
    float hap0_high_qual_read_depth, hap0_total_read_depth; 
    CNVInterval * ptr_itv; 
    int l;

    tid = start_pos = end_pos = 0;
    all_hap_high_qual_read_depth = all_hap_total_read_depth = hap1_high_qual_read_depth = hap1_total_read_depth = hap2_high_qual_read_depth = hap2_total_read_depth = hap0_high_qual_read_depth = hap0_total_read_depth = 0.0;

    line = new char[LOCAL_LINE_MAX];
    hap_type_read_depth_fp = gzopen(hap_type_read_depth_file.c_str(), "r");
    if (Z_NULL == hap_type_read_depth_fp) 
    {
        fprintf(stderr, "ERROR! Failed to open file for reading: %s\n", hap_type_read_depth_file.c_str());
        exit(1);
    }

    fprintf(stderr, "processing hap type read depth file: %s\n", hap_type_read_depth_file.c_str()); 


    int i;
    i = 0;

    while (gzgets(hap_type_read_depth_fp, line, LOCAL_LINE_MAX))
    {
        if (line[0] == '#') { continue; }
        sscanf(line, "%d\t%d\t%d%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", &tid, &start_pos, &end_pos, &all_hap_high_qual_read_depth, &all_hap_total_read_depth, &hap1_high_qual_read_depth, &hap1_total_read_depth, &hap2_high_qual_read_depth, &hap2_total_read_depth, &hap0_high_qual_read_depth, &hap0_total_read_depth);

        start_idx = start_pos / bin_size;
        l = end_pos - start_pos;
        ptr_itv = &wg_2d_interval_vector[tid][start_idx];

        ptr_itv->depth[0] += hap0_total_read_depth * (double)l / (double)bin_size;
        ptr_itv->depth[1] += hap1_total_read_depth * (double)l / (double)bin_size;
        ptr_itv->depth[2] += hap2_total_read_depth * (double)l / (double)bin_size;
       

        ptr_itv->high_qual_depth[0] += hap0_high_qual_read_depth  * (double)l / (double)bin_size;
        ptr_itv->high_qual_depth[1] += hap1_high_qual_read_depth  * (double)l / (double)bin_size;
        ptr_itv->high_qual_depth[2] += hap2_high_qual_read_depth  * (double)l / (double)bin_size;
      
            
        i++;
        if (i % 1000000 == 0)
        {
            std::cerr << "processed " << i << " lines in read depth file." << std::endl;
        }
    }

    if (debug)
    {
        std::string out_file = "wg_2d_interval_vector.txt";
        FILE * out_fp;
        out_fp = fopen(out_file.c_str(), "w");
        if (out_fp == NULL)
        {
            fprintf(stderr, "ERROR! failed to open file: %s\n", out_file.c_str());
            exit(1);
        }
        for (int tid = 0; tid < wg_2d_interval_vector.size(); tid++)
        {
            for (auto itv : wg_2d_interval_vector[tid])
            {
                itv.output_line(out_fp, chr_info);
                fprintf(out_fp, "\n");
            }
        }
    }
    gzclose(hap_type_read_depth_fp);

    delete [] line;
    return 0; 
}

int read_blacklist_file(const std::string & blacklist_bed_file, Settings & global_settings)
{
    std::vector <Interval> blacklist_interval_vector;

    get_interval_vector_from_bed_file(blacklist_bed_file, global_settings.chr_info, blacklist_interval_vector);

    fprintf(stderr, "number of intervals in blacklist file = %d\n", blacklist_interval_vector.size());
    global_settings.black_list_cr_interval_tree = generate_cr_interval_tree(blacklist_interval_vector);

    return 0;
}



int call_small_deletions(const std::string hap_type_read_depth_file, const std::string weird_reads_file, const std::string bcd22_file, const std::string faidx_file, const std::string blacklist_bed_file, const std::string out_file)
{

    Settings global_settings;
    std::vector <std::vector <CNVInterval> > wg_2d_interval_vector; 
    std::vector <CNVCall> wg_candidate_cnv_region_vector;
    std::vector <CNVCall> wg_r1_merged_cnv_region_vector;
    std::vector <CNVCall> wg_r2_merged_cnv_region_vector;


    global_settings.chr_info = get_chr_info(faidx_file.c_str());
    global_settings.bin_size = 100; 

    read_blacklist_file(blacklist_bed_file, global_settings); 

    generate_whole_genome_interval_vector(wg_2d_interval_vector, global_settings.chr_info, global_settings.bin_size);

    calculate_depth_for_each_interval(hap_type_read_depth_file, wg_2d_interval_vector, global_settings.chr_info, global_settings.bin_size);
    
    calculate_depth_distributions(wg_2d_interval_vector, global_settings.depth_quantiles_vector, global_settings.hiqh_qual_depth_quantiles_vector);
    
    get_candidate_regions(wg_2d_interval_vector, wg_candidate_cnv_region_vector, global_settings);

    output_cnv_candidate_regions(wg_candidate_cnv_region_vector, "candidate_cnv_regions.txt", global_settings.chr_info);

    fprintf(stderr, "merging candidatge cnv regions (round 1)\n");
    round1_merge_candidate_regions(global_settings, wg_candidate_cnv_region_vector, wg_r1_merged_cnv_region_vector);

    round2_merge_candidate_regions(global_settings, wg_2d_interval_vector, wg_r1_merged_cnv_region_vector, wg_r2_merged_cnv_region_vector, true);

    return 0;
}

int main (int argc, char * argv[])
{
    std::string usage; 
    std::string hap_type_read_depth_file, weird_reads_file, bcd22_file, faidx_file, blacklist_bed_file, out_file;

    usage = "Usage: call_small_deletions <in.hap_type_read_depth_file> <in_weird_reads_file> <in.bcd22> <faidx_file> <blacklist.bed> <out_file>\n";

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
    weird_reads_file         = argv[2];
    bcd22_file               = argv[3];
    faidx_file               = argv[4];
    blacklist_bed_file       = argv[5];
    out_file                 = argv[6];

    call_small_deletions(hap_type_read_depth_file, weird_reads_file, bcd22_file, faidx_file, blacklist_bed_file, out_file);

    return 0;
}