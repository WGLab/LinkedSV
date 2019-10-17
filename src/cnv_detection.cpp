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

int read_gap_region_file(const std::string & gap_region_bed_file, CNVSettings & global_settings)
{
    get_interval_vector_from_bed_file(gap_region_bed_file, global_settings.chr_info, global_settings.gap_region_interval_vector);

    if (error_level > 1)
    {
        fprintf(stderr, "number of intervals in gap_region file = %d\n", global_settings.gap_region_interval_vector.size());
    }
    return 0;
}


double calculate_distance_to_next (double wg_sd_depth, const std::vector <std::vector <double>> & wg_2d_depth_vector, std::vector <MergedCNVInterval> & merged_cnv_vector_of_1chr, int32_t idx)
{
    double zscore;
    double n1, n2;
    int32_t next_idx;
    double mean1, mean2;

    next_idx = merged_cnv_vector_of_1chr[idx].next_idx;

    n1 = merged_cnv_vector_of_1chr[idx].end_idx - merged_cnv_vector_of_1chr[idx].start_idx;
    n2 = merged_cnv_vector_of_1chr[next_idx].end_idx - merged_cnv_vector_of_1chr[next_idx].start_idx;

    mean1 = merged_cnv_vector_of_1chr[idx].mean_depth;
    mean2 = merged_cnv_vector_of_1chr[next_idx].mean_depth;
    
    //double r1, r2;
    //r1 = var1 / n1;
    //r2 = var2 / n2;

  
    zscore = (mean1 - mean2) / wg_sd_depth * sqrt (n1 * n2) / sqrt(n1 + n2);
   
    if (zscore < 0)
    {
        zscore = -zscore;
    }
    return zscore;
}


int merge_two_intervals(double wg_sd_depth, const std::vector <std::vector <double>> & wg_2d_depth_vector, std::vector <MergedCNVInterval> & merged_cnv_vector_of_1chr, int32_t idx)
{
    if (idx < 0 || idx >= merged_cnv_vector_of_1chr.size()-1)
    {
        return 0;
    }

    int32_t next_idx;

    next_idx = merged_cnv_vector_of_1chr[idx].next_idx;
  
    if (next_idx < 0 || next_idx >= merged_cnv_vector_of_1chr.size())
    {
        return 0;
    }

    if (merged_cnv_vector_of_1chr[next_idx].is_active == false)
    {
        fprintf(stderr, "ERROR! next_idx is not active!\n");
        exit(1);
    }
    if (merged_cnv_vector_of_1chr[idx].is_active == false)
    {
        fprintf(stderr, "ERROR! idx is not active!\n");
        exit(1);
    }

    // begin merge

    merged_cnv_vector_of_1chr[next_idx].is_active = false;

    merged_cnv_vector_of_1chr[idx].end_idx = merged_cnv_vector_of_1chr[next_idx].end_idx;
    merged_cnv_vector_of_1chr[idx].next_idx = merged_cnv_vector_of_1chr[next_idx].next_idx;
    merged_cnv_vector_of_1chr[merged_cnv_vector_of_1chr[idx].next_idx].prev_idx = idx;
    merged_cnv_vector_of_1chr[idx].calculate_mean_depth(wg_2d_depth_vector);
    merged_cnv_vector_of_1chr[idx].calculate_variance(wg_2d_depth_vector);
 
    if (merged_cnv_vector_of_1chr[idx].next_idx >= merged_cnv_vector_of_1chr.size())
    {
        merged_cnv_vector_of_1chr[idx].distance_to_next = 1e30;
    }else{
        merged_cnv_vector_of_1chr[idx].distance_to_next = calculate_distance_to_next(wg_sd_depth, wg_2d_depth_vector, merged_cnv_vector_of_1chr, idx);
    }


    return 0;
    
}


int detect_large_cnv_from1chr (CNVSettings & global_settings, const std::vector <std::vector <double>> & wg_2d_depth_vector, size_t tid, std::vector <MergedCNVInterval> & merged_cnv_vector_of_1chr)
{
    MergedCNVInterval merged_cnv_interval;
    std::multimap <double, int32_t> delta_map;
    static int64_t n_merged_intervals = 0;
    int chrom_copy_number;

    chrom_copy_number = 2;
  
    if (global_settings.special_chrom_tid_map.count(tid) > 0)
    {
        MergedCNVInterval whole_chr_cnv_interval;
        double whole_chr_nongap_mean_depth;
        if (global_settings.special_chrom_tid_map[tid] == 0 || global_settings.special_chrom_tid_map[tid] == 3)
        {
            return 0;
        }

        whole_chr_cnv_interval.tid = tid;
        whole_chr_cnv_interval.start_idx = 0;
        whole_chr_cnv_interval.end_idx = wg_2d_depth_vector[tid].size();
        whole_chr_cnv_interval.is_active = true;
        whole_chr_cnv_interval.distance_to_next = 1e30;
        whole_chr_cnv_interval.calculate_mean_depth(wg_2d_depth_vector);
        double gap_region_depth = 0;
        int32_t gap_bin_cnt;
        gap_bin_cnt = 0;
        for (Interval & itv: global_settings.gap_region_interval_vector)
        {
            if (itv.tid != tid) { continue; }
            int32_t start_idx = itv.start_pos / global_settings.bin_size; 
            int32_t end_idx = itv.end_pos / global_settings.bin_size;
            for (int32_t idx = start_idx; idx < end_idx; idx++)
            {
                gap_region_depth += wg_2d_depth_vector[tid][idx];
                gap_bin_cnt++;
            }
        }
        whole_chr_nongap_mean_depth = whole_chr_cnv_interval.mean_depth * (whole_chr_cnv_interval.end_idx - whole_chr_cnv_interval.start_idx) - gap_region_depth;
        whole_chr_nongap_mean_depth /= whole_chr_cnv_interval.end_idx - whole_chr_cnv_interval.start_idx - gap_bin_cnt;

        if (whole_chr_nongap_mean_depth >= 0.75 * global_settings.wg_mean_depth)
        {
            chrom_copy_number = 2;
        }
        else if (whole_chr_nongap_mean_depth >= 0.25 * global_settings.wg_mean_depth && whole_chr_nongap_mean_depth < 0.75 * global_settings.wg_mean_depth)
        {
            chrom_copy_number = 1;
        }
        else if (whole_chr_nongap_mean_depth < 0.25 * global_settings.wg_mean_depth)
        {
            chrom_copy_number = 0;
        }

        if (global_settings.special_chrom_tid_map[tid] == 1) // chrY
        {
            if (chrom_copy_number > 1)
            {
                chrom_copy_number = 1;
            }
        }
        else if (global_settings.special_chrom_tid_map[tid] == 2) // chrX
        {
            if (chrom_copy_number < 1){
                chrom_copy_number = 1;
            }
        }
        if (error_level > 1)
        {
            fprintf(stderr, "tid=%d, whole_chr_nongap_mean_depth=%.2f, chrom_copy_number=%d\n", tid, whole_chr_nongap_mean_depth, chrom_copy_number);
        }
    }

    global_settings.chrom_copy_number_vector[tid] = chrom_copy_number;
    
    merged_cnv_vector_of_1chr.clear();
    int32_t n_active_bin = 0;
    // init vector
    for (int32_t idx = 0; idx < wg_2d_depth_vector[tid].size(); idx++)
    {
        merged_cnv_interval.tid = tid;
        merged_cnv_interval.start_idx = idx;
        merged_cnv_interval.end_idx = idx + 1;
		merged_cnv_interval.prev_idx = idx - 1;
        merged_cnv_interval.next_idx = idx + 1;
        merged_cnv_interval.is_active = true;
        merged_cnv_interval.calculate_mean_depth(wg_2d_depth_vector);
        merged_cnv_interval.variance = 0.0;
        merged_cnv_interval.distance_to_next = 1e30;
        
        merged_cnv_vector_of_1chr.push_back(merged_cnv_interval);
    }

    // merge elements with the same depth
    if (error_level > 1)
    {
        fprintf(stderr, "merge adjacent bins with nearly equal depth\n");
    }

    int32_t idx, next_idx, prev_idx;
    idx = 0;
    double depth_change;
    while (idx < merged_cnv_vector_of_1chr.size() && merged_cnv_vector_of_1chr[idx].next_idx < merged_cnv_vector_of_1chr.size())
    {
        if (merged_cnv_vector_of_1chr[idx].is_active == false) { continue; }
        while(1)
        {
            depth_change = merged_cnv_vector_of_1chr[idx].mean_depth - merged_cnv_vector_of_1chr[merged_cnv_vector_of_1chr[idx].next_idx].mean_depth;
            if (depth_change > 0.5 || depth_change < -0.5) { break; }
            next_idx = merged_cnv_vector_of_1chr[idx].next_idx;
            if (next_idx >= merged_cnv_vector_of_1chr.size() || next_idx < 0 ) {break;}
            if (merged_cnv_vector_of_1chr[next_idx].is_active == false)
            {
                fprintf(stderr, "ERROR! merged_cnv_vector_of_1chr[next_idx].is_active == false. idx=%d, next_idx=%d\n", idx, next_idx);
                exit(1);
            }
            merged_cnv_vector_of_1chr[next_idx].is_active = false;
            merged_cnv_vector_of_1chr[idx].end_idx = merged_cnv_vector_of_1chr[next_idx].end_idx;
            merged_cnv_vector_of_1chr[idx].next_idx = merged_cnv_vector_of_1chr[next_idx].next_idx;
            merged_cnv_vector_of_1chr[merged_cnv_vector_of_1chr[idx].next_idx].prev_idx = idx;
            if (depth_change != 0.0)
            {
                merged_cnv_vector_of_1chr[idx].calculate_mean_depth(wg_2d_depth_vector);
            }
            
            n_merged_intervals++;
            if (merged_cnv_vector_of_1chr[idx].next_idx >=  merged_cnv_vector_of_1chr.size()) {break;}
            if (merged_cnv_vector_of_1chr[idx].next_idx < 0) {break;}
        }
        idx = merged_cnv_vector_of_1chr[idx].next_idx;
    }


    if (error_level > 1)
    {
        fprintf(stderr, "tid = %d, merged_cnv_vector_of_1chr.size()=%d, merged %d bins with equal depth\n", tid, merged_cnv_vector_of_1chr.size(), n_merged_intervals);
    }

	std::multimap <double, int32_t>::iterator iter;
    
    // build tree
    
    if (error_level > 1)
    {
        fprintf(stderr, "building tree\n");
    }

    n_active_bin = 0;
    for (int32_t idx = 0; idx < merged_cnv_vector_of_1chr.size() - 1; idx++)
    {
        if (merged_cnv_vector_of_1chr[idx].is_active == false) { continue; }
        merged_cnv_vector_of_1chr[idx].calculate_variance(wg_2d_depth_vector);
        n_active_bin++;
        merged_cnv_vector_of_1chr[idx].distance_to_next = calculate_distance_to_next(global_settings.wg_sd_depth, wg_2d_depth_vector, merged_cnv_vector_of_1chr, idx);
        delta_map.insert({merged_cnv_vector_of_1chr[idx].distance_to_next, idx});
    }

    if (error_level > 1 )
    {
        fprintf(stderr, "delta_map.size()=%d, n_active_bin=%d\n", delta_map.size(), n_active_bin);
        fprintf(stderr, "delta_map.begin()->first=%.6f\n", delta_map.begin()->first);
    }

    auto map_iter = delta_map.begin();
    
    while (1)
    {
        map_iter = delta_map.begin();
        if (map_iter->first > global_settings.alpha) {  break; }
        if (delta_map.size() == 0) { break; }

        idx = map_iter->second;
        next_idx = merged_cnv_vector_of_1chr[idx].next_idx;
        prev_idx = merged_cnv_vector_of_1chr[idx].prev_idx;
       
        if (idx >= merged_cnv_vector_of_1chr.size() - 1 || next_idx >= merged_cnv_vector_of_1chr.size())
        {
            delta_map.erase(map_iter);
            continue; 
        }

        if (idx < 0 || next_idx < 0)
        {
            fprintf(stderr, "ERROR! idx or next_idx < 0. idx=%d, next_idx=%d\n", idx, next_idx);
            delta_map.erase(map_iter);
            continue; 
        }

        if (merged_cnv_vector_of_1chr[idx].is_active == false)
        {
            delta_map.erase(map_iter);
            continue; 
        }

        if (merged_cnv_vector_of_1chr[next_idx].is_active == false)
        {
            fprintf(stderr, "ERROR! next_idx is not active! next_idx=%d\n", next_idx);
            delta_map.erase(map_iter);
        }

        if (merged_cnv_vector_of_1chr[idx].distance_to_next != map_iter->first)  // should be skipped 
        {
            delta_map.erase(map_iter);
            continue; 
        }

        merge_two_intervals(global_settings.wg_sd_depth, wg_2d_depth_vector, merged_cnv_vector_of_1chr, idx);
        delta_map.erase(map_iter);
        delta_map.insert({merged_cnv_vector_of_1chr[idx].distance_to_next, idx});
        if (prev_idx >= 0)
        {
            double new_distance_to_next;
            new_distance_to_next = calculate_distance_to_next(global_settings.wg_sd_depth, wg_2d_depth_vector, merged_cnv_vector_of_1chr, prev_idx);
            if (new_distance_to_next != merged_cnv_vector_of_1chr[prev_idx].distance_to_next)
            {
                merged_cnv_vector_of_1chr[prev_idx].distance_to_next = new_distance_to_next;
                delta_map.insert({merged_cnv_vector_of_1chr[prev_idx].distance_to_next, prev_idx}); 
            }
        }
        
        n_merged_intervals ++;
        if (error_level > 1 && n_merged_intervals % 100000 == 0)
        {
            fprintf(stderr, "total number of merged intervals=%d\n", n_merged_intervals);
            fprintf(stderr, "delta_map.size()=%d\n", delta_map.size());
            fprintf(stderr, "delta_map.begin()->first=%.10f\n", delta_map.begin()->first);
        }
    }
    
    for (size_t idx = 0; idx < merged_cnv_vector_of_1chr.size(); idx++)
    {
        merged_cnv_vector_of_1chr[idx].filter = "LOW_QUAL";
        if (merged_cnv_vector_of_1chr[idx].is_active == false) { continue; }
        merged_cnv_vector_of_1chr[idx].calculate_genotype(wg_2d_depth_vector, global_settings.wg_mean_depth, global_settings.wg_sd_depth, chrom_copy_number);
        if (merged_cnv_vector_of_1chr[idx].cnv_type == 0) { continue; }
        if (merged_cnv_vector_of_1chr[idx].gt_score < global_settings.min_large_cnv_llr) { continue; }
        if (merged_cnv_vector_of_1chr[idx].length(global_settings.bin_size) < global_settings.min_large_cnv_size) { continue; }
        merged_cnv_vector_of_1chr[idx].filter = "PASS";
    }

    return 0;
}


int get_special_chrom_tid_maps(const CHR_INFO * chr_info, std::map<int32_t, int32_t> & special_chrom_tid_map)
{
    int32_t n_chr = chr_info->chrname_list->size;
    char * ctg;
    for (int32_t tid = 0; tid < n_chr; tid++)
    {
        ctg = chr_info->chrname_list->data_list[tid];
        if (strcmp(ctg, "chrX") == 0 || strcmp(ctg, "X") == 0 || strcmp(ctg, "chrx") == 0 || strcmp(ctg, "x") == 0 )
        {
            special_chrom_tid_map[tid]  = 2; // chrX
        }
        if (strcmp(ctg, "chrY") == 0 || strcmp(ctg, "Y") == 0 || strcmp(ctg, "chry") == 0 || strcmp(ctg, "y") == 0)
        {
            special_chrom_tid_map[tid]  = 1; // chrY
        }
        if (strcmp(ctg, "chrM") == 0 || strcmp(ctg, "chrm") == 0 || strcmp(ctg, "MT") == 0 || strcmp(ctg, "M") == 0 )
        {
            special_chrom_tid_map[tid]  = 3; // mitochondria
        }
        if (strcmp(ctg, "hs37d5") == 0 || strcmp(ctg, "hs38d1") == 0 )
        {
            special_chrom_tid_map[tid]  = 0; // blacklist
        }
        if (strstr(ctg, "chrUn_") != NULL || strstr(ctg, "_random") != NULL || strstr(ctg, "chrEBV") != NULL ) // hg19 and hg38
        {
            special_chrom_tid_map[tid]  = 0; // blacklist
        }
        if (strstr(ctg, "GL000") != NULL || strstr(ctg, "NC_007605") != NULL  ) // b37
        {
            special_chrom_tid_map[tid]  = 0; // blacklist
        }
    }

    return 0;
    
}


int cnv_detection(const std::string hap_type_read_depth_file, const std::string faidx_file, const std::string gap_region_bed_file, const std::string out_file, double alpha, double min_LRR, int32_t min_size)
{
    CNVSettings global_settings;
    std::vector <std::vector <CNVInterval> > wg_2d_interval_vector; 

    std::vector <std::vector <double>> wg_2d_depth_vector;
    std::vector <double> v;
    double depth, max_depth;
    FILE * out_fp;


    global_settings.alpha = alpha; // 40
    global_settings.bin_size = 100;
    global_settings.min_large_cnv_size = min_size; // 2e5
    global_settings.min_large_cnv_llr = min_LRR; // 200

    global_settings.chr_info = get_chr_info(faidx_file.c_str());
    
    get_special_chrom_tid_maps(global_settings.chr_info, global_settings.special_chrom_tid_map);

    read_gap_region_file(gap_region_bed_file, global_settings); 

    generate_whole_genome_interval_vector(wg_2d_interval_vector, global_settings.chr_info, global_settings.bin_size);
  
    get_depth_for_each_interval(hap_type_read_depth_file, wg_2d_interval_vector, global_settings.chr_info, global_settings.bin_size);
   
    calculate_depth_distributions (wg_2d_interval_vector, global_settings.depth_quantiles_vector, global_settings.hiqh_qual_depth_quantiles_vector);

    global_settings.wg_mean_depth = global_settings.depth_quantiles_vector[3].q[500];
    global_settings.wg_sd_depth = global_settings.wg_mean_depth;

    if (error_level > 1)
    {
        fprintf(stderr, "whole genome average depth = %.2f\n", global_settings.wg_mean_depth);
    }

    max_depth = global_settings.wg_mean_depth * 5;

    for (size_t tid = 0; tid < wg_2d_interval_vector.size(); tid++)
    {
        wg_2d_depth_vector.push_back(v);
        global_settings.chrom_copy_number_vector.push_back(0);
        for (size_t idx = 0; idx < wg_2d_interval_vector[tid].size(); idx++)
        {
            depth = wg_2d_interval_vector[tid][idx].depth[3];
            if (depth > max_depth){ depth = max_depth; }
            wg_2d_depth_vector[tid].push_back(depth);
        }
    }

    for (Interval & itv: global_settings.gap_region_interval_vector)
    {
        int32_t start_idx = itv.start_pos / global_settings.bin_size;
        int32_t end_idx = itv.end_pos / global_settings.bin_size;
        for (int32_t idx = start_idx; idx < end_idx; idx++)
        {
            if (itv.tid >= wg_2d_depth_vector.size()) { continue; }
            if (idx >= wg_2d_depth_vector[itv.tid].size()) { continue; }
            wg_2d_depth_vector[itv.tid][idx] = global_settings.wg_mean_depth;
        }
    }
    
    std::vector <std::vector <MergedCNVInterval>> wg_merged_cnv_vector;
    for (size_t tid = 0; tid < wg_2d_depth_vector.size(); tid++)
    {
        std::vector <MergedCNVInterval> merged_cnv_vector_of_1chr;
        detect_large_cnv_from1chr(global_settings, wg_2d_depth_vector, tid, merged_cnv_vector_of_1chr);
        wg_merged_cnv_vector.push_back(merged_cnv_vector_of_1chr);
    }

    out_fp = fopen(out_file.c_str(), "a");
    {
        for (size_t tid = 0; tid < wg_merged_cnv_vector.size(); tid++)
        {
            for (size_t idx = 0; idx < wg_merged_cnv_vector[tid].size(); idx++)
            {
                if (wg_merged_cnv_vector[tid][idx].is_active == false) { continue; }
                if (wg_merged_cnv_vector[tid][idx].filter == "PASS") 
                wg_merged_cnv_vector[tid][idx].output_bedpe(global_settings.chr_info, global_settings.bin_size, out_fp);
            }
        }
    }

    fclose(out_fp);


    return 0;    
}


int main (int argc, char * argv[])
{
    std::string usage; 
    std::string hap_type_read_depth_file, faidx_file, gap_region_bed_file, out_file;
    double alpha;
    double min_LRR;
    int32_t min_size;

    usage = "Usage: cnv_detection <in.hap_type_read_depth_file> <faidx_file> <gap_region.bed> <out_file> <alpha> <min_LRR> <min_size>\n";

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
    faidx_file               = argv[2];
    gap_region_bed_file      = argv[3];
    out_file                 = argv[4];
    alpha                    = atof(argv[5]);
    min_LRR                  = atof(argv[6]);
    min_size                 = atoi(argv[7]);

    FILE * fp;
    fp = fopen(out_file.c_str(), "w"); 
    fclose(fp);

    cnv_detection(hap_type_read_depth_file, faidx_file, gap_region_bed_file, out_file, alpha, min_LRR, min_size);

    return 0;
}
