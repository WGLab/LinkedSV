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

int generate_whole_genome_interval_vector(std::vector <std::vector <CNVInterval>> & wg_2d_interval_vector, const CHR_INFO * chr_info, int bin_size)
{
    int n_chr;
    int tid;
    char * chrname;
    int chr_length;
    int n_bin;
    
    CNVInterval itv;
    
    n_chr = chr_info->chr_length_list->size;
    if (error_level > 1)
    {
        fprintf(stderr, "generating whole genome intervals\n");
    }
    

    std::vector <CNVInterval> one_chr_interval_vector;
    for (tid = 0; tid < n_chr; tid++)
    {
        wg_2d_interval_vector.push_back(one_chr_interval_vector);
    }

    for (tid = 0; tid < n_chr; tid++)
    {
        chr_length = chr_info->chr_length_list->data_list[tid]; 
        chrname = chr_info->chrname_list->data_list[tid];
        itv.tid = tid; 
        n_bin = chr_length / bin_size + 10;
        wg_2d_interval_vector[tid].reserve(n_bin);
        for (int start_pos = 0; start_pos < chr_length + bin_size; start_pos += bin_size)
        {
            itv.start_pos = start_pos;
            itv.end_pos = start_pos + bin_size; 
            wg_2d_interval_vector[tid].push_back(itv);
        }
    }
    
    return 0;

}


int calculate_depth_distributions(const std::vector <std::vector <CNVInterval> > & wg_2d_interval_vector, std::vector <QuantileNumbers> & out_depth_quantiles, std::vector <QuantileNumbers> & out_high_qual_depth_quantiles)
{
    int depth;

    std::vector <std::vector <int> > depth_count_2d_vector;
    std::vector <std::vector <int> > hiqh_qual_depth_count_2d_vector;


    for (int hap_type = 0; hap_type < 4; hap_type ++)
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
                if (wg_2d_interval_vector[tid][idx].depth[3] < 0.001) //skip regions where total depth is 0
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

        if (error_level > 1)
        {
            fprintf(stderr, "\n\n\nhap_type =%d\n", hap_type);
            fprintf(stderr, "total depth quantiles:\n");
            out_depth_quantiles[hap_type].print(stderr, 0);
            fprintf(stderr, "\nhigh qual depth quantiles:\n");
            out_high_qual_depth_quantiles[hap_type].print(stderr, 0);
            fprintf(stderr, "\n\n");
        }
    }
    
    return 0;   
}



int get_depth_for_each_interval(std::string hap_type_read_depth_file, std::vector <std::vector <CNVInterval> > & wg_2d_interval_vector, const CHR_INFO * chr_info, int bin_size)
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

    if (error_level > 1)
    {
        fprintf(stderr, "processing hap type read depth file: %s\n", hap_type_read_depth_file.c_str()); 
    }
    
    int i;
    i = 0;

    while (gzgets(hap_type_read_depth_fp, line, LOCAL_LINE_MAX))
    {
        if (line[0] == '#') { continue; }
        sscanf(line, "%d\t%d\t%d%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", &tid, &start_pos, &end_pos, &all_hap_high_qual_read_depth, &all_hap_total_read_depth, &hap1_high_qual_read_depth, &hap1_total_read_depth, &hap2_high_qual_read_depth, &hap2_total_read_depth, &hap0_high_qual_read_depth, &hap0_total_read_depth);

        start_idx = start_pos / bin_size;
        l = end_pos - start_pos;
        ptr_itv = &wg_2d_interval_vector[tid][start_idx];

        ptr_itv->depth[0] += hap0_total_read_depth * (float)l / (float)bin_size;
        ptr_itv->depth[1] += hap1_total_read_depth * (float)l / (float)bin_size;
        ptr_itv->depth[2] += hap2_total_read_depth * (float)l / (float)bin_size;
        ptr_itv->depth[3] += all_hap_total_read_depth * (float)l / (float)bin_size;

        ptr_itv->high_qual_depth[0] += hap0_high_qual_read_depth  * (float)l / (float)bin_size;
        ptr_itv->high_qual_depth[1] += hap1_high_qual_read_depth  * (float)l / (float)bin_size;
        ptr_itv->high_qual_depth[2] += hap2_high_qual_read_depth  * (float)l / (float)bin_size;
        ptr_itv->high_qual_depth[3] += all_hap_high_qual_read_depth  * (float)l / (float)bin_size;
            
        i++;
        if (error_level > 1 && i % 10000000 == 0)
        {
            std::cerr << "read " << i << " lines in read depth file." << std::endl;
        }
    }

    gzclose(hap_type_read_depth_fp);

    delete [] line;
    return 0; 
}