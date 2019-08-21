#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <unordered_map>

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <thread>
#include <algorithm>
#include <limits>

#include <zlib.h>

#include "cluster_reads.h"
#include "tk.h"


const int line_max = 4096;
const int MAX_INNER_SIZE = 100 * 1000; 
const int MAX_GAP_DISTANCE = 100 * 1000;
const int MAX_READ_LENGTH = 500;

int frag_id;
std::unordered_map <std::string, int> weird_read_ids_map; 
std::unordered_map <std::string, int> all_read_ids_map;
int64_t abs_int64_t(int64_t a);

class ReadPairIndex{
public: 
    int r1;
    int r2;
    ReadPairIndex(){
        r1 = -1;
        r2 = -1;
    }
};

int32_t inner_size (Bcd21 & aln1, Bcd21 & aln2)
{
    int isize; // actually not insert size, but the gap length of the two reads
    if (aln1.is_read1() == aln2.is_read1() || aln1.is_read2() == aln2.is_read2())
    {
        return -(1<<30); // aln1 and aln2 are the same read  (not pairs)
    }

    if (aln1.map_direction() == aln2.map_direction())
    {
        return -(1<<30);
    }

    if (aln1.tid != aln2.tid)
    {
        return -(1<<30);
    }

    if (aln1.map_direction() == 1 && aln1.start > aln2.start)
    {
        return -(1<<30);
    }

    if (aln2.map_direction() == 1 && aln2.start > aln1.start)
    {
        return -(1<<30);
    }

    if (aln1.map_direction() == 1)
    {
        isize = aln2.start - aln1.end;
    }

    if (aln2.map_direction() == 1)
    {
        isize = aln1.start - aln2.end;
    }
    return isize;
}



int32_t is_weird_read_pair(Bcd21 & aln1, Bcd21 & aln2, int inner_size_cutoff)
{
    int isize; 
    int32_t ret; 
    isize = 0;
    ret = 0; 
    if (aln1.is_read1() == aln2.is_read1() || aln1.is_read2() == aln2.is_read2())
    {
        return 0; // aln1 and aln2 are the same read  (not pairs)
    }
    if (abs_int64_t(aln1.start - aln2.end) < inner_size_cutoff && abs_int64_t(aln2.start - aln1.end) < inner_size_cutoff)
    {
        return 0; // supporting SV is too small 
    }


    if (aln1.tid != aln2.tid)
    {
        ret |= 0x1; // different chrom
    } else
    {
        if (aln1.is_read1() != aln2.is_read1() && aln1.map_direction() == aln2.map_direction())
        {
            ret |= 0x2; // same mapping direction 
        }

        if (aln1.map_direction() == 1 && aln2.map_direction() == -1 && aln1.start > aln2.start)
        {
            ret |= 0x4; // temdam dup
        }

        if (aln2.map_direction() == 1 && aln1.map_direction() == -1 && aln2.start > aln1.start)
        {
            ret |= 0x4; // temdam dup
        }

        isize = inner_size(aln1, aln2);
        if (isize > inner_size_cutoff){ 
            ret |= 0x8; // large insert size
        }
    }

    return ret; 
}

int convert1line2bcd21(const char * line, Bcd21 & bcd21, char * bcd, char * read_id, char * cigar_string)
{
    sscanf(line, "%d\t%d\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n", &bcd21.tid, &bcd21.start, &bcd21.end, &bcd21.mapq, bcd, &bcd21.hptype, read_id, &bcd21.flag, &bcd21.n_left_clip, &bcd21.n_right_clip, &bcd21.inner_size, &bcd21.mate_tid, &bcd21.mate_pos, cigar_string);
    bcd21.bcd = bcd;
    bcd21.read_id = read_id;
    bcd21.cigar_string = cigar_string;
    bcd21.is_weird_read = 0; 
    return 0;
}

bool sort_bcd21_by_pos (const Bcd21 & a, const Bcd21 & b){ 
    return a.key_start() < b.key_start();
}

bool sort_bcd21_by_read_id (const Bcd21 & a, const Bcd21 & b){ 
    return a.read_id < b.read_id;
}

int64_t abs_int64_t(int64_t a)
{
    if (a > 0){
        return a;
    }else{
        return (int64_t) 0 - a;
    }
}



Bcd22 convert_bcd21vector_to_bcd22 (Settings & global_settings, const std::vector <Bcd21> & fragment_bcd21_vector)
{
    // fragment_bcd21_vector should be sorted by map position. 
    Bcd22 frm;

    frm.tid       = fragment_bcd21_vector[0].tid;
    frm.start     = fragment_bcd21_vector[0].start;  
    frm.end       = fragment_bcd21_vector.back().end;  
    frm.bcd       = fragment_bcd21_vector[0].bcd;
    frm.num_reads = fragment_bcd21_vector.size();
    frm.frag_id   = frag_id;
    frag_id++;

    
    frm.hp0 = frm.hp1 = frm.hp2 = 0;
    frm.map_pos = "";
    frm.map_qual = "";
    frm.num_good_reads = 0;
    for (auto  & bcd21 : fragment_bcd21_vector)
    {
        if (bcd21.hptype == 0) {
            frm.hp0++;
        }else if (bcd21.hptype == 1){
            frm.hp1++;
        }else if (bcd21.hptype == 2){
            frm.hp2++;
        }
        frm.map_pos += std::to_string(bcd21.start) + "," + std::to_string(bcd21.end) + ";";

        if (bcd21.mapq < global_settings.min_mapq)
        {
            frm.map_qual += "0";
        }else{
            frm.map_qual += "1";
            frm.num_good_reads++;
        }
    }
    if (frm.map_pos.size()> 0){
        frm.map_pos.pop_back();
    }

    frm.left_weird_reads_output = frm.right_weird_reads_output = frm.other_weird_reads_output = ""; 
    for (int i = 0; i < fragment_bcd21_vector.size(); i++)
    {
        if (fragment_bcd21_vector[i].is_weird_read != 0)
        {
             
            if (fragment_bcd21_vector[i].start - frm.start < global_settings.inner_size_cutoff)
            {
                frm.n_left_weird_reads += 1; 
                frm.left_weird_reads_output += fragment_bcd21_vector[i].output_read_info(); 
            } else if (frm.end - fragment_bcd21_vector[i].end < global_settings.inner_size_cutoff) {
                frm.n_right_weird_reads += 1;
                frm.right_weird_reads_output += fragment_bcd21_vector[i].output_read_info(); 
            } else {
                frm.other_weird_reads_output += fragment_bcd21_vector[i].output_read_info();
            }
            
        }
    }
    

    if (frm.left_weird_reads_output == "")
    {
        frm.left_weird_reads_output = ".";
    }
    if (frm.right_weird_reads_output == "")
    {
        frm.right_weird_reads_output = ".";
    }
    if (frm.other_weird_reads_output == "")
    {
        frm.other_weird_reads_output = ".";
    }
  
    if (frag_id % 1000000 == 0)
    {
        fprintf(stderr, "grouped %d fragments\n", frag_id);
    }

    return frm; 

}

int calculate_inner_size_and_gap_distance(std::vector <Bcd21> & fragment_bcd21_vector, std::vector <int> & inner_size_count_vector, std::vector <int> & gap_distance_count_vector)
{
    // fragment_bcd21_vector should be sorted by map position 
    std::map <std::string, ReadPairIndex> read_idx_map; 
    int flag;
    ReadPairIndex invalid_index; 
    invalid_index.r1 = -1;
    invalid_index.r2 = -1;

    if (fragment_bcd21_vector.size() < 2){
        return 0;
    }

    for (int i = 0; i < fragment_bcd21_vector.size(); i++)
    {
        if (fragment_bcd21_vector[i].flag & (0x4 + 0x100 + 0x200 + 0x400 + 0x800))
        { 
            continue; // ignore non-primary alignments
        }
        if (read_idx_map.count(fragment_bcd21_vector[i].read_id) == 0){
            read_idx_map[fragment_bcd21_vector[i].read_id] = invalid_index;
        }

        if (fragment_bcd21_vector[i].is_read1()){
            read_idx_map[fragment_bcd21_vector[i].read_id].r1 = i;
        }else{
            read_idx_map[fragment_bcd21_vector[i].read_id].r2 = i;
        }
    }

    int r1_idx, r2_idx; 
    int isize = 0; 
    int shifted_isize = 0;
    for (auto & x: read_idx_map)
    {
        if (x.second.r1 < 0 || x.second.r2 < 0)
        {
            continue; // at least one read is not found. 
        }
        r1_idx = x.second.r1; 
        r2_idx = x.second.r2;
        if (is_weird_read_pair(fragment_bcd21_vector[r1_idx], fragment_bcd21_vector[r2_idx], MAX_INNER_SIZE * 2))
        {
            continue; 
        }
        isize = inner_size(fragment_bcd21_vector[r1_idx], fragment_bcd21_vector[r2_idx]);
        shifted_isize = isize + 2 * MAX_READ_LENGTH;
        if (shifted_isize > MAX_INNER_SIZE - 1){  shifted_isize = MAX_INNER_SIZE - 1; }

        if (shifted_isize > 0)
        {
            inner_size_count_vector[shifted_isize] += 1;
        }
    }

    int gap_distance;
    for (int i = 1; i < fragment_bcd21_vector.size(); i++)
    {
        if (fragment_bcd21_vector[i].read_id == fragment_bcd21_vector[i-1].read_id)
        {
            continue; 
        }
        gap_distance = fragment_bcd21_vector[i].start - fragment_bcd21_vector[i-1].end;
        if (gap_distance > MAX_GAP_DISTANCE){ gap_distance = MAX_GAP_DISTANCE; }
        if (gap_distance > 0)
        {
            gap_distance_count_vector[gap_distance] += 1; 
        }

    }
    return 0;
}

int first_round_bcd21_to_bcd22_file (Settings & global_settings, std::string in_bcd21_file, std::string out_bcd22_file, int length_cut, std::vector <int> & inner_size_count_vector, std::vector <int> & gap_distance_count_vector) 
{

    int min_num_good_reads;
    gzFile in_bcd21_fp;
    FILE * out_bcd22_fp;
    std::string bcd22_header;
    
    if (global_settings.is_wgs) {
        min_num_good_reads = 6;
    } else {
        min_num_good_reads = 3;
    }

    if (global_settings.user_defined_min_num_good_reads_per_fragment > 0){
        min_num_good_reads = global_settings.user_defined_min_num_good_reads_per_fragment;
    }

    global_settings.min_num_good_reads_per_fragment = min_num_good_reads;

    in_bcd21_fp = gzopen(in_bcd21_file.c_str(), "r");
	if (Z_NULL == in_bcd21_fp) {
		fprintf(stderr, "ERROR! Failed to open file for reading: %s\n", in_bcd21_file.c_str());
		exit(1);
	}

    out_bcd22_fp = fopen(out_bcd22_file.c_str(), "w");
    if (NULL == out_bcd22_fp){
        fprintf(stderr, "ERROR! Failed to open file for writing: %s\n", out_bcd22_file.c_str());
		exit(1);
    }
    
    bcd22_header = "#tid\tfrag_start\tfrag_end\tfrag_length\tfrag_barcode\tfrag_ID\tnum_reads\thptype0\thptype1\thptype2\tmap_pos\tmap_qual\tnum_left_weird_reads\tnum_right_weird_reads\tleft_weird_reads_info\tright_weird_reads_info\tother_weird_reads_info\n";
    bcd22_header += "#gap_distance_cut_off=" + std::to_string(length_cut) + "\n";

    fprintf(out_bcd22_fp, "%s", bcd22_header.c_str());
    
   
    char * line = new char[line_max];;
    char * new_bcd = new char[1024]; 
    char * new_read_id = new char[1024];
    char * cigar_string = new char [4096];
    std::vector <Bcd21> fragment_bcd21_vector;

    std::string curr_bcd = "";
    int bcd_cnt = 0;
    int curr_tid = -1;
    frag_id = 0;

    fragment_bcd21_vector.reserve(100);
    
    Bcd22 frm; 
    while (gzgets(in_bcd21_fp, line, line_max))
    {
        if (line[0] == '#') { continue; }
        Bcd21 new_bcd21; 
		convert1line2bcd21(line, new_bcd21, new_bcd, new_read_id, cigar_string);
        if (new_bcd21.flag & (4 + 256 + 1024 + 2048)){ continue; }

        if (fragment_bcd21_vector.size() == 0 || (new_bcd21.bcd == fragment_bcd21_vector.back().bcd and new_bcd21.key_start() - fragment_bcd21_vector.back().key_end() < length_cut) ) {
            fragment_bcd21_vector.push_back(new_bcd21);
        }else{
            std::sort (fragment_bcd21_vector.begin(), fragment_bcd21_vector.end(), sort_bcd21_by_pos);
            frm = convert_bcd21vector_to_bcd22 (global_settings, fragment_bcd21_vector); 
     
            if (frm.num_good_reads >= min_num_good_reads)
            {
                fprintf(out_bcd22_fp, "%s\n", frm.output().c_str());
                calculate_inner_size_and_gap_distance(fragment_bcd21_vector, inner_size_count_vector, gap_distance_count_vector); 
            }
            fragment_bcd21_vector.clear();
            fragment_bcd21_vector.push_back(new_bcd21);
        }
    }
 
    if (fragment_bcd21_vector.size()> 0)
    {
        std::sort (fragment_bcd21_vector.begin(), fragment_bcd21_vector.end(), sort_bcd21_by_pos);
        frm = convert_bcd21vector_to_bcd22(global_settings, fragment_bcd21_vector);

        if (frm.num_good_reads >= min_num_good_reads)
        {
            fprintf(out_bcd22_fp, "%s\n", frm.output().c_str());
            calculate_inner_size_and_gap_distance(fragment_bcd21_vector, inner_size_count_vector, gap_distance_count_vector); 
        }
    }

    delete [] line;
    delete [] new_bcd;
    delete [] new_read_id;
    delete [] cigar_string;
    gzclose(in_bcd21_fp);
    fclose(out_bcd22_fp);

    return 0;
}



int estimate_distribution(Settings & global_settings, std::vector<int> & inner_size_count_vector, std::vector<int> &gap_distance_count_vector)
{

    QuantileNumbers inner_size_quantiles; 
    QuantileNumbers gap_distance_quantiles;

    calculate_distribution_from_count_vector(inner_size_count_vector, inner_size_quantiles);
    calculate_distribution_from_count_vector(gap_distance_count_vector, gap_distance_quantiles);

    double median, rstd; 
    median = inner_size_quantiles.q[500]; 
    rstd = inner_size_quantiles.q[841] - median;
    median -= MAX_READ_LENGTH * 2;

    global_settings.inner_size_cutoff = inner_size_quantiles.q[990] - 2 * MAX_READ_LENGTH;
    global_settings.gap_distance_cutoff = gap_distance_quantiles.q[990];

    using namespace std;

    cerr << "mean and sd of inner size: " << median << ", " << rstd <<  endl;
    cerr << "quantiles of inner size:" << endl;
    inner_size_quantiles.print(stderr, 2 * MAX_READ_LENGTH); 
    cerr << endl;

    cerr << "quantiles of gap distance:" << endl;
    gap_distance_quantiles.print(stderr, 0.0); 
    cerr << endl;

    cerr << "inner size cut-off = " << global_settings.inner_size_cutoff << endl;
    cerr << "gap distance cut-off = " << global_settings.gap_distance_cutoff << endl; 

    return 0; 

}


int find_weird_reads_from_bcd21_vector(Settings & global_settings, std::vector <Bcd21> & onebarcode_bcd21_vector, int start_idx, int end_idx)
{

    int weird_read_flag = 0; 
    for (int i = start_idx; i < end_idx; i++)
    {
        if (onebarcode_bcd21_vector[i].not_primary_alignment() || onebarcode_bcd21_vector[i].mapq < 10)
        {
            continue;
        }
        for (int j = start_idx +1; j < end_idx; j++)
        {
            if (onebarcode_bcd21_vector[j].not_primary_alignment() || onebarcode_bcd21_vector[j].mapq < 10)
            {
                continue;
            }
            weird_read_flag = is_weird_read_pair(onebarcode_bcd21_vector[i], onebarcode_bcd21_vector[j], global_settings.inner_size_cutoff);
            onebarcode_bcd21_vector[i].is_weird_read |= weird_read_flag; 
            onebarcode_bcd21_vector[j].is_weird_read |= weird_read_flag;
        }
    }
    return 0;
}

int process_all_bcd21_of_one_barcode(Settings & global_settings, std::vector <Bcd21> & onebarcode_bcd21_vector, std::vector <Bcd22> & output_frm_vector, FILE * weird_read_fp)
{
    // find all weird reads
    // group fragments (only split read pairs that support deletion )
    // output results

    int start_idx, end_idx; 

    if (onebarcode_bcd21_vector.size() < 2)
    {
        return 0; 
    }

    std::sort (onebarcode_bcd21_vector.begin(), onebarcode_bcd21_vector.end(), sort_bcd21_by_read_id);

    start_idx = 0;
    end_idx = 1;
    while (end_idx < onebarcode_bcd21_vector.size())
    {
        if (onebarcode_bcd21_vector[end_idx].read_id == onebarcode_bcd21_vector[start_idx].read_id)
        {
            end_idx++; 
        }else  {
            find_weird_reads_from_bcd21_vector(global_settings, onebarcode_bcd21_vector, start_idx, end_idx); // all alignments of one pair
            start_idx = end_idx;
            end_idx++;
        }
    }
    if (end_idx > start_idx)
    {
        find_weird_reads_from_bcd21_vector(global_settings, onebarcode_bcd21_vector, start_idx, end_idx); 
    }


    for (size_t i = 0; i < onebarcode_bcd21_vector.size(); i++)
    {
        if (onebarcode_bcd21_vector[i].is_weird_read == 0x8 && onebarcode_bcd21_vector[i].inner_size < global_settings.gap_distance_cutoff && onebarcode_bcd21_vector[i].inner_size > -global_settings.gap_distance_cutoff )
        {
            fprintf(weird_read_fp, "%s\n", onebarcode_bcd21_vector[i].output_Bcd21().c_str());
        }
    }

    std::sort (onebarcode_bcd21_vector.begin(), onebarcode_bcd21_vector.end(), sort_bcd21_by_pos);
    std::vector<Bcd21> one_fragment_bcd21_vector; 
    one_fragment_bcd21_vector.reserve(100);
    one_fragment_bcd21_vector.push_back(onebarcode_bcd21_vector[0]);

    bool cut = false;
    Bcd22 frm; 
    for (int i = 1; i < onebarcode_bcd21_vector.size(); i++)
    {
        if ( ! (onebarcode_bcd21_vector[i].flag & (4 + 256 + 512 + 1024 + 2048)) )
        {
            global_settings.total_num_reads ++; 
        }
        if (onebarcode_bcd21_vector[i].is_weird_read)
        {
            global_settings.total_num_weird_reads ++; 
        }

        if (onebarcode_bcd21_vector[i].read_id == onebarcode_bcd21_vector[i-1].read_id && inner_size(onebarcode_bcd21_vector[i], onebarcode_bcd21_vector[i-1]) > global_settings.inner_size_cutoff) 
        {
            cut = true;
        } else if (onebarcode_bcd21_vector[i].key_start() - onebarcode_bcd21_vector[i-1].key_end() > global_settings.gap_distance_cutoff)
        {
            cut = true;
        } else {
            cut = false;
        }
        if (cut && one_fragment_bcd21_vector.size() > 0)
        {
            frm = convert_bcd21vector_to_bcd22 (global_settings, one_fragment_bcd21_vector);
            if(frm.num_good_reads > global_settings.min_num_good_reads_per_fragment)            {
                output_frm_vector.push_back(frm);
            }
            one_fragment_bcd21_vector.clear();
            one_fragment_bcd21_vector.push_back(onebarcode_bcd21_vector[i]);
        }else{
            one_fragment_bcd21_vector.push_back(onebarcode_bcd21_vector[i]);
        }
    }

    if (one_fragment_bcd21_vector.size() > 0)
    {
        frm = convert_bcd21vector_to_bcd22 (global_settings, one_fragment_bcd21_vector);
        if(frm.num_good_reads > global_settings.min_num_good_reads_per_fragment) {
            output_frm_vector.push_back(frm);
        }
    }

    return 0; 

}


int second_round_bcd21_to_bcd22_file (Settings & global_settings, std::string in_bcd21_file, std::string out_bcd22_file)
{
    
    int min_num_good_reads;
    gzFile in_bcd21_fp;
    FILE * out_bcd22_fp;
    FILE * weird_reads_fp;
    std::string bcd22_header;
    
    if (global_settings.is_wgs) {
        min_num_good_reads = 6;
    } else {
        min_num_good_reads = 3;
    }

    if (global_settings.user_defined_min_num_good_reads_per_fragment > 0){
        min_num_good_reads = global_settings.user_defined_min_num_good_reads_per_fragment;
    }

    global_settings.min_num_good_reads_per_fragment = min_num_good_reads;


    in_bcd21_fp = gzopen(in_bcd21_file.c_str(), "r");
	if (Z_NULL == in_bcd21_fp) {
		fprintf(stderr, "ERROR! Failed to open file for reading: %s\n", in_bcd21_file.c_str());
		exit(1);
	}

    weird_reads_fp = fopen(global_settings.weird_reads_file.c_str(), "w");
    if (NULL == weird_reads_fp)
    {
        fprintf(stderr, "ERROR! Failed to open file for writing: %s\n", global_settings.weird_reads_file.c_str());
		exit(1);
    }

    out_bcd22_fp = fopen(out_bcd22_file.c_str(), "w");
    if (NULL == out_bcd22_fp){
        fprintf(stderr, "ERROR! Failed to open file for writing: %s\n", out_bcd22_file.c_str());
		exit(1);
    }
    
    bcd22_header = "#tid\tfrag_start\tfrag_end\tfrag_length\tfrag_barcode\tfrag_ID\tnum_reads\thptype0\thptype1\thptype2\tmap_pos\tmap_qual\tnum_left_weird_reads\tnum_right_weird_reads\tleft_weird_reads_info\tright_weird_reads_info\tother_weird_reads_info\n";
    bcd22_header += "#gap_distance_cut_off=" + std::to_string(global_settings.gap_distance_cutoff) + "\n";
    bcd22_header += "#inner_size_cutoff=" + std::to_string(global_settings.inner_size_cutoff) + "\n";

    fprintf(out_bcd22_fp, "%s", bcd22_header.c_str());
    
    char * line = new char[line_max];;
    char * new_bcd = new char[1024]; 
    char * new_read_id = new char[1024];
    char * cigar_string = new char [4096];
    std::vector <Bcd21> fragment_bcd21_vector;

    std::string curr_bcd = "";
    int bcd_cnt = 0;
    int curr_tid = -1;
    frag_id = 0;

    fragment_bcd21_vector.reserve(100);
    
    Bcd22 frm; 
    std::vector <Bcd22> output_frm_vector;
    output_frm_vector.reserve(10);
    
    while (gzgets(in_bcd21_fp, line, line_max))
    {
        if (line[0] == '#') { continue; }
        Bcd21 new_bcd21; 
		convert1line2bcd21(line, new_bcd21, new_bcd, new_read_id, cigar_string);
        if (new_bcd21.flag & (4 + 1024)){ continue; } // ingnore unmapped reads and duplicated reads

        if (fragment_bcd21_vector.size() == 0 || new_bcd21.bcd == fragment_bcd21_vector.back().bcd) {
            fragment_bcd21_vector.push_back(new_bcd21);
        }else{
            output_frm_vector.clear();
            process_all_bcd21_of_one_barcode(global_settings, fragment_bcd21_vector, output_frm_vector, weird_reads_fp);
            for (auto & frm : output_frm_vector)
            {
                fprintf(out_bcd22_fp, "%s\n", frm.output().c_str());
            }
            fragment_bcd21_vector.clear();
            fragment_bcd21_vector.push_back(new_bcd21);
        }
    }

    if (fragment_bcd21_vector.size() > 0)
    {
        output_frm_vector.clear();
        process_all_bcd21_of_one_barcode(global_settings, fragment_bcd21_vector, output_frm_vector, weird_reads_fp);
        for (auto & frm : output_frm_vector)
        {
            fprintf(out_bcd22_fp, "%s\n", frm.output().c_str());
        }
        
    }

    delete [] line;
    delete [] new_bcd;
    delete [] new_read_id;
    delete [] cigar_string; 
    
    gzclose(in_bcd21_fp);
    fclose(out_bcd22_fp);
    fclose(weird_reads_fp);
    
    return 0;

}

int cluster_reads(Settings & global_settings)
{

    int length_cut;
    /* first round */
    std::string tmpbcd22_file;
    tmpbcd22_file = global_settings.bcd22_file + ".tmp";
    length_cut = 50 * 1000; // initial length cut value


    std::vector <int> inner_size_count_vector (MAX_INNER_SIZE+1, 0); 
    std::vector <int> gap_distance_count_vector(MAX_GAP_DISTANCE+1, 0);

    first_round_bcd21_to_bcd22_file (global_settings, global_settings.bcd21_file, tmpbcd22_file, length_cut, inner_size_count_vector, gap_distance_count_vector);

    estimate_distribution (global_settings, inner_size_count_vector, gap_distance_count_vector); 

    fprintf(stderr, "started second round\n");
    second_round_bcd21_to_bcd22_file (global_settings, global_settings.bcd21_file, global_settings.bcd22_file);

    fprintf(stderr, "total number of reads is: %d\n", global_settings.total_num_reads);
    fprintf(stderr, "total number of weird reads is: %d\n", global_settings.total_num_weird_reads);

    return 0;
}

int usage(FILE * fp)
{
    fprintf (fp, "Usage: cluster_reads <in.bcd21> <output_bcd22> <output_weird_read_file> <is_wgs> <user_defined_min_num_good_reads_per_fragment> <min_mapq>\n");
    return 0;
}

int main (int argc, char * argv[])
{
    if (argc < 7){
        usage(stderr);
        return 1;
    }

    Settings global_settings; 

    global_settings.bcd21_file = argv[1];
    global_settings.bcd22_file = argv[2];
    global_settings.weird_reads_file = argv[3];
    global_settings.is_wgs = atoi(argv[4]);
    global_settings.min_mapq = atoi(argv[5]);
    global_settings.user_defined_min_num_good_reads_per_fragment = atoi(argv[6]);
    
    cluster_reads(global_settings);

    return 0;

}
