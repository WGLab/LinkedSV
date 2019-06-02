#include <iostream>
#include <vector>

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <zlib.h>

#include "cluster_reads.h"

const int MAX_PATH_LENGTH = 8192;
const LINE_MAX = 4096;

inline int convert1line2bcd21(const char * line, Bcd21 & bcd21, char * bcd, char * read_id, char * cigar_string)
{

    sscanf(line, "%d\t%d\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%t%s\n", &bcd21.tid, &bcd21.start, &bcd21.pos, &bcd21.mapq, bcd, &bcd21.hptype, read_id, &bcd21.flag, &bcd21.n_left_clip, &bcd21.n_right_clip, &bcd21.insert_size, &bcd21.mate_tid, &bcd21.mate_pos, cigar_string);
    bcd21.bcd = bcd;
    bcd21.read_id = read_id;
    bcd21.cigar_string = cigar_string;
    return 0;    
}

inline int convert_bcd21vector_to_bcd22 (std::vector <std::string> & fragment_bcd21_vector, int frag_id, Bcd22 & bcd22, bool is_fast_mode)
{
    return 0;
}

int bcd21_to_bcd22_file (Settings global_settings, std::string in_bcd21_file, std::string out_bcd22_file, int length_cut, bool is_fast_mode = true) 
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

    in_bcd21_fp = fopen(in_bcd21_file.c_str(), "r");

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
    
   
    bcd22_header = "#tid\tfrag_start\tfrag_end\tfrag_length\tfrag_barcode\tfrag_ID\tnum_reads\thptype0\thptype1\thptype2\tmap_pos\tmap_qual\tnum_left_weird_reads\tnum_right_weird_reads\tleft_weird_reads_info\tright_weird_reads_info\tother_weird_reads_info; gap_distance_cut_off=";
    bcd22_header += std::to_string(length_cut) + "\n";

    fprintf(out_bcd22_fp, bcd22_header.c_str());
    
   
    char * line = new char[LINE_MAX];;
    char * new_bcd = new char[1024]; 
    char * new_read_id = new char[1024];
    char * cigar_string = new char [4096];
    std::vector <Bcd21> fragment_bcd21_vector;

    std::string curr_bcd = "";
    int bcd_cnt = 0;
    int curr_tid = -1;
    int frag_id = 0;
    

    fragment_bcd21_vector.reserve(100);
    
    while (gzgets(bcd21_fp, line, LINE_MAX))
    {
        if (line[0] == '#') { continue; }
        Bcd21 new_bcd21; 
        convert1line_to_bcd21(line, new_bcd21, new_bcd, new_read_id);
        if (new_bcd21.flag & (256 + 1024 + 2048)){ continue; }

        if (fragment_bcd21_vector.size() == 0 || (new_bcd21.bcd == fragment_bcd21_vector.back().bcd and new_bcd21.key_start() - fragment_bcd21_vector.back().key_end() < length_cut) ) {
            fragment_bcd21_vector.push_back(new_bcd21);
        }else{
            Bcd22 bcd22; 
            bcd22 = convert_bcd21vector_to_bcd22 (fragment_bcd21_vector, frag_id, bcd22, is_fast_mode); 
            
            frag_id += 1;
            fragment_bcd21_vector.clear();
            fragment_bcd21_vector.push_back(new_bcd21);
        }

 
		if (total_read_pos_list->size == 0)
		{
			curr_tid = new_tid;
			strcpy(curr_bcd, new_bcd);
			append_int_list(total_read_pos_list, avg_pos);
			if (mapq >= mapq_cutoff) {
				append_int_list(high_mapq_read_pos_list, avg_pos);
			}
			continue;
		}
		if (curr_tid != new_tid || strcmp(curr_bcd, new_bcd) != 0) {
			process_one_barcode_in_one_chr(curr_tid, total_read_pos_list, wg_total_bcd_depth_list, bin_size, idx_list_total);
    }

    delete[] line;
    delete [] new_bcd;
    delete [] new_read_id; 
    delete [] cigar_string; 

    /*

    weird_readname_dict = dict()

    if with_weird_reads == True:
        tid2chrname_list, chrname2tid_dict = get_chrnames(args.faidx_file)
        fid = 1
        weird_reads_file = args.weird_reads_file + '.split%d' % (fid) 
        args.temp_file_list.append(weird_reads_file)
        myprint('getting weird read names from file: %s' % weird_reads_file)
        weird_readname_dict = get_weird_readname_dict (chrname2tid_dict, weird_reads_file)
        myprint('finished getting weird read names')

    bcd21_fp = my_utils.gzopen(bcd21_file, 'r')

    bcd22_fp = open(bcd22_file, 'w')
    bcd22_header = 
    bcd22_fp.write(bcd22_header)

    fragment_bcd21_list = list()

    frag_id = 0
    while 1:
        line = bcd21_fp.readline()
        if not line: break
        if line[0] == '#': continue
        line = line.strip().split(tab)
        new_bcd21_term = Bcd21(line)
        bcd = new_bcd21_term.bcd
        read_id = new_bcd21_term.read_id

        if with_weird_reads == True and bcd in split_barcode_set:

            bcd22 = convert_bcd21list_to_bcd22(args, fragment_bcd21_list, frag_id, is_fast_mode, weird_readname_dict)
            if n_good_reads(bcd22) >= min_num_good_reads: bcd22_fp.write(bcd22.output() + endl)
            frag_id += 1
            fragment_bcd21_list = list()
            fragment_bcd21_list.append(new_bcd21_term)

            fid += 1
            split_barcode_set.remove(bcd)
            weird_reads_file = args.weird_reads_file + '.split%d' % (fid) 
            args.temp_file_list.append(weird_reads_file) 
            del weird_readname_dict
            weird_readname_dict = get_weird_readname_dict (chrname2tid_dict, weird_reads_file)
            continue

        if len(fragment_bcd21_list) == 0 or ( new_bcd21_term.bcd == fragment_bcd21_list[-1].bcd and new_bcd21_term.key_start() - fragment_bcd21_list[-1].key_end() < length_cut):
            fragment_bcd21_list.append(new_bcd21_term)
        else:
            bcd22 = convert_bcd21list_to_bcd22 (args, fragment_bcd21_list, frag_id, is_fast_mode, weird_readname_dict)
            if n_good_reads(bcd22) >= min_num_good_reads: bcd22_fp.write(bcd22.output() + endl)
            frag_id += 1
            fragment_bcd21_list = list()
            fragment_bcd21_list.append(new_bcd21_term)

            if frag_id % 1000000 == 0: myprint('grouped %d fragments' % frag_id)

    if len(fragment_bcd21_list) > 0:
        bcd22 = convert_bcd21list_to_bcd22 (args, fragment_bcd21_list, frag_id, is_fast_mode, weird_readname_dict)
        if n_good_reads(bcd22) >= min_num_good_reads: bcd22_fp.write(bcd22.output() + endl)

    bcd21_fp.close()
    bcd22_fp.close()

    if with_weird_reads == True: del weird_readname_dict

    gc.collect()

    return
    */
    return 0;
}

int cluster_reads(Settings global_settings)
{

    int length_cut;
    bool is_fast_mode;


    /* first round */
    std::string tmpbcd22_file;
    tmpbcd22_file = global_settings.bcd22_file + ".tmp";
    length_cut = 50 * 1000; // initial length cut value
    is_fast_mode = true;

    bcd21_to_bcd22_file (global_settings.bcd21_file, tmpbcd22_file, length_cut, is_wgs, is_fast_mode);

    delete [] tmpbcd22_file ;
    return 0;
}

int usage(FILE * fp)
{
    fprintf (fp, "Usage: cluster_reads <in.bcd21> <output_bcd22> <is_wgs> <user_defined_min_num_good_reads_per_fragment>\n");
    return 0;
}

int main (int argc, char * argv[])
{
    if (argc < 4){
        usage(stderr);
        return 1;
    }

    Settings global_settings; 

    char * in_bcd21_file;
    char * out_bcd22_file;
    int is_wgs;
    int min_num_good_reads;
    
    global_settings.bcd21_file = argv[1];
    global_settings.bcd22_file = argv[2];
    global_settings.is_wgs = atoi(argv[3]);
    global_settings.user_defined_min_num_good_reads_per_fragment = atoi(argv[4]);
    
    cluster_reads(global_settings);

    return 0;

}