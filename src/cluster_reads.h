# include <iostream>
# include <vector>
# include <string>
# include <unordered_map>
# include <stdint.h>
# include <stdio.h>

const int64_t FIX_LENGTH = (int64_t)(1e10);

class Bcd21{
public:
    int32_t tid, start, end, mapq, flag, n_left_clip, n_right_clip, inner_size, mate_tid, mate_pos, hptype;
    std::string bcd, read_id, cigar_string;
    int32_t is_weird_read;
    Bcd21(){
        is_weird_read = 0; 
    }

    int64_t key_start () const
    {
        return (int64_t) tid * FIX_LENGTH + (int64_t) start;
    }

    int64_t key_end () const
    {
        return (int64_t) tid * FIX_LENGTH + (int64_t) end;
    }

    bool not_primary_alignment() const 
    {
        if (flag & (0x100 + 0x200 + 0x400 + 0x800)){
            return true;
        }else{
            return false;
        }
    }

    int map_direction() const 
    {
        if (flag & 0x10){
            return -1; 
        }else{
            return 1;
        }
    }

    int is_read1() const 
    {
        if (flag & 0x40){
            return true;
        } else {
            return false;
        }
    }

    int is_read2() const 
    {
        if (flag & 0x80){
            return true;
        } else {
            return false;
        }
    }
    std::string support_sv_type() const
    {
        using namespace std;
        string outstring; 
        outstring = "";
        if (is_weird_read & 0x1)
        {
            outstring += "TRA;";
        }else if (is_weird_read & 0x2){
            outstring += "INV;";
        }else if (is_weird_read & 0x4){
            outstring += "DUP;";
        }else if (is_weird_read & 0x8){
            outstring += "DEL;";
        }
        return outstring;

    }
    std::string output_read_info() const
    {
        using namespace std; 

        std::string outstring;
        outstring = "[";
        outstring += read_id + "]" + to_string (tid) + "]" + to_string(start) + "]" + to_string(end) + "]";
        outstring += to_string(mapq) + "]" + to_string (hptype) + "]" + to_string(flag) + "]" + to_string(inner_size) + "]";
        outstring += to_string(mate_tid) + "]" + to_string (mate_pos) + "]";

        return outstring; 
    }

    std::string output_Bcd21() const
    {
        using namespace std; 

        std::string outstring;
        
        outstring = to_string (tid) + "\t" + to_string(start) + "\t" + to_string(end) + "\t" + to_string(mapq) + "\t" + bcd + "\t"; 
        outstring += to_string (hptype) + "\t" + read_id + "\t" + to_string(flag) + "\t" + to_string(n_left_clip) + "\t" + to_string(n_right_clip) + "\t";
        outstring += to_string(inner_size) + "\t" + to_string(mate_tid) + "\t" + to_string (mate_pos) + "\t" + cigar_string;

        return outstring; 
    }
};

class Bcd22{
public: 
    int32_t tid, start, end, frag_id, num_reads, hp0, hp1, hp2, n_left_weird_reads, n_right_weird_reads; 
    std::string bcd, map_pos, map_qual; 
    std::string left_weird_reads_output, right_weird_reads_output, other_weird_reads_output;

    // members that will not be output
    int32_t num_good_reads;

    Bcd22()
    {
        tid = start = end = frag_id = num_reads = -1;
        hp0 = hp1 = hp2 = 0;
        n_left_weird_reads = n_right_weird_reads = 0;
        bcd = map_pos = map_qual = "";
        left_weird_reads_output = right_weird_reads_output = other_weird_reads_output = "";
        num_good_reads = 0;
    }

    int32_t length () const
    {
        return end - start;
    }

    int64_t key_start () const
    {
        return (int64_t) tid * FIX_LENGTH + (int64_t) start;
    }

    int64_t key_end () const
    {
        return (int64_t) tid * FIX_LENGTH + (int64_t) end;
    }

    std::string output() const
    {
        using namespace std; 
        std::string outstring;
        outstring = "";
        outstring += to_string(tid) + "\t" + to_string(start) + "\t" + to_string(end) + "\t" + to_string(length()) + "\t";
        outstring += bcd + "\t" + to_string(frag_id) + "\t" + to_string(num_reads) + "\t";
        outstring += to_string(hp0) + "\t" + to_string(hp1) + "\t" + to_string(hp2) + "\t" + map_pos + "\t" + map_qual + "\t";
        outstring += to_string(n_left_weird_reads) + "\t" + to_string(n_right_weird_reads) + "\t" + left_weird_reads_output + "\t" + right_weird_reads_output + "\t" + other_weird_reads_output;
        return outstring;
    }

};



class Settings {
public:
	int is_wgs;
	std::string output_dir;
	int user_defined_min_num_good_reads_per_fragment;
	int min_num_good_reads_per_fragment;  //min_num_good_reads_per_fragment that is actually used 
	std::string bcd21_file;
	std::string bcd22_file;
	std::string weird_reads_file;
	int min_mapq;
	int inner_size_cutoff; 
    int gap_distance_cutoff;
    int total_num_reads; 
    int total_num_weird_reads;
    
 
	Settings() {
		is_wgs = 1;
		output_dir = "";
		bcd21_file = "";
		bcd22_file = "";
		weird_reads_file = "";
        total_num_reads = 0; 
        total_num_weird_reads = 0; 
		user_defined_min_num_good_reads_per_fragment = -1;
		min_num_good_reads_per_fragment = 6; // default of wgs
		min_mapq = 20;

        if (is_wgs){
            gap_distance_cutoff = 10000;
        }else{
            gap_distance_cutoff = 20000;
        }
	}

};