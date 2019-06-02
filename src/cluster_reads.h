# include <iostream>
# include <vector>
# include <stdint.h>



const int64_t FIX_LENGTH = (int64_t)(1e10);


class Bcd21{
public:
    int tid;
    int start, end;
    int mapq;
    int hptype; // 3 valide values: 0, 1, 2
    std::string bcd;
    std::string read_id;
    int flag;
    int n_left_clip;
    int n_right_clip;
    int insert_size;
    int mate_tid;
    int mate_pos;
    std::string cigar_string; 

    int64_t key_start (){
        return (int64_t) tid * FIX_LENGTH + (int64_t) start;
    }
    int64_t key_end (){
        return (int64_t) tid * FIX_LENGTH + (int64_t) end;
    }

};

class Bcd22{
    
}

class Settings{
public:
    int is_wgs;
    std::string output_dir; 
    int user_defined_min_num_good_reads_per_fragment; 
    int min_num_good_reads_per_fragment;  //min_num_good_reads_per_fragment that is actually used 
    std::string bcd21_file; 
    std::string bcd22_file;
    Settings(){
        is_wgs = 1;
        output_dir = "";
        bcd21_file = "";
        bcd22_file = "";
        user_defined_min_num_good_reads_per_fragment = -1;
        min_num_good_reads_per_fragment = 6; // default of wgs

    }
}


/*
class Bcd21:

    def __init__(self, attr_list):
        #barcode    tid start_post  end_pos map_qual    hap_type    ReadID  flag    n_left_clip n_right_clip    insert_size mate_tid    mate_pos
        self.tid, self.start, self.end, self.mapq, self.bcd, self.hptype, self.read_id, self.flag, self.n_left_clip, self.n_right_clip, self.insert_size, self.mate_tid, self.mate_pos = attr_list[0:13]
        self.tid = int(self.tid)
        self.start = int(self.start)
        self.end = int(self.end)
        self.mapq = int(self.mapq)
        self.hptype = int(self.hptype)
        self.flag = int(self.flag)
        self.n_left_clip = int(self.n_left_clip)
        self.n_right_clip = int(self.n_right_clip)
        self.insert_size = int(self.insert_size)
        self.mate_tid = int(self.mate_tid)
        self.mate_pos = int(self.mate_pos)

    def key_start(self):
        return self.tid * FIX_LENGTH + self.start

    def key_end(self):
        return self.tid * FIX_LENGTH + self.end

    def output_read_info (self):
        output = '[%s]%d]%d]%d]' % (self.read_id, self.tid, self.start, self.end) 
        output += '%d]%d]%d]%d]' % (self.mapq, self.hptype, self.flag, self.insert_size)
        output += '%d]%d]' % (self.mate_tid, self.mate_pos)
        return output 

*/
