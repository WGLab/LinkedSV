#!/usr/bin/env python

from short_read_info import *
from my_utils import *

class AlignedRead:
    def __init__(self, attr_list):

        self.read_id, self.flag, self.chrm, self.left_pos, self.mapq, self.cigar_string, self.tid = attr_list[0:7]

        self.flag = int(self.flag)
        self.left_pos = int(self.left_pos)
        self.mapq = int(self.mapq)
        self.tid = int(self.tid)

        self.right_pos = None
        self.cigar_operation_list = None
        self.cigar_operation_length_list = None
        self.n_cigar = None
        self.endR_pos = None
        self.cigar_analysis()
        self.get_right_ref_pos() 
        self.get_endR_pos()

    def output_core(self):
        outstring = '%s\t%d\t%s\t%d\t%s\t%d' % (self.read_id, self.flag, self.chrm, self.left_pos, self.cigar_string, self.n_cigar)
        return outstring
        
    def output_core_at(self):
        outstring = '%s@%d@%s@%d@%s' % (self.read_id, self.flag, self.chrm, self.left_pos, self.cigar_string)
        return outstring

    def map_orientation(self):
        if self.flag & 0x10:
            return -1
        else:
            return 1 

    def cigar_analysis(self):
        p1 = p2 = 0 
        number_list = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
        n_cigar = 0 
        self.cigar_operation_list = list()
        self.cigar_operation_length_list = list()
        for i in range(0, len(self.cigar_string)):
            c = self.cigar_string[i]
            if c in number_list: continue
            p2 = i 
            cigar_operation = self.cigar_string[p2]
            self.cigar_operation_list.append(cigar_operation)
            cigar_operation_length = int(self.cigar_string[p1:p2])
            self.cigar_operation_length_list.append(cigar_operation_length)
            p1 = p2 + 1 
            p2 = p2 + 1 
    
        self.n_cigar = len(self.cigar_operation_list)

    def get_right_ref_pos(self):
        if self.cigar_operation_list == None:
            self.cigar_analysis()

        ref_align_length = 0
        affect_ref_align_cigar = set(['M', 'D', 'N', '=', 'X'])
        for i in range(0, self.n_cigar):
            if self.cigar_operation_list[i] in affect_ref_align_cigar:
                ref_align_length += self.cigar_operation_length_list[i]
        
        self.right_pos = self.left_pos + ref_align_length - 1 
        
    def has_left_clip(self):

        if self.n_cigar == 0 or self.n_cigar == None: 
            return False

        if self.cigar_operation_list[0] == 'S' or self.cigar_operation_list[0] == 'H':
            return True
        else:
            return False

    def has_right_clip(self): 
        if self.n_cigar == 0 or self.n_cigar == None: 
            return False

        if self.cigar_operation_list[-1] == 'S' or self.cigar_operation_list[-1] == 'H':
            return True
        else:
            return False

    def has_clip(self):
        return (self.has_left_clip() or self.has_right_clip())

    def get_endR_pos(self):
        if self.flag & 0x10:
            self.endR_pos = self.left_pos 
        else:
            self.endR_pos = self.right_pos

    def major_clip_side(self):

        if self.has_left_clip() == True and self.has_right_clip() == False: 
            return 'left'
        elif self.has_left_clip() == False and self.has_right_clip() == True: 
            return 'right'
        else:
            left_clip_length = self.cigar_operation_length_list[0]
            right_clip_length = self.cigar_operation_length_list[-1]
            if left_clip_length > right_clip_length:
                return 'left'
            elif left_clip_length < right_clip_length:
                return 'right'
            else:
                randnumber = random.randint(1,2)
                if randnumber == 1:
                    return 'left'
                else:
                    return 'right'

class PairedEndSupport:

    def __init__(self, aligned_read1 = None, aligned_read2 = None):

        self.aligned_read1 = None
        self.aligned_read2 = None
        self.endtype1 = None
        self.endtype2 = None

        if aligned_read1 != None:
            self.aligned_read1 = aligned_read1
            self.aligned_read2 = aligned_read2

            if self.aligned_read1.map_orientation() == 1:
                self.endtype1 = 'R_end'
            else: 
                self.endtype1 = 'L_end'

            if self.aligned_read2.map_orientation() == 1:
                self.endtype2 = 'R_end'
            else: 
                self.endtype2 = 'L_end'

    def output(self):

        chr1 = self.aligned_read1.chrm
        chr2 = self.aligned_read2.chrm

        if self.endtype1 == 'R_end':
            pos1 = self.aligned_read1.endR_pos
        else:
            pos1 = self.aligned_read1.left_pos

        if self.endtype2 == 'R_end':
            pos2 = self.aligned_read2.endR_pos
        else:
            pos2 = self.aligned_read2.left_pos

        min_mapq = min(self.aligned_read1.mapq, self.aligned_read2.mapq)

        outstring = '%s\t%d\t%s\t%d\t%s\t%s\t%d\t%s\t%s\t%s' % (chr1, pos1, chr2, pos2, self.endtype1, self.endtype2, min_mapq, 'PE', self.aligned_read1.output_core_at(), self.aligned_read2.output_core_at())

        return outstring

class SplitReadSupport:

    def __init__(self, aligned_read1 = None, aligned_read2 = None): 
        self.aligned_read1 = None
        self.aligned_read2 = None
        self.endtype1 = None
        self.endtype2 = None
        self.major_clip_pos1 = None
        self.major_clip_pos2 = None

        if aligned_read1 != None:
            self.aligned_read1 = aligned_read1
            self.aligned_read2 = aligned_read2

            if self.aligned_read1.major_clip_side() == 'left': 
                self.endtype1 = 'L_end'
            else:
                self.endtype1 = 'R_end'

            if self.aligned_read2.major_clip_side() == 'left': 
                self.endtype2 = 'L_end'
            else:
                self.endtype2 = 'R_end'

            if self.aligned_read1.major_clip_side() == 'left':
                self.major_clip_pos1 = self.aligned_read1.left_pos
            else:
                self.major_clip_pos1 = self.aligned_read1.right_pos

            if self.aligned_read2.major_clip_side() == 'left':
                self.major_clip_pos2 = self.aligned_read2.left_pos
            else:
                self.major_clip_pos2 = self.aligned_read2.right_pos

    def output(self):

        chr1 = self.aligned_read1.chrm
        chr2 = self.aligned_read2.chrm

        pos1 = self.major_clip_pos1
        pos2 = self.major_clip_pos2

        min_mapq = min(self.aligned_read1.mapq, self.aligned_read2.mapq)

        outstring = '%s\t%d\t%s\t%d\t%s\t%s\t%d\t%s\t%s\t%s' % (chr1, pos1, chr2, pos2, self.endtype1, self.endtype2, min_mapq, 'SR', self.aligned_read1.output_core_at(), self.aligned_read2.output_core_at())

        return outstring
        

class SamCore:
    def __init__(self, attr_list):
        self.read_id, self.flag, self.chrm, self.pos, self.mapq, self.cigar, self.tid = attr_list[0:7]
        self.flag = int(self.flag)
        self.pos = int(self.pos)
        self.mapq = int(self.mapq)

    def key(self):
        return  FIX_LENGTH * self.tid + self.pos

    def output(self):
        outstring = '%s\t%d\t%s\t%d\t%d\t%s' % (self.read_id, self.flag, self.chrm, self.pos, self.mapq, self.cigar)
        return outstring
        
def main():

    args, dbo_args, endpoint_args = parse_user_arguments()

    extract_weird_reads(args, dbo_args, endpoint_args) 

    return

def extract_weird_reads(args, dbo_args, endpoint_args): 

    ## sorting bam by read name ##
    task = 'sorting bam by read name'
    cmd = '%s %s | %s sort -m 2G -@ %d -n -o %s -' % (args.output_bam_coreinfo, args.bam, args.samtools, args.n_thread, args.sortn_bam)
    if check_file_exists(args.sortn_bam):
        myprint('File: %s existed, skipped %s' % (args.sortn_bam, task))
    else:
        myprint(task)
        myprint(cmd)
        os.system(cmd)

    task = 'extracting read info'
    n_compress_threads = args.n_thread - 1
    if n_compress_threads < 1: n_compress_threads = 1
    cmd = "%s view %s | cut -f 1-6 | pigz --fast --processes %d  > %s " % (args.samtools, args.sortn_bam, n_compress_threads, args.sortn_bam_core_file)
    if args.run_from_begining == False and check_file_exists(args.sortn_bam_core_file):
        myprint('File: %s existed, skipped %s' % (args.sortn_bam_core_file, task))
    else:
        myprint(task)
        myprint(cmd)
        os.system(cmd)

    ## removing temp bam file ##
    if args.rm_temp_files and check_file_exists(args.sortn_bam_core_file):
        if os.path.exists(args.sortn_bam): os.remove(args.sortn_bam)

    # extracting weird reads from core file #
    task = 'extracting weird reads'
    args.weird_reads_file
    if args.run_from_begining == False and check_file_exists(args.weird_reads_file):
        myprint('File: %s existed, skipped %s' % (args.weird_reads_file, task))
    else:
        myprint(task)

        max_distance = 1500
        min_mapq = args.min_mapq

        alt_chr_name_set = args.alt_chr_name_set
        tid2chrname_list, chrname2tid_dict = get_chrnames (args.faidx_file)

        if args.sortn_bam_core_file[-2:] == 'gz':
            sortn_bam_coreinfo_fp = gzip.open(args.sortn_bam_core_file, 'r')
        else:
            sortn_bam_coreinfo_fp = open(args.sortn_bam_core_file, 'r')

        out_fp = open(args.weird_reads_file, 'w')

        current_read_id = ''
        sam_core_list = list()

        while 1:
            line = sortn_bam_coreinfo_fp.readline()
            if not line: break
            line = line.strip().split(tab)
            if line[2] in alt_chr_name_set: continue
            if line[2] not in chrname2tid_dict: continue
            tid = chrname2tid_dict[line[2]]
            sam_core = SamCore(line + [tid])
            #if sam_core.mapq < min_mapq: continue

            if sam_core.read_id == current_read_id or current_read_id == '': 
                sam_core_list.append(sam_core)
                if current_read_id == '': current_read_id = sam_core.read_id
            else:
                weird_read_info = find_weird_reads(sam_core_list, max_distance)

                if weird_read_info != '': out_fp.write(weird_read_info)

                current_read_id = sam_core.read_id
                sam_core_list = list()
                sam_core_list.append(sam_core)
        out_fp.close()

    myprint ('finished extracting weird reads')
    ## removing bam_coreinfo file ##
    if args.rm_temp_files and check_file_exists(args.weird_reads_file):
        if os.path.exists(args.sortn_bam_core_file): os.remove(args.sortn_bam_core_file)

    return

def find_weird_reads(sam_core_list, max_distance):
    if len(sam_core_list) == 2 and abs(sam_core_list[0].key() - sam_core_list[1].key()) < max_distance: 
        return '' 

    if len(sam_core_list) < 2: return '' 

    sam_core_list.sort(key=lambda x: x.key())

    return_string = ''
    for i in range(1, len(sam_core_list)):
        d = abs(sam_core_list[i].key() - sam_core_list[i-1].key())
        support_read_info_string = get_support_read_info(sam_core_list[i], sam_core_list[i-1])
        if support_read_info_string != '': 
            return_string += support_read_info_string + endl

    return return_string

def get_support_read_info(sam_core1, sam_core2):

    if sam_core1.key() > sam_core2.key():
        tmp = sam_core1
        sam_core1 = sam_core2 
        sam_core2 = tmp

    attr_list1 = [sam_core1.read_id, sam_core1.flag, sam_core1.chrm, sam_core1.pos, sam_core1.mapq, sam_core1.cigar, sam_core1.tid]
    attr_list2 = [sam_core2.read_id, sam_core2.flag, sam_core2.chrm, sam_core2.pos, sam_core2.mapq, sam_core2.cigar, sam_core2.tid]

    aligned_read1 = AlignedRead(attr_list1)
    aligned_read2 = AlignedRead(attr_list2)

    pe_support = None
    sr_support = None

    if aligned_read1.flag & 0x40 != aligned_read2.flag & 0x40:
        pe_support = PairedEndSupport(aligned_read1, aligned_read2)
        return pe_support.output()
    else:
        if aligned_read1.has_clip() and aligned_read2.has_clip():
            sr_support = SplitReadSupport(aligned_read1, aligned_read2)
            return sr_support.output()
        else: 
            return ''

if __name__ == '__main__':
    main()

