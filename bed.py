#!/usr/bin/env python

tab  = '\t'
endl = '\n'

class Bed:
    def __init__(self, chrm, start, end, chrname2tid = None):
        self.chrm = chrm
        self.start = int(start)
        self.end = int(end) 
        if chrname2tid != None and self.chrm in chrname2tid:
            self.tid = chrname2tid[self.chrm]
        else:
            self.tid = -1

    def output(self):
        return "%s\t%d\t%d" %(self.chrm, self.start, self.end)

def read_bed_file(in_bed_file, chrname2tid = None):

    in_bed_fp = open(in_bed_file, 'r')
    lines = list(in_bed_fp)
    in_bed_fp.close()

    bed_list = list()
    for line in lines:
        line = line.strip().split(tab)
        bed = Bed(line[0], line[1], line[2], chrname2tid)
        bed_list.append(bed)

    return bed_list

def calculate_bed_length(in_bed_file):

    total_length = 0 
    bed_list = read_bed_file(in_bed_file) 
    for bed in bed_list:
        total_length += bed.end - bed.start

    return total_length
    

def output_bed_list(bed_list, out_bed_file):
    
    out_bed_fp = open(out_bed_file, 'r')
    for bed in bed_list:
        out_bed_fp.write(bed.output() + endl)

    out_bed_fp.close()

    return


def bed2tidbed_file(in_bed_file, chrname2tid_dict, out_tidbed_file):
    
    bed_list = read_bed_file(in_bed_file)
    out_tidbed_fp = open(out_tidbed_file, 'w')
    for bed in bed_list:
        if bed.chrm not in chrname2tid_dict: continue
        tid = chrname2tid_dict[bed.chrm]
        out_tidbed_fp.write('%d\t%d\t%d\n' % (tid, bed.start, bed.end))

    out_tidbed_fp.close()

