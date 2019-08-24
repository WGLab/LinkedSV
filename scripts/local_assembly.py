#!/usr/bin/env python

import os
import sys
import my_utils
import time
import multiprocessing

tab = '\t'
endl = '\n'
arg = sys.argv[1:]

usage = 'python ' + __file__ + ' ' + '<input_bam_file> <out_dir> <out_del_call_file> <n_threads>  <ref_fasta> <fermikit_dir> <samtools> <bedtools> '


argc = 0
for i in range(0, len(usage)):
    if usage[i]  == '<':
        argc += 1

class Interval:
    def __init__(self, chrom, start_pos, end_pos):
        self.chrom = chrom
        self.start_pos = int(start_pos)
        self.end_pos = int(end_pos)
        self.out_dir = ''
        self.shell_file = ''
        self.shell_cmds = ''
        self.region_bed_file = ''


def main():

    if len(arg) < argc: 
        print (usage)
        sys.exit()

    input_bam_file = os.path.abspath(arg.pop(0))
    out_dir        = os.path.abspath(arg.pop(0))
    out_del_call_file     = arg.pop(0)
    n_threads      = int(arg.pop(0))
    ref_fasta_file = os.path.abspath(arg.pop(0))  
    fermikit_dir   = os.path.abspath(arg.pop(0)) 
    samtools       = os.path.abspath(arg.pop(0))
    bedtools       = os.path.abspath(arg.pop(0)) 
    window_size    = int(2e5)
    faidx_file     = ref_fasta_file + '.fai'
    max_depth      = 500

    small_deletion_dection_by_local_assembly(samtools, bedtools, fermikit_dir, input_bam_file, ref_fasta_file, faidx_file, out_dir, out_del_call_file, n_threads, window_size, max_depth)

    return

def small_deletion_dection_by_local_assembly(samtools, bedtools, fermikit_dir, input_bam_file, ref_fasta_file, faidx_file, out_dir, out_del_call_file, n_threads, window_size, max_depth, rm_temp_files = 1):

    
    if os.path.exists(faidx_file) == False:
        cmd = '%s faidx %s' % (samtools, ref_fasta_file)
        my_utils.myprint(cmd)
        os.system(cmd)

    if os.path.exists(faidx_file) == False:
        my_utils.myprint ('ERROR! The index file of the reference fasta file does not exist!' )
        sys.exit()

    cmd = 'mkdir -p %s' % out_dir
    my_utils.myprint(cmd)
    os.system('mkdir -p %s' % out_dir)
    tid2chrname_list, chrname2tid_dict = my_utils.get_chrnames(faidx_file)
    chr_len_list = my_utils.get_chr_length(faidx_file)

    overlap_length = int(window_size/10)
    interval_list  = generate_interval_list(chr_len_list, tid2chrname_list, chrname2tid_dict, window_size, overlap_length)

    process_list = list()
    out_combined_vcf_file_list = list()
    for i in range(0, n_threads):
        out_combined_vcf_file = os.path.join(out_dir, 'assembly_raw_variants.%d.txt' % i)
        out_combined_vcf_file_list.append(out_combined_vcf_file)
        t = multiprocessing.Process(target=small_deletion_dection_from_interval_list, args=(i, n_threads, samtools, bedtools, fermikit_dir, input_bam_file, ref_fasta_file, out_dir, window_size, max_depth, interval_list, out_combined_vcf_file))
        process_list.append(t)
        t.start()

    for t in process_list:
        t.join()
    

    all_processes_out_combined_vcf_file = os.path.join(out_dir, 'local_assembly_raw_variants.txt')
    cmd = 'cat ' 
    for out_combined_vcf_file in out_combined_vcf_file_list:
        cmd += ' %s ' % out_combined_vcf_file

    cmd += ' > %s ' % all_processes_out_combined_vcf_file
    my_utils.myprint(cmd)
    os.system(cmd)


    extract_del_from_vcf_file(all_processes_out_combined_vcf_file, out_del_call_file)

    if rm_temp_files:
        for out_combined_vcf_file in out_combined_vcf_file_list:
            os.remove(out_combined_vcf_file)
        os.remove(all_processes_out_combined_vcf_file)

    return

def small_deletion_dection_from_interval_list(thread_id, n_threads, samtools, bedtools, fermikit_dir, input_bam_file, ref_fasta_file, out_dir, window_size, max_depth, interval_list, out_combined_vcf_file):


    out_combined_vcf_fp = open(out_combined_vcf_file, 'w')
    out_combined_vcf_fp.write('')
    out_combined_vcf_fp.close()

    for region_id in range(0, len(interval_list)):
        if region_id % n_threads != thread_id: continue
        itv = interval_list[region_id]
        process1region(samtools, bedtools, fermikit_dir, ref_fasta_file, input_bam_file, out_dir, itv, region_id, window_size, max_depth, 1, out_combined_vcf_file)


    return
 
def process1region(samtools, bedtools, fermikit_dir, ref_fasta_file, input_bam_file, out_dir, itv, region_id, window_size, max_depth, n_threads_for_one_process, out_combined_vcf_file): 

    curr_out_dir = os.path.join(out_dir, 'region_%06d' % (region_id))
    out_bam_file = os.path.join(curr_out_dir, 'region_%06d.bam' % region_id)
    out_all_fastq_file = os.path.join(curr_out_dir, 'region_%06d.all.fastq' % region_id)
    region_bed_file = os.path.join(curr_out_dir, 'region_%06d.bed' % region_id)
    region_fasta_file = os.path.join(curr_out_dir, 'region_%06d.fasta' % region_id)

    interval = '%s:%d-%d' % (itv.chrom, itv.start_pos+1, itv.end_pos)

    cmd = 'mkdir -p %s' % curr_out_dir
    my_utils.myprint(cmd)
    os.system(cmd)
    time.sleep(0.05)
    
    if os.path.exists(curr_out_dir) == False:
        os.system(cmd)
        time.sleep(1)

    if os.path.exists(curr_out_dir) == False:
        my_utils.myprint('Failed to creat directory: %s' % curr_out_dir)
        cmd = 'rm -rf %s' % curr_out_dir
        my_utils.myprint(cmd)
        os.system(cmd)
        return

    cmd = extract_bam_region(samtools, input_bam_file, interval, out_bam_file, n_threads_for_one_process)
    my_utils.myprint(cmd)
    os.system(cmd)
    cmd = index_bam(samtools, out_bam_file)
    my_utils.myprint(cmd)
    os.system(cmd)
    cmd = bam_to_1fastq(samtools, out_bam_file, out_all_fastq_file)
    my_utils.myprint(cmd)
    os.system(cmd)

    fastq_file_size = os.path.getsize(out_all_fastq_file)
    if fastq_file_size > window_size * max_depth * 2 or fastq_file_size < 200: 
        cmd = 'rm -r %s' % curr_out_dir
        my_utils.myprint(cmd)
        os.system(cmd)
        return

    region_bed_fp = open(region_bed_file, 'w')
    region_bed_fp.write('%s\t%d\t%d\n' % (itv.chrom, itv.start_pos, itv.end_pos))
    region_bed_fp.close()

    cmd = extract_ref_region(bedtools, ref_fasta_file, region_bed_file, region_fasta_file)
    my_utils.myprint(cmd)
    os.system(cmd)


    out_prefix = os.path.join(curr_out_dir, 'region_%06d.all_hap' % region_id)

    fermikit_variant_calling(fermikit_dir, samtools, n_threads_for_one_process, region_fasta_file, window_size, out_all_fastq_file, curr_out_dir, out_prefix)

    indel_call_file = out_prefix + '.flt.vcf'
    sv_call_file    = out_prefix + '.sv.vcf'

    cmd = 'gunzip --force %s.gz' % indel_call_file
    my_utils.myprint(cmd)
    os.system(cmd)
    cmd = 'gunzip --force %s.gz' % sv_call_file
    my_utils.myprint(cmd)
    os.system(cmd)

    cmd = 'cat %s %s >> %s' % (indel_call_file, sv_call_file, out_combined_vcf_file)
    my_utils.myprint(cmd)
    os.system(cmd)

    cmd = 'rm -r %s' % curr_out_dir
    my_utils.myprint(cmd)
    os.system(cmd)

    
    return


def fermikit_variant_calling(fermikit_dir, samtools, n_threads_for_one_process, region_fasta_file, window_size, input_fastq_file, curr_out_dir, out_prefix):


    out_mak_file = os.path.join(curr_out_dir, '%s.mak' % out_prefix)
    assembly_contigs_file = os.path.join(curr_out_dir, '%s.mag.gz' % out_prefix)

    cmd = '%s/bwa index %s' % (fermikit_dir, region_fasta_file)
    my_utils.myprint(cmd)
    os.system(cmd)
    cmd = 'perl %s/fermi2.pl unitig -s %s -l 151 -t %d -p %s %s > %s\n\n' % (fermikit_dir, window_size, n_threads_for_one_process, out_prefix,  input_fastq_file, out_mak_file)
    my_utils.myprint(cmd)
    os.system(cmd)

    cmd = 'make -f %s\n\n' % out_mak_file
    my_utils.myprint(cmd)
    os.system(cmd)

    cmd = 'perl %s/run-calling -t %d %s %s | sh \n\n' % (fermikit_dir, n_threads_for_one_process, region_fasta_file, assembly_contigs_file)
    my_utils.myprint(cmd)
    os.system(cmd)
 
    return 

def extract_bam_region(samtools, input_bam, interval, output_bam, n_threads_for_one_process):

    cmd = 'time %s view -1 -hb -@ %d %s %s > %s' % (samtools, n_threads_for_one_process, input_bam, interval, output_bam)
    return cmd

def index_bam(samtools, input_bam):

    cmd = 'time %s index %s' % (samtools, input_bam)
    return cmd

def bam_to_1fastq(samtools, input_bam, out_fastq):

    cmd = 'time %s fastq %s > %s' % (samtools, input_bam, out_fastq)
    return cmd 

def split_bam_by_hap_type(split_hap_type_bam, in_bam_file, out_hap0_bam_file, out_hap1_bam_file, out_hap2_bam_file, out_unmapped_bam_file):

    cmd = 'time %s %s %s %s %s %s ' % (split_hap_type_bam, in_bam_file, out_hap0_bam_file, out_hap1_bam_file, out_hap2_bam_file, out_unmapped_bam_file)
    return cmd

def extract_ref_region(bedtools, ref_fasta_file, bed_file, out_fasta_file):

    cmd = '%s getfasta -fi %s -bed %s -fo %s' % (bedtools, ref_fasta_file, bed_file, out_fasta_file )

    return cmd

        
def generate_interval_list(chr_len_list, tid2chrname_list, chrname2tid_dict, window_size, overlap_length):
    
    interval_list = list()

    step_length = window_size - overlap_length

    for tid in range(0, len(chr_len_list)):
        chrom = tid2chrname_list[tid]
        chr_len = chr_len_list[tid]
        for start_pos in range(0, chr_len, step_length):
            end_pos = start_pos + window_size
            if end_pos > chr_len: 
                end_pos = chr_len
                break
            itv = Interval(chrom, start_pos, end_pos)
            interval_list.append(itv)
            if end_pos == chr_len:
                break

    return interval_list 


def extract_del_from_vcf_file(in_vcf_file, out_file):
    
    in_vcf_fp = open(in_vcf_file, 'r')
    out_fp    = open(out_file, 'w')
    min_del_size = 50
    while 1:
        line = in_vcf_fp.readline().strip()
        if not line: break
        if line[0] == '#' : continue

        items = line.split('\t')

        chrom1 = items[0]
        pos1 = int(items[1])
        ref_allele = items[3]
        alt_allele = items[4]
        flt = items[6]
        info = items[7]

        sv_type = ''
        sv_size = 0

        pos2 = -1

        if '[' in alt_allele or ']' in alt_allele: continue

        ref_chr, ref_start_end = chrom1.split(':')
        ref_start, ref_end = ref_start_end.split('-')
        ref_start = int(ref_start)
        chrom1 = ref_chr
        pos1 += ref_start


        if len(ref_allele) > min_del_size and len(ref_allele) - len(alt_allele) > min_del_size:
            sv_type = 'DEL'
            sv_size = len(ref_allele) - len(alt_allele)
            pos2 = pos1 + sv_size
        else:
            for ele in info.split(';'):
                key = ele.split('=')[0]
                if key == 'SVTYPE':
                    sv_type = ele.split('=')[1]
                elif key == 'SVLEN':
                    sv_size = abs(int(ele.split('=')[1]))
                elif key == 'END' and pos2 == -1:
                    pos2 = int(ele.split('=')[1]) + ref_start

        if sv_type != 'DEL': continue
        
        chrom2 = chrom1
        flt = 'PASS'
        score = 200.0

        out_item = '%s\t%d\t%d\t%s\t%d\t%d\t' % (chrom1, pos1, pos1+1, chrom2, pos2, pos2+1)
        out_item += '%s\t%s\t%d\t%.2f\tPRECISE;SVMETHOD=local_assembly\n' % (sv_type, flt, sv_size, score)

        out_fp.write(out_item)

    in_vcf_fp.close()
    out_fp.close()

    return
 


if __name__ == '__main__':
    main()
