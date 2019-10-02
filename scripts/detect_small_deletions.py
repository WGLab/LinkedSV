#!/usr/bin/env python

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import sys

try:
    from scripts import local_assembly
except ImportError:
    import local_assembly

try:
    from scripts import cluster_weird_reads
except ImportError:
    import cluster_weird_reads

try:
    from scripts import svtk
except ImportError:
    import svtk

try:
    from scripts import my_utils
except ImportError:
    import my_utils

tab  = '\t'
endl = '\n'
arg = sys.argv[1:]

usage = 'python ' + __file__ + ' ' + '<input_bam_file> <out_dir> <out_del_call_file> <n_threads> <ref_fasta_file> <fermikit_dir> <samtools> <bedtools> <in_weird_reads_file> <weird_reads_cluster_file> <call_small_deletions> <cal_hap_read_depth_from_bcd21> <bcd21_file> <bcd22_file> <hap_type_read_depth_file> <gap_region_bed_file>'
argc  = 16


def main():


    if len(arg) < argc:
        print (usage)
        sys.exit()

    input_bam_file =os.path.abspath(arg.pop(0))
    out_dir =os.path.abspath(arg.pop(0))
    out_del_call_file =os.path.abspath(arg.pop(0))
    n_threads =int(arg.pop(0))
    ref_fasta_file =os.path.abspath(arg.pop(0))
    fermikit_dir =os.path.abspath(arg.pop(0))
    samtools =os.path.abspath(arg.pop(0))
    bedtools =os.path.abspath(arg.pop(0))
    in_weird_reads_file =os.path.abspath(arg.pop(0))
    weird_reads_cluster_file =os.path.abspath(arg.pop(0))
    call_small_deletions =os.path.abspath(arg.pop(0))
    cal_hap_read_depth_from_bcd21 =os.path.abspath(arg.pop(0))
    bcd21_file =os.path.abspath(arg.pop(0))
    bcd22_file =os.path.abspath(arg.pop(0))
    hap_type_read_depth_file =os.path.abspath(arg.pop(0))
    gap_region_bed_file =os.path.abspath(arg.pop(0))


    detect_small_deletions(input_bam_file, out_dir, out_del_call_file, n_threads, ref_fasta_file, fermikit_dir, samtools, bedtools, in_weird_reads_file, weird_reads_cluster_file, call_small_deletions, cal_hap_read_depth_from_bcd21, bcd21_file, bcd22_file, hap_type_read_depth_file, gap_region_bed_file)

    return


def detect_small_deletions(input_bam_file, out_dir, out_del_call_file, n_threads, ref_fasta_file, fermikit_dir, samtools, bedtools, in_weird_reads_file, weird_reads_cluster_file, call_small_deletions, cal_hap_read_depth_from_bcd21, bcd21_file, bcd22_file, hap_type_read_depth_file, gap_region_bed_file, rm_temp_files = 1):

    faidx_file     = ref_fasta_file + '.fai'
    window_size    = int(2e5)
    max_depth      = 500
    bin_size       = 100
    mapq_cutoff    = 20
    local_assembly_out_file = os.path.join(out_dir, 'local_assembly.del.bedpe')
    short_reads_del_call_file = os.path.join(out_dir, 'discordant_read_pairs.del.bedpe')

    local_assembly.small_deletion_dection_by_local_assembly(samtools, bedtools, fermikit_dir, input_bam_file, ref_fasta_file, faidx_file, out_dir, local_assembly_out_file, n_threads, window_size, max_depth, rm_temp_files)
    

    cluster_weird_reads.cluster_weird_reads(in_weird_reads_file, weird_reads_cluster_file, faidx_file)

    cmd = '%s %s %s %s %d %d' %  (cal_hap_read_depth_from_bcd21, bcd21_file, hap_type_read_depth_file, faidx_file, bin_size, mapq_cutoff )
    my_utils.myprint(cmd)
    os.system(cmd)

    
    cmd = '%s %s %s %s %s %s %s' % (call_small_deletions, hap_type_read_depth_file, weird_reads_cluster_file, bcd22_file, faidx_file, gap_region_bed_file, short_reads_del_call_file)
    my_utils.myprint(cmd)
    os.system(cmd)
    


    tid2chrname_list, chrname2tid_dict = my_utils.get_chrnames(faidx_file)
    merge_sv_calls(local_assembly_out_file, short_reads_del_call_file, out_del_call_file, tid2chrname_list, chrname2tid_dict)

    if rm_temp_files:
        os.remove(local_assembly_out_file)
        os.remove(short_reads_del_call_file)
 

    return

def merge_sv_calls(local_assembly_out_file, short_reads_del_call_file, out_del_call_file, tid2chrname_list, chrname2tid_dict):
    
    local_assembly_del_call_list = svtk.read_sv_bedpe_file(local_assembly_out_file, chrname2tid_dict)
    local_assembly_del_call_list = svtk.remove_redundantsv(local_assembly_del_call_list)
    short_reads_del_call_list    = svtk.read_sv_bedpe_file(short_reads_del_call_file, chrname2tid_dict)
    merged_del_call_list         = svtk.merge2svcallset(local_assembly_del_call_list, short_reads_del_call_list)
    merged_del_call_list.sort(key = lambda svcall: svcall.key1())
    svtk.output_sv_list(merged_del_call_list, out_del_call_file)

    return


    
    
    
if __name__ == '__main__':
    main()




