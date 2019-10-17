#!/usr/bin/env python

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
try:
    from scripts.arguments import *
except ImportError:
    from arguments import *


tab  = '\t'
endl = '\n'
arg = sys.argv[1:]



def main():
    
    args, dbo_args, endpoint_args = parse_user_arguments()

    detect_sv_from_short_reads(args, dbo_args, endpoint_args)

    return 

def detect_sv_from_short_reads(args, dbo_args, endpoint_args):


    out_short_read_sv_file = args.short_reads_sv_call_file ## combination of local assembly, discordant read pairs and large CNV

    ### calculate read depth ###
    
    bin_size = 100
    mapq_cutoff = 20

    if os.path.exists(args.hap_type_read_depth_file) == False or os.path.getsize(args.hap_type_read_depth_file) == 0:
        
        cmd = '%s %s %s %s %d %d' %  (args.cal_hap_read_depth_from_bcd21, endpoint_args.bcd21_file, args.hap_type_read_depth_file, args.faidx_file, bin_size, mapq_cutoff )
        my_utils.myprint('running command:' + cmd)
        os.system(cmd)

    ### CNV  detection ### 
    cmd = '%s %s %s %s %s %d %d %d' % (args.cnv_detection, args.hap_type_read_depth_file, args.faidx_file, args.gap_region_bed_file, args.cnv_call_file, 40, 200, 500000)
    my_utils.myprint('running command:' + cmd)
    os.system(cmd)
    
    ### small deletion detection from paired-end reads ###

    window_size    = int(2e5)
    max_depth      = 500
    bin_size       = 100
    mapq_cutoff    = 20

    local_assembly_out_file = os.path.join(args.out_dir, 'local_assembly.del.bedpe')
    short_reads_del_call_file = os.path.join(args.out_dir, 'discordant_read_pairs.del.bedpe')
    
    cluster_weird_reads.cluster_weird_reads(args.weird_reads_file, args.weird_reads_cluster_file, args.faidx_file)

    cmd = '%s %s %s %s %s %s %s' % (args.small_deletion_detection, args.hap_type_read_depth_file, args.weird_reads_cluster_file, endpoint_args.bcd22_file, args.faidx_file, args.gap_region_bed_file, short_reads_del_call_file)
    my_utils.myprint('running command:' + cmd)
    os.system(cmd)

    rm_temp_files = 1
    local_assembly.small_deletion_dection_by_local_assembly(args.samtools, args.bedtools, args.fermikit_dir, args.input_bam, args.ref_fa, args.faidx_file, args.out_dir, local_assembly_out_file, args.n_thread, window_size, max_depth, rm_temp_files)

    ### merge call files ###

    merge_sv_calls(local_assembly_out_file, short_reads_del_call_file, out_short_read_sv_file, args.tid2chrname, args.chrname2tid)
    
    rm_temp_files = 1 
    if rm_temp_files: 
        os.remove(local_assembly_out_file)
        os.remove(short_reads_del_call_file)

    return

def merge_sv_calls(local_assembly_out_file, short_reads_del_call_file, out_short_read_sv_file, tid2chrname_list, chrname2tid_dict):
    
    local_assembly_del_call_list = svtk.read_sv_bedpe_file(local_assembly_out_file, chrname2tid_dict)
    local_assembly_del_call_list = svtk.remove_redundantsv(local_assembly_del_call_list)
    short_reads_del_call_list    = svtk.read_sv_bedpe_file(short_reads_del_call_file, chrname2tid_dict)
    merged_del_call_list         = svtk.merge2svcallset(local_assembly_del_call_list, short_reads_del_call_list)
    merged_del_call_list.sort(key = lambda svcall: svcall.key1())
    svtk.output_sv_list(merged_del_call_list, out_short_read_sv_file)
    
    return


    
    
    
if __name__ == '__main__':
    main()




