#!/usr/bin/env python

import argparse
import copy
import os
import sys
import my_utils

class global_parameter:
    def __init__(self, parser_args):

        self.root_dir = os.path.join(os.path.split(os.path.abspath(__file__))[0], '../')

        self.input_bam  = os.path.abspath(parser_args.bam)
        self.out_dir    = os.path.abspath(parser_args.out_dir)
        self.ref_fa     = os.path.abspath(parser_args.ref) 

        self.ref_version = parser_args.ref_version
        self.gap_region_bed_file   = parser_args.gap_region_bed
        self.black_region_bed_file = parser_args.black_region_bed

        self.black_region_2d_file  = ''
        self.low_mapq_region_bed_file = ''

        self.min_mapq = 20
        self.n_thread = parser_args.n_thread

        self.samtools = parser_args.samtools
        self.bedtools = parser_args.bedtools
        self.is_wgs             = parser_args.is_wgs
        self.target_region_bed  = parser_args.target_region
        if parser_args.save_temp_files == True:
            self.rm_temp_files      = False
        else:
            self.rm_temp_files      = True

        self.germline_mode      = parser_args.germline_mode

        self.temp_file_list = list()

        self.run_from_begining = True
        self.only_method1      = False
        self.only_method2      = True

        self.bam_name = os.path.split(self.input_bam)[1]
        self.out_prefix = os.path.join(self.out_dir, self.bam_name)
        self.bam = self.input_bam

        self.sortn_bam = self.out_prefix + '.sortn.bam' 
        self.sortbx_bam = self.out_prefix + '.sortbx.bam' 

        self.sortn_bam_core_file = self.out_prefix + '.sortn.bam.coreinfo.gz' 
        self.weird_reads_file  = self.out_prefix + '.weird_reads.txt' 
        
        self.alt_chr_name_set = set()
        self.alt_tid_set = set()

        self.faidx_file = self.ref_fa + '.fai'

        self.extract_barcode = os.path.join(self.root_dir, 'bin/extract_barcode_info')
        self.remove_sparse_nodes = os.path.join(self.root_dir, 'bin/remove_sparse_nodes');
        self.output_bam_coreinfo = os.path.join(self.root_dir, 'bin/output_bam_coreinfo');
        self.cal_read_depth_from_bcd21 = os.path.join(self.root_dir, 'bin/cal_read_depth_from_bcd21');
        self.cal_barcode_depth_from_bcd21 = os.path.join(self.root_dir, 'bin/cal_barcode_depth_from_bcd21');
        self.cal_centroid_from_read_depth = os.path.join(self.root_dir, 'bin/cal_centroid_from_read_depth');
        self.cal_twin_win_bcd_cnt = os.path.join(self.root_dir, 'bin/cal_twin_win_bcd_cnt');
        self.grid_overlap = os.path.join(self.root_dir, 'bin/grid_overlap');
        self.split_weird_reads_program = os.path.join(self.root_dir, 'scripts/split_weird_reads_file.py')
        self.cal_expected_overlap_value = os.path.join(self.root_dir, 'scripts/cal_expected_overlap_value.py')
        self.cal_2d_overlapping_barcodes = os.path.join(self.root_dir, 'bin/cal_2d_overlapping_barcodes')
        self.pigz = os.path.join(self.root_dir, 'bin/pigz')
        self.cluster_reads = os.path.join(self.root_dir, 'bin/cluster_reads')

        self.alt_ctg_file  = os.path.join(self.root_dir, 'black_lists/alternative_contigs.txt')

        self.user_defined_min_frag_length = parser_args.min_fragment_length
        self.user_defined_gap_distance_cut_off = parser_args.gap_distance_cut_off
        self.user_defined_min_reads_in_fragment = parser_args.min_reads_in_fragment

        self.bcd_file  = self.out_prefix + '.bcd' 
        self.args_file = self.out_prefix + '.arguments' 
        self.num_split_bcd21_file = self.out_prefix  + '.num_split_bcd21_file'

        self.stat_file = self.out_prefix + '.barcode_statistics' 

        self.genome_length = None
        self.target_region_length = None
        
        self.global_distribution_file = self.out_prefix + '.fragment_statistics'

        self.node33_file = self.out_prefix + '.node33'
        self.node55_file = self.out_prefix + '.node55'
        self.node53_file = self.out_prefix + '.node53'
        self.node35_file = self.out_prefix + '.node35'
        self.node33_candidate_file = self.out_prefix + '.node33.candidates'
        self.node55_candidate_file = self.out_prefix + '.node55.candidates'
        self.node53_candidate_file = self.out_prefix + '.node53.candidates'
        self.node35_candidate_file = self.out_prefix + '.node35.candidates'
        self.n_node33 = None

        self.node_cluster33_file = self.out_prefix + '.node_cluster33'
        self.node_cluster55_file = self.out_prefix + '.node_cluster55'
        self.node_cluster53_file = self.out_prefix + '.node_cluster53'
        self.node_cluster35_file = self.out_prefix + '.node_cluster35'

        self.bk_cand_pair_file = self.out_prefix + '.bk_cand_pairs'
        self.quantified_bk_pair_file = self.out_prefix + '.qbkpair.bedpe'
        self.refinedbedpe_file = self.out_prefix + '.qbkpair.refined.bedpe'
        self.merged_bedpe_file = self.out_prefix + '.raw_svcalls.bedpe'
        self.filter_bedpe_file = self.out_prefix + '.filtered_svcalls.bedpe'

        self.read_depth_file = self.out_prefix + '.read_depth.txt'
        self.hap_type_read_depth_file = self.out_prefix + '.hap_depth.txt'
        self.chrname2tid = None
        self.tid2chrname = None
        self.global_distribution_calculated = False

        ### statistics ###
        self.length_empirical_pmf = None
        self.length_empirical_cdf = None
        self.fragment_length_lmda = None
        self.read_per_bp_genome = None
        self.read_per_bp_ontarget = None
        self.read_per_bp_offtarget = None
        self.num_reads_genome = None
        self.num_reads_ontarget = None
        self.num_reads_offtarget = None
        self.mean_num_fragment_per_bcd = None
        self.median_num_fragment_per_bcd = None
        self.total_num_fragment = None
        self.total_num_bcd = None

        self.median_isize = 300 
        self.max_isize_cutoff = 500

        self.rd1_length = 85
        self.rd2_length = 100
        
        self.gap_distance_pmf = None
        self.gap_distance_cdf = None
        self.gap_distance750 = None
        self.gap_distance_lmda = None

        self.gap_distance_cutoff = 10000
        self.gap_distance500 = 200
        self.gap_distance900 = 3000
        self.gap_distance950 = 5000
        self.gap_distance990 = 7000
        self.gap_distance999 = 10000 
        self.median_fragment_length = 50000
        self.mean_fragment_length = 72000   

        if parser_args.min_supp_barcodes > 5:
            self.min_support_fragments = parser_args.min_supp_barcodes
        elif parser_args.min_supp_barcodes > 0 and parser_args.min_supp_barcodes <= 5:
            self.min_support_fragments = 5 
        else:
            self.min_support_fragments = 10

    def copy(self): 
        return copy.deepcopy(self)
    

class dbo_parameter:
    def __init__(self, global_args):

        ### output files ###
        self.prefix  = global_args.out_prefix
        
        self.bcd11_file = self.prefix + '.bcd11' # twin window file 
        self.bcd12_file = self.prefix + '.bcd12' # twin window file plus centroid
        self.bcd13_file = self.prefix + '.bcd13' # twin window file plus predicted num of barcodes
        
        ### parameters for analysis ###
        self.bin_size = 100 
        self.min_bcd_num = 50 



class endpoint_parameter:
    def __init__(self, global_args):

        ### output files ###
        self.prefix = global_args.out_prefix + '.endpoints'
        self.bcd21_file = global_args.out_prefix + '.bcd21.gz' # sorting results 
        self.low_mapq_bcd21_file = global_args.out_prefix + '.low_mapq.bcd21.gz' # sorting results 
        self.tmpbcd22_file = global_args.out_prefix + '.tmpbcd22'
        self.bcd21_file_of_target_region = global_args.out_prefix + '.on_target.bcd21.gz'
        self.bcd22_file = global_args.out_prefix + '.bcd22' # result of clustering (fragment file)
        self.bcd23_file = global_args.out_prefix + '.bcd23' # candidate fragments (long) 
        self.bcd24_file = global_args.out_prefix + '.bcd24' # candidate fragments (short) 
        self.barcode_cov_file = global_args.out_prefix + '.barcode_cov.bed'
        self.high_cov_file = global_args.out_prefix + '.high_cov.bed'
        self.bk_file = self.prefix + '.bk.bed'
        self.endpoints_density_file = global_args.out_prefix + '.endpoints_density' 

        self.peak_file = self.prefix + '.bk.peak'
        self.peak2d_file = self.prefix + '.bk.peak2d'
        self.cluster_file = self.prefix + '.cluster'
        self.sv_cand_file = self.prefix + '.sv_cand.bedpe'
        self.one_frm_sv_cand_file = self.prefix + '.one_frm_sv_cand.bedpe'

        ### parameters for analysis ###
        #self.win_size = 1000
        self.min_distance_for_peak_calling = 100

        self.bin_len_for_calculate_endpoint_density = 100

        self.bin_len_for_dbscan = 1000
        self.min_frag_length = 5000 # will be updated when running

        #self.min_shared_frag = 3
        self.min_support_barcode_cnt = 3
        self.min_interval_length = 12000 # for dbscan

    def copy(self):
        return copy.deepcopy(self)

def check_arguments(args):

    if os.path.exists(args.input_bam) == False:
        my_utils.myprint("ERROR! input bam file (%s) does not exist!" %(args.bam))
        sys.exit()

    if os.path.exists(args.out_dir) == False:
        os.system('mkdir -p ' + args.out_dir)
        if os.path.exists(args.out_dir) == False:
            my_utils.myprint("ERROR! can not create output directory: %s" %(args.out_dir))
            sys.exit()

    if os.path.exists(args.ref_fa) == False:
        my_utils.myprint("ERROR! reference FASTA file (%s) does not exist!" % (args.ref_fa))
        sys.exit()

    if os.path.exists(args.faidx_file) == False:
        my_utils.myprint("Index file of reference FASTA file does not exist, creating index using samtools ...")
        cmd = args.samtools + ' faidx ' + args.ref_fa
        os.system(cmd)
        if os.path.exists(args.faidx_file) == False:
            my_utils.myprint("ERROR! Cannot generate reference index file!")
            sys.exit()

    if os.path.exists(args.extract_barcode) == False:
        my_utils.myprint("ERROR! extract_barcode does not exist!")
        sys.exit()

    if args.is_wgs == False and os.path.exists(args.target_region_bed) == False:
        my_utils.myprint ("ERROR! target region bed is required if --targeted is specified. If you don't have this file, please specify --wgs instead")
        sys.exit()


    if args.ref_version == 'hg19':
        args.gap_region_bed_file      = os.path.join(args.root_dir, 'black_lists/hg19_gap.bed') 
        args.black_region_bed_file    = os.path.join(args.root_dir, 'black_lists/hg19_black_list.bed')
        args.black_region_2d_file     = os.path.join(args.root_dir, 'black_lists/hg19.2D.blacklist.gz')
        args.low_mapq_region_bed_file = os.path.join(args.root_dir, 'black_lists/hg19_low_mapq_regions.bed')

    elif args.ref_version == 'b37':
        args.gap_region_bed_file      = os.path.join(args.root_dir, 'black_lists/b37_gap.bed') 
        args.black_region_bed_file    = os.path.join(args.root_dir, 'black_lists/b37_black_list.bed')
        args.black_region_2d_file     = os.path.join(args.root_dir, 'black_lists/b37.2D.blacklist.gz')
        args.low_mapq_region_bed_file = os.path.join(args.root_dir, 'black_lists/b37_low_mapq_regions.bed')

    elif args.ref_version == 'hg38':
        args.gap_region_bed_file      = os.path.join(args.root_dir, 'black_lists/hg38_gap.bed') 
        args.black_region_bed_file    = os.path.join(args.root_dir, 'black_lists/hg38_black_list.bed')
        args.black_region_2d_file     = os.path.join(args.root_dir, 'black_lists/hg38.2D.blacklist.gz')
        args.low_mapq_region_bed_file = os.path.join(args.root_dir, 'black_lists/hg38_low_mapq_regions.bed')



    return

def parse_user_arguments():

    parser = argparse.ArgumentParser(description='Detection of SVs from linked-read sequencing data')
    ### required arguments ###
    parser.add_argument('-i', '--bam', required = True, metavar = 'input.phased_possorted_bam.bam', type = str, help = 'input bam file (should be the phased_possorted_bam.bam generated by Longranger')
    parser.add_argument('-d', '--out_dir', required = True, metavar = 'output_directory', type = str, help = 'output directory')
    parser.add_argument('-r', '--ref', required = True, metavar = 'ref.fa', type = str, help ='reference FASTA file')

    ### optional arguments ###
    parser.add_argument('-v', '--ref_version', required = False, metavar = 'version', type = str, default = '', help ='version of reference fasta file. Current supported versions are: hg19, b37, hg38')
    parser.add_argument('--gap_region_bed', required = False, metavar = 'BED', type = str, default = '', help ='reference gap region in bed format, required if --ref_version is not specified')
    parser.add_argument('--black_region_bed', required = False, metavar = 'BED', type = str, default = '', help ='black region in bed format, required if --ref_version is not specified')
    parser.add_argument('-t', '--n_thread', required = False, metavar = 'num_thread', type = int, default = 1, help ='number of threads (default: 4)')
    parser.add_argument('--min_fragment_length', metavar = 'INT', required = False, type = int, default = -1, help ='minimal fragment length considered for SV calling')
    parser.add_argument('--min_reads_in_fragment', metavar = 'INT', required = False, type = int, default = -1, help ='minimal number of confidently mapped reads in one fragment')
    parser.add_argument('--min_supp_barcodes', metavar = 'INT', required = False, type = int, default = 10, help ='minimal number of shared barcodes between two SV breakpoints (default: 10)')
    parser.add_argument('--samtools', required = False, metavar = 'path/to/samtools', type = str, default = 'samtools', help ='path to samtools (default: find in environmental path)')
    parser.add_argument('--bedtools', required = False, metavar = 'path/to/bedtools', type = str, default = 'bedtools', help ='path to bedtools (default: find in environmental path)')
    parser.add_argument('--wgs', dest='is_wgs', action='store_true', help='the input is whole-genome sequencing data')
    parser.add_argument('--targeted', dest='is_wgs', action='store_false', help='the input is targeted region sequencing data (such as WES)')
    parser.add_argument('--germline_mode', dest='germline_mode', action='store_true', help='detect germline SVs')
    parser.add_argument('--somatic_mode', dest='germline_mode', action='store_false', help='detect somatic SVs (with low variant allele frequency)')
    parser.set_defaults(germline_mode = True)
    parser.set_defaults(is_wgs = True)
    parser.add_argument('--target_region', required = False, metavar = 'BED', type = str, default = '', help ='bed file of target regions (required if --targeted is specified)')
    parser.add_argument('--gap_distance_cut_off', required = False, metavar = 'INT', type = int, default = -1, help ='max distance between two reads in a HMW DNA molecule (default: automatically determined)')
    parser.add_argument('--save_temp_files', dest='save_temp_files', action='store_true', help='Do not remove intermediate files after the run. Use in debug mode. Default: False')

    input_args = parser.parse_args()

    args = global_parameter(input_args)

    args.tid2chrname, args.chrname2tid = my_utils.get_chrnames(args.faidx_file)

    args.alt_chr_name_set = my_utils.read_alternative_contig_file(args.alt_ctg_file)

    args.alt_tid_set = my_utils.get_alternative_tid_set(args.alt_ctg_file, args.faidx_file)

    check_arguments(args)

    dbo_args = dbo_parameter(args)

    endpoint_args = endpoint_parameter(args)

    return args, dbo_args, endpoint_args

def output_arguments2file(args, dbo_args, endpoint_args):
    out_fp = open(args.args_file, 'w')
    for key in args.__dict__:    
        value = args.__dict__[key]
        if type(value) is not list and type(value) is not dict and type(value) is not set:
            print >> out_fp, key, '\t', value

    for key in dbo_args.__dict__:    
        value = dbo_args.__dict__[key]
        if type(value) is not list and type(value) is not dict and type(value) is not set:
            print >> out_fp, key, '\t', value

    for key in endpoint_args.__dict__:    
        value = endpoint_args.__dict__[key]
        if type(value) is not list and type(value) is not dict and type(value) is not set:
            print >> out_fp, key, '\t', value

    out_fp.close()
    return
