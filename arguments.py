#!/usr/bin/env python

import argparse
import copy
import os
import sys
import my_utils

class global_parameter:
    def __init__(self, parser_args):

        self.input_bam  = os.path.abspath(parser_args.bam)
        self.out_dir    = os.path.abspath(parser_args.out_dir)
        self.ref_fa     = os.path.abspath(parser_args.ref) 

        self.min_mapq = parser_args.min_mapq
        self.n_thread = parser_args.n_thread

        self.samtools = parser_args.samtools
        self.bedtools = parser_args.bedtools

        self.is_wgs             = parser_args.is_wgs
        self.all_to_all         = parser_args.all_to_all
        self.target_region_bed  = parser_args.target_region

        self.debug             = parser_args.debug
        self.run_from_begining = parser_args.run_from_begining
        self.only_method1      = parser_args.only_method1
        self.only_method2      = parser_args.only_method2

        self.bam_name = os.path.split(self.input_bam)[1]
        self.out_prefix = os.path.join(self.out_dir, self.bam_name)
        self.bam = self.input_bam

        self.faidx_file = self.ref_fa + '.fai'

        self.root_dir = os.path.split(os.path.abspath(__file__))[0]
        self.extract_barcode = os.path.join(self.root_dir, 'extract_barcode')
        self.sort_barcode = os.path.join(self.root_dir, 'sort_barcode')

        self.user_defined_min_frag_length = parser_args.min_fragment_length

        self.bcd_file = self.out_prefix + '.bcd' 
        self.args_file = self.out_prefix + '.arguments' 
        self.bcd_file_of_target_region = self.out_prefix + '.on_target.bcd'

        self.stat_file = self.out_prefix + '.barcode_statistics' 

        self.min_sv_length = 20000

        self.genome_length = None
        self.target_region_length = None
        
        self.global_distribution_file = self.out_prefix + '.fragment_statistics'

        self.node33_file = self.out_prefix + '.node33'
        self.node55_file = self.out_prefix + '.node55'
        self.node53_file = self.out_prefix + '.node53'
        self.node35_file = self.out_prefix + '.node35'
        self.n_node33 = None

        self.node_cluster33_file = self.out_prefix + '.node_cluster33'
        self.node_cluster55_file = self.out_prefix + '.node_cluster55'
        self.node_cluster53_file = self.out_prefix + '.node_cluster53'
        self.node_cluster35_file = self.out_prefix + '.node_cluster35'

        self.bk_cand_pair_file = self.out_prefix + '.bk_cand_pairs'
        self.quantified_bk_pair_file = self.out_prefix + '.qbkpair.bedpe'
        self.refinedbedpe_file = self.out_prefix + '.qbkpair.refined.bedpe'
        self.merged_bedpe_file = self.out_prefix + '.svcalls.bedpe'

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

        self.min_support_fragments = 3

    def copy(self): 
        return copy.deepcopy(self)
    

class dbo_parameter:
    def __init__(self, global_args):

        ### output files ###
        self.prefix  = global_args.out_prefix + '.dbo'
        self.bk_file = self.prefix + '.bk.bed'
        self.peak_file = self.prefix + '.bk.peak' 
        self.peak2d_file = self.prefix + '.bk.peak2d'
        #self.nonadjcent_empirical_distribution_file = global_args.out_prefix + '.nonadjacent_distribution.txt'
        #self.twin_empirical_dist_file = global_args.out_prefix +'.twin_empirical_dist.txt'
        
        self.chr_bcd_file_list = list()
        self.bcd11_file_list = list()
        self.bcd12_file_list = list()
        self.bcd13_file_list = list()
        for tid in range(0, len(global_args.tid2chrname)):
            chr_bcd_file = global_args.out_prefix + '.tid_' + str(tid) + '.bcd'
            bcd11_file   = global_args.out_prefix + '.tid_' + str(tid) + '.bcd11'
            bcd12_file   = global_args.out_prefix + '.tid_' + str(tid) + '.bcd12'
            bcd13_file   = global_args.out_prefix + '.tid_' + str(tid) + '.bcd13'
            self.chr_bcd_file_list.append(chr_bcd_file)
            self.bcd11_file_list.append(bcd11_file)
            self.bcd12_file_list.append(bcd12_file)
            self.bcd13_file_list.append(bcd13_file)
        
        ### parameters for analysis ###
        self.bin_size = 100 
        if global_args.is_wgs:
            self.win_size = 100 
        else:
            self.win_size = 500 

        #self.num_random_bk_pairs = 20000 ### number of random bk pairs for calculating nonadjcent_empirical_distribution
        self.min_bcd_num = 50 
        self.min_dist_random_bk = 20000 ### min distance between random bk pairs
        self.max_dist_random_bk = 100000 ### min distance between random bk pairs
        self.min_distance_for_peak_calling = 10

        self.min_shared_frag = 5 ### min shared fragments between two clusters to be considered as a same cluster group
        self.min_support_barcode_cnt = 10 # min support barcode count for two clusters to be considered as a paired breakpoint loci

    def copy(self):
        return copy.deepcopy(self)

class endpoint_parameter:
    def __init__(self, global_args):

        ### output files ###
        self.prefix = global_args.out_prefix + '.endpoints'
        self.bcd21_file = global_args.out_prefix + '.bcd21' # sorting results 
        self.tmpbcd22_file = global_args.out_prefix + '.tmpbcd22'
        self.bcd22_file = global_args.out_prefix + '.bcd22' # result of clustering (fragment file)
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

    if os.path.exists(args.sort_barcode) == False:
        my_utils.myprint ("ERROR! sort_barcode does not exist!")
        sys.exit()

    if args.is_wgs == False and os.path.exists(args.target_region_bed) == False:
        my_utils.myprint ("ERROR! target region bed is required if --targeted is specified. If you don't have this file, please specify --wgs instead")
        sys.exit()


    return

def parse_user_arguments():

    parser = argparse.ArgumentParser(description='Detection of SVs from linked-read sequencing data')
    ### required arguments ###
    parser.add_argument('-i', '--bam', required = True, metavar = 'input.sorted.bam', type = str, help = 'sorted bam file as input')
    parser.add_argument('-d', '--out_dir', required = True, metavar = 'out_directory', type = str, help = 'output directory')
    parser.add_argument('-r', '--ref', required = True, metavar = 'ref.fasta', type = str, help ='reference FASTA file')

    ### optional arguments ###
    parser.add_argument('-q', '--min_mapq', required = False, metavar = 'min_map_qual', type = int, default = 20, help ='minimal map quality of reads used for analysis (default: 20)')
    parser.add_argument('-t', '--n_thread', required = False, metavar = 'num_thread', type = int, default = 1, help ='number of threads (default: 1)')

    parser.add_argument('--min_fragment_length', metavar = 'min_fragment_length', required = False, type = int, default = -1, help ='minimal fragment length considered for SV calling')

    parser.add_argument('--samtools', required = False, metavar = 'path/to/samtools', type = str, default = 'samtools', help ='path to samtools (default: find in environmental path)')
    parser.add_argument('--bedtools', required = False, metavar = 'path/to/bedtools', type = str, default = 'bedtools', help ='path to bedtools (default: find in environmental path)')

    parser.add_argument('--wgs', dest='is_wgs', action='store_true', help='the input is whole-genome sequencing data')
    parser.add_argument('--targeted', dest='is_wgs', action='store_false', help='the input is targeted region sequencing data (such as WES)')
    parser.set_defaults(is_wgs = True)
    parser.add_argument('--target_region', required = False, metavar = 'target_region.bed', type = str, default = '', help ='bed file of target regions (required if --targeted is specified)')

    ## debug arguments ###
    parser.add_argument('--debug', dest='debug', action='store_true', help='run in debug mode')
    parser.set_defaults(debug = False)
    parser.add_argument('--run_from_begining', dest = 'run_from_begining', action='store_true', default=False, help ='run from begining (default: false)')
    parser.add_argument('--only_method1', dest = 'only_method1', action='store_true', default = False, help ='only run method 1 (default: False)')
    parser.add_argument('--only_method2', dest = 'only_method2', action='store_true', default = False, help ='only run method 2 (default: False)')
    parser.add_argument('--all_to_all', dest = 'all_to_all', action='store_true', default = False, help ='compare all genomic loci')

    input_args = parser.parse_args()

    args = global_parameter(input_args)
    args.tid2chrname, args.chrname2tid = my_utils.get_chrnames(args.faidx_file)
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
