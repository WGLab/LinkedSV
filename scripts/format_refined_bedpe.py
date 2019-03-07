#!/usr/bin/env python

import os
import sys
from bedpe import *

tab  = '\t'
endl = '\n'
arg = sys.argv[1:]

usage = 'python ' + __file__ + ' ' + '<refined.bedpe> <formatted_bedpe>'
argc  = 2 

def main():
    if len(arg) < argc:
        print usage
        sys.exit()

    in_bedpe_file = arg.pop(0)
    out_bedpe_file = arg.pop(0)

    in_bedpe_fp = open(in_bedpe_file, 'r')
    out_bedpe_fp = open( out_bedpe_file, 'w')

    while 1:
        line = in_bedpe_fp.readline()
        if not line: break
        line = line.strip().split(tab)
        refined_bedpe = RefinedQuantBKCand(line) 
        outstring = '%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%s\t%s\t%f\t%f\t%f\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%s\t%s\t%d\t%d\t%d\t%d' % (refined_bedpe.chrm1, refined_bedpe.refined_bk1, refined_bedpe.refined_bk1+1, refined_bedpe.chrm2, refined_bedpe.refined_bk2, refined_bedpe.refined_bk2+1, refined_bedpe.svtype, refined_bedpe.svlength, refined_bedpe.num_fragment_support, refined_bedpe.endtype1, refined_bedpe.endtype2, refined_bedpe.score, refined_bedpe.type_score, refined_bedpe.refine_score, refined_bedpe.n_pe_support, refined_bedpe.n_sr_support, refined_bedpe.n_support_reads, refined_bedpe.logp_nosv_one_mol, refined_bedpe.logp_nosv_two_mol, refined_bedpe.logp_sv_one_mol, refined_bedpe.logp_sv_two_mol, refined_bedpe.endtype1_logp, refined_bedpe.endtype2_logp, refined_bedpe.start1_logp, refined_bedpe.end1_logp, refined_bedpe.start2_logp, refined_bedpe.end2_logp, refined_bedpe.support_frm_ids1, refined_bedpe.support_frm_ids2, refined_bedpe.support_barcodes, refined_bedpe.tid1, refined_bedpe.tid2, refined_bedpe.start1, refined_bedpe.start2) 
        
        out_bedpe_fp.write(outstring + endl)

    in_bedpe_fp.close()
    out_bedpe_fp.close()


if __name__ == '__main__':
    main()
