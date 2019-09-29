#!/usr/bin/env python

import numpy as np
import math

import bisect
from sklearn import  linear_model
from sklearn.metrics import r2_score
from scipy.stats import norm

try:
    from scripts.my_utils import *
except ImportError:
    from my_utils import *


tab  = '\t'
endl = '\n'
arg = sys.argv[1:]

usage = 'python ' + __file__ + ' ' + '<in.bcd12> <output_file (in.bcd13)> <is_wgs (1 or 0)>'
argc  = 3 

def main():

    #args, dbo_args, endpoint_args = parse_user_arguments()
    #cal_expected_overlap_bcd_cnt(args, dbo_args, endpoint_args)

    if len(arg) < argc:
        print (usage)
        sys.exit()

    bcd12_file = arg.pop(0)
    out_file   = arg.pop(0)
    is_wgs = int(arg.pop(0))

    cal_expected_overlap_bcd_cnt(bcd12_file, out_file, is_wgs)

    return

def cal_expected_overlap_bcd_cnt(bcd12_file, out_file, is_wgs):

    if is_wgs != 1 and is_wgs != 0:
        myprint('ERROR! is_wgs should be 1 or 0')
        sys.exit()


    m1_all_list, m2_all_list, m0_all_list, d_all_list, n_lines_bcd12 = read_bcd12_file (bcd12_file, 0, 0, 10000000, 10000000)

    bk_list = range(0, 101)
    m1_mean, m1_std_r, m1_std_l = fit_mean_std(m1_all_list) 
    m2_mean, m2_std_r, m2_std_l = fit_mean_std(m2_all_list) 

    myprint('m1_mean, m1_std_r, m1_std_l: %f, %f, %f' % (m1_mean, m1_std_r, m1_std_l))
    myprint('m2_mean, m2_std_r, m2_std_l: %f, %f, %f' % (m2_mean, m2_std_r, m2_std_l))

    if is_wgs:
        min_m1_value = max(m1_mean - 3 * m1_std_l, m1_mean / 3.0)
        min_m2_value = max(m2_mean - 3 * m2_std_l, m2_mean / 3.0)
        max_m1_value = min(m1_mean + 3 * m1_std_r, m1_mean * 3.0)
        max_m2_value = min(m2_mean + 3 * m2_std_r, m2_mean * 3.0) 
    else:
        min_m1_value = m1_mean / 3.0 
        min_m2_value = m2_mean / 3.0
        max_m1_value = m1_mean * 3.0
        max_m2_value = m2_mean * 3.0

    myprint('min_m1_value, min_m2_value, max_m1_value, max_m2_value: %.f, %.f, %.f, %.f' % (min_m1_value, min_m2_value, max_m1_value, max_m2_value))


    del m1_all_list, m2_all_list, m0_all_list, d_all_list
    m1_all_list, m2_all_list, m0_all_list, d_all_list, n_lines_bcd12 = read_bcd12_file (bcd12_file, min_m1_value, min_m2_value, max_m1_value, max_m2_value)

    myprint('number of windows passed filtering: %d (%.2f %%)' % ( len(m1_all_list), len(m1_all_list)/float(n_lines_bcd12) * 100.0 ))

    logm0_array = np.array(np.log(m0_all_list))
    logm1_array = np.log(m1_all_list)
    logm2_array = np.log(m2_all_list)
    d_array     = np.array(d_all_list)

    myprint ('fitting model')

    x = np.array([logm1_array, logm2_array, d_array], np.float64).transpose()
    y = logm0_array

    regr = linear_model.LinearRegression()
    regr.fit(x, y)
    y_predict = regr.predict(x)
    m0_predict = np.exp(y_predict)

    myprint('finished fitting model')

    a = regr.coef_[0]
    b = regr.coef_[1]
    alpha = regr.coef_[2]
    logn = 0 - regr.intercept_ 

    Y_list = list()
    for i in range(0, len(m1_all_list)):
        m0 = m0_all_list[i]
        m1 = m1_all_list[i]
        m2 = m2_all_list[i]
        d = d_all_list[i]
        predicted_m0 = math.exp(a * math.log(m1) + b * math.log(m2) + alpha * d - logn)
        if m0 == 0: m0 = 0.2
        Y = math.log(predicted_m0 / float(m0), 2)
        Y_list.append(Y)

    Y_mean, Y_std_r, Y_std_l = fit_mean_std(Y_list)
    myprint('Y_mean, Y_std_r, Y_std_l: %f, %f, %f' % (Y_mean, Y_std_r, Y_std_l))
    n_window_pairs = len(Y_list) 
    bk_list = list()
    c = 1.05
    max_n = int(math.log(len(Y_list), c))
    for i in range(1, max_n):
        bk = 100.0 - 100.0 / math.pow(c,i)
        bk_list.append(bk)

    Y_percentile_array = np.percentile(Y_list, bk_list)

    del m1_all_list, m2_all_list, m0_all_list, d_all_list
    

    predict_overlap_bcd_cnt(bcd12_file, out_file, a, b, alpha, logn, min_m1_value, min_m2_value, max_m1_value, max_m2_value, Y_percentile_array, bk_list, c, Y_mean, Y_std_r, Y_std_l, n_window_pairs)



    return 

def predict_overlap_bcd_cnt(bcd12_file, out_file, a, b, alpha, logn, min_m1_value, min_m2_value, max_m1_value, max_m2_value, Y_percentile_array, bk_list, c, Y_mean, Y_std_r, Y_std_l, n_window_pairs): 
    
    myprint ('calculating expected overlap barcode count')
    header = '#tid\twin2_start_pos\twin2_end_pos\tn_bcd_win1\tn_bcd_win2\tn_ovl_bcd\tpredicted_n_ovl_bcd\tminus_log10_empirical_p_value\tminus_log10_p_value'

    out_fp = open(out_file, 'w')
    out_fp.write(header + endl)

    bcd12_fp = open(bcd12_file, 'r')
    n_bcd12_line = 0
    n_predicted_line = 0

    while 1:
        line = bcd12_fp.readline()
        if not line: break
        if line[0] == '#': continue
        n_bcd12_line += 1

        if n_bcd12_line % 10000000 == 0: myprint('processed %d window pairs' % n_bcd12_line)

        line = line.strip().split(tab)

        m1 = float(line[3])
        m2 = float(line[4])
        m0 = float(line[5])

        centroid1 = float(line[6])
        centroid2 = float(line[7])

        if m1 < min_m1_value or m2 < min_m2_value or centroid1 < 0 or centroid2 < 0 or m1 > max_m1_value or m2 > max_m2_value:
            predicted_n_ovl_bcd = -1
            minus_log10_empirical_p_value = 0.0
            minus_log10_p_value = 0.0
            ret_string = tab.join(line[0:6]) + '\t%d\t%.4f\t%.4f\n' % (predicted_n_ovl_bcd, minus_log10_empirical_p_value, minus_log10_p_value)
            out_fp.write(ret_string)
            continue

        d = centroid2 - centroid1
        predicted_m0 = math.exp(a * math.log(m1) + b * math.log(m2) + alpha * d - logn)
        predicted_n_ovl_bcd = int(predicted_m0 + 0.5)
        if predicted_n_ovl_bcd > min(m1, m2): predicted_n_ovl_bcd = min(m1, m2)
        if m0 < 1: m0 = 0.2
        Y = math.log(predicted_n_ovl_bcd / float(m0), 2)
        empirical_p_value = cal_empirical_p_value(Y, Y_percentile_array, bk_list, c) 
        minus_log10_empirical_p_value = 0 - math.log(empirical_p_value, 10)

        zscore = (Y - Y_mean) / Y_std_r
        p_value = norm.sf(zscore, 0, 1) 

        if p_value < 1e-100: 
            p_value = 1e-100
            minus_log10_p_value = 100.0
        else:
            minus_log10_p_value = 0 - math.log(p_value, 10)

        ret_string = tab.join(line[0:6]) + '\t%d\t%.4f\t%.4f\n' % (predicted_n_ovl_bcd, minus_log10_empirical_p_value, minus_log10_p_value)
        out_fp.write(ret_string)
        n_predicted_line += 1

    bcd12_fp.close()

    out_fp.close()
    myprint ('finished writing to file: %s' % out_file)
    myprint ('%d (%.2f %%) windows have valid predictions' % (n_predicted_line, float(n_predicted_line) / float(n_bcd12_line) * 100.0))

    return

def cal_empirical_p_value(Y, Y_percentile_array, bk_list, c):

    idx = bisect.bisect(Y_percentile_array, Y) - 1
    if idx < 0: idx = 0
    if idx > len(Y_percentile_array) - 1: idx = len(Y_percentile_array) - 1
    p_value = (100.0 - bk_list[idx] ) / 100.0
    return p_value

def read_bcd12_file(bcd12_file, min_m1_value, min_m2_value, max_m1_value, max_m2_value):

    myprint ('reading bcd12 file: %s' % (bcd12_file))

    m1_all_list = list()
    m2_all_list = list()
    m0_all_list = list()
    d_all_list  = list()
    
    n_lines_bcd12 = 0

    bcd12_fp = gzopen(bcd12_file, 'r')
    while 1:
        line = bcd12_fp.readline()
        if not line: break
        if line[0] == '#': continue
        n_lines_bcd12 += 1
        line = line.strip().split(tab)

        m1 = float(line[3])
        m2 = float(line[4])
        m0 = float(line[5])

        if m1 < min_m1_value or m2 < min_m2_value or m0 < 1: continue
        if m1 > max_m1_value or m2 > max_m2_value: continue

        centroid1 = float(line[6])
        centroid2 = float(line[7])

        if centroid1 < 0 or centroid2 < 0: continue

        d = centroid2 - centroid1

        m1_all_list.append(m1)
        m2_all_list.append(m2)
        m0_all_list.append(m0)

        d_all_list.append(d)

    bcd12_fp.close()
    
    myprint ('finished reading bcd12 files')

    return m1_all_list, m2_all_list, m0_all_list, d_all_list, n_lines_bcd12

def fit_mean_std(data_list):

    if len(data_list) < 3:
        myprint('ERROR! Only %d elements in data_list!' % len(data_list) )
        sys.exit()

    mean = np.median(data_list)
    r_std = np.percentile(data_list, 84.135) - np.percentile(data_list, 50.0)
    l_std = np.percentile(data_list, 50.0) - np.percentile(data_list, 15.865)
    #std2 = (np.percentile(data_list, 97.725) - mean / 2.0
    #std = (std1 + std2) / 2.0
    return mean, r_std, l_std 









if __name__ == '__main__':
    main()
