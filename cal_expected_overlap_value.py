#!/usr/bin/env python

import numpy as np

from my_utils import *
from sklearn import  linear_model
from sklearn.metrics import r2_score

def main():

    args, dbo_args, endpoint_args = parse_user_arguments()
    cal_expected_overlap_bcd_cnt(args, dbo_args, endpoint_args)
    return

def read_all_bcd12_files(bcd12_file_list):

    gc.enable()
    myprint ('reading bcd12 files')
    m1_all_list = list()
    m2_all_list = list()
    m0_all_list = list()
    d_all_list = list()

    for bcd12_file in bcd12_file_list:
        myprint ('current bcd12 file: %s' % bcd12_file)
        bcd12_fp = open(bcd12_file, 'r')
        while 1:
            line = bcd12_fp.readline()
            if not line:
                break
            line = line.strip().split(tab)
            m1_all_list.append(float(line[1]))
            m2_all_list.append(float(line[2]))
            m0_all_list.append(float(line[3]))
            d_all_list.append(float(line[4]) + float(line[5]))

        bcd12_fp.close()
    
    myprint ('finished reading bcd12 files')
    return m1_all_list, m2_all_list, m0_all_list, d_all_list

def cal_expected_overlap_bcd_cnt(args, dbo_args, endpoint_args):

    is_wgs = args.is_wgs
    bcd12_file_list = dbo_args.bcd12_file_list
    bcd13_file_list = dbo_args.bcd13_file_list
    
    m1_all_list, m2_all_list, m0_all_list, d_all_list = read_all_bcd12_files(bcd12_file_list)
    min_m1_value = 50
    min_m2_value = 50
    myprint ('filtering windows... min_m1_value=%.2f, min_m2_value=%.2f' %(min_m1_value, min_m2_value))
    m1list = list()
    m2list = list()
    m0list = list()
    dlist = list()

    for i in range(0, len(m1_all_list)): 
        if m1_all_list[i] < min_m1_value or m2_all_list[i] < min_m1_value or m0_all_list[i] < 1:
            continue

        m1list.append(m1_all_list[i])
        m2list.append(m2_all_list[i])
        m0list.append(m0_all_list[i])
        dlist.append(d_all_list[i])

    del m1_all_list, m2_all_list, m0_all_list, d_all_list

    m0array = np.array(m0list)
    darray  = np.array(dlist)

    logm1xm2array = np.log(m1list) + np.log(m2list)

    myprint ('fitting model')

    a = np.array([logm1xm2array, darray])
    x = np.asmatrix(a).transpose()

    b = np.asarray([m0array])
    y = np.asmatrix(b).transpose()
    logy = np.log(y)

    regr = linear_model.LinearRegression()
    regr.fit(x, logy)
    logy_predict = regr.predict(x)
    y_predict = np.exp(logy_predict) 

    myprint ('finished fitting model')
    myprint ('Coefficients: \n')
    myprint (regr.coef_)
    r2 = r2_score(logy, logy_predict)
    myprint ('R squired = %.4f' % (r2))
    
    predict_overlap_bcd_cnt(bcd12_file_list, regr, bcd13_file_list, min_m1_value, min_m2_value) 

    myprint ('finished outputing bcd13 files')

    return min_m1_value, min_m2_value

def predict_overlap_bcd_cnt(bcd12_file_list, regr, bcd13_file_list, min_m1_value, min_m2_value):
    
    myprint ('calculating expected overlap barcode count...')
    for i in range(0, len(bcd12_file_list)):
        bcd12_file = bcd12_file_list[i]
        bcd13_file = bcd13_file_list[i]
        myprint ('current bcd12 file: %s, current output file: %s' %(bcd12_file, bcd13_file))
        bcd12_fp = open(bcd12_file, 'r')
        bcd13_fp = open(bcd13_file, 'w')

        pos_list = list()
        m1list = list()
        m2list = list()
        m0list = list()
        d1list = list()
        d2list = list()
        dlist = list()
        
        while 1:
            line = bcd12_fp.readline()
            if not line:
                break
            pos, m1, m2, m0, d1, d2 = line.strip().split(tab)

            pos = int(pos)
            m1 = float(m1) + 1e-6
            m2 = float(m2) + 1e-6
            m0 = float(m0)
            d1 = float(d1)
            d2 = float(d2)

            pos_list.append(pos)
            m1list.append(m1)
            m2list.append(m2)
            m0list.append(m0)
            d1list.append(d1)
            d2list.append(d2)
            dlist.append(d1+d2)


        m0array = np.array(m0list)
        darray = np.array(dlist)

        logm1xm2array = np.log(m1list) + np.log(m2list)

        a = np.array([logm1xm2array, darray])
        x = np.asmatrix(a).transpose()

        predict_logm0array = regr.predict(x) 
        predict_m0_array = np.exp(predict_logm0array)
        
        for i in range(0, len(pos_list)):
            predict_m0 = predict_m0_array[i]
            min_m1_m2 = min(m1list[i], m2list[i])
            if predict_m0 > min_m1_m2:
                predict_m0 = float(min_m1_m2)
            if predict_m0 < 0:
                predict_m0 = 0.0

            m0 = m0list[i] 
            if m0 == 0:
                m0 = 0.5
            if min_m1_m2 < 50:
                ratio = 1
            elif min_m1_m2 > 100:
                ratio = (predict_m0) / (m0)
            else:
                ratio1 = (predict_m0) / (m0) 
                ratio2 = 1
                ratio = (2 * min_m1_m2 - 100) *  ratio1 + (200 - 2 * min_m1_m2) * ratio2
                ratio = ratio / 100.0

            output = '%d\t%.f\t%.f\t%.f\t%.2f\t%.2f\t%.2f\t%.2f\n' %(pos_list[i], m1list[i]-1e-6, m2list[i]-1e-6, m0list[i], d1list[i], d2list[i], predict_m0, ratio)
            bcd13_fp.write(output)

        bcd12_fp.close()
        bcd13_fp.close()
             

    myprint ('finished output bcd13 files')
    return


if __name__ == '__main__':
    main()
