#!/usr/bin/env python

import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

try:
    from scripts import my_utils
except ImportError:
    import my_utils

tab  = '\t'
endl = '\n'
arg = sys.argv[1:]

usage = 'python ' + __file__ + ' ' + '<input.2d_overlapping_barcodes.txt> <out_dir> <bin_size>'
argc  = 3 

def plot_2d_overlapping_barcodes(in_2d_values_file, target_region_bedpe_list, bin_size, max_ovl_num, out_dir, out_prefix):

    my_utils.make_dir(out_dir)

    in_2d_values_fp = my_utils.gzopen(in_2d_values_file, 'r')

    region_title = ''
    xmin, xmax, ymin, ymax = (0, 0, 0, 0)
    xsize = 0
    ysize = 0
    current_x = 0
    current_y = 0

    ovl_2d_array = np.zeros(shape = (xsize, ysize), dtype=np.int32)
    
    while 1:
        line = in_2d_values_fp.readline()
        if not line: break
        line = line.strip()
        if line[0:2] == '##':
            if (xsize and ysize): plot_one_bedpe(out_dir, target_region_bedpe_list, out_prefix, region_title, ovl_2d_array, xmin, xmax, ymin, ymax, bin_size) 
            region_title = line[2:]
            xmin, xmax, ymin, ymax = (0, 0, 0, 0)
            xsize = 0
            ysize = 0
            current_x = 0
            current_y = 0
            continue
        if line[0:5] == '#xmin':
            line = line.split('=')
            xmin,xmax,ymin,ymax = line[1].split(',')
            xmin = int(xmin)
            xmax = int(xmax)
            ymin = int(ymin)
            ymax = int(ymax)
            xsize = int((xmax - xmin)/bin_size)
            ysize = int((ymax - ymin)/bin_size)
            ovl_2d_array = np.zeros(shape = (xsize, ysize), dtype=np.int32)
            current_x = 0
            current_y = 0
            continue
        line = line.split(tab)
        for current_y in range(0, len(line)):
            ovl_2d_array[current_x][current_y] = int(line[current_y])
            if ovl_2d_array[current_x][current_y] > max_ovl_num: 
                ovl_2d_array[current_x][current_y] = max_ovl_num
        current_x += 1
   
    if (xsize and ysize): plot_one_bedpe(out_dir, target_region_bedpe_list, out_prefix, region_title, ovl_2d_array, xmin, xmax, ymin, ymax, bin_size) 
    
    in_2d_values_fp.close()
    return

def find_sv_id(xchr, xstart, xend, ychr, ystart, yend, target_region_bedpe_list):

    sv_id = 'UNK'
    svtype = 'UNK'
    xstart = int(xstart)
    xend = int(xend)
    ystart = int(ystart)
    yend = int(yend)

    for i in range(0, len(target_region_bedpe_list)):
        bedpe1 = target_region_bedpe_list[i]
        if xchr == bedpe1.chrm1 and xstart == bedpe1.start1 and xend == bedpe1.end1 and ychr == bedpe1.chrm2 and ystart == bedpe1.start2 and yend == bedpe1.end2:
            sv_id = bedpe1.sv_id
            svtype = bedpe1.svtype
            break
        if ychr == bedpe1.chrm1 and ystart == bedpe1.start1 and yend == bedpe1.end1 and xchr == bedpe1.chrm2 and xstart == bedpe1.start2 and xend == bedpe1.end2:
            sv_id = bedpe1.sv_id
            svtype = bedpe1.svtype
            break

    return sv_id, svtype

def plot_one_bedpe(out_dir, target_region_bedpe_list, out_prefix, region_title, ovl_2d_array, xmin, xmax, ymin, ymax, bin_size): 

    tp_ovl_2d_array = ovl_2d_array.transpose()

    xlab, ylab = region_title.split(';') 
    xlab = xlab.strip()
    ylab = ylab.strip()
    xchr, x_start_end = xlab.split(':')
    ychr, y_start_end = ylab.split(':')

    xstart, xend = x_start_end.split('-')
    ystart, yend = y_start_end.split('-')

    sv_id, svtype = find_sv_id(xchr, xstart, xend, ychr, ystart, yend, target_region_bedpe_list)
    if sv_id == 'UNK' or svtype == 'UNK':
        sv_id = '%s_%s_%s.%s_%s_%s' % (xchr, xstart, xend, ychr, ystart, yend)
        svtype = 'unknown_sv_type'

    xsize = int((xmax - xmin) / bin_size)
    ysize = int((ymax - ymin) / bin_size)

    pd_ovl_2d_array = pd.DataFrame(tp_ovl_2d_array)
    xticks_dict = dict()
    yticks_dict = dict()

    for i in range(0, xsize):
        xticks_dict[i] = xmin + i * bin_size 

    for i in range(0, ysize):
        yticks_dict[i] = ymin + i * bin_size 

    pd_ovl_2d_array = pd_ovl_2d_array.rename(columns=xticks_dict, index = yticks_dict)
    
    cmrmap_r = cm.get_cmap('brg_r', 1000)
    cmrmap_r_colors = cmrmap_r(np.linspace(0, 1, 1000))
    r50, g50, b50, a50 = cmrmap_r_colors[500]
    r = np.linspace(1.0, r50, 500)
    g = np.linspace(1.0, g50, 500)
    b = np.linspace(1.0, b50, 500)
    a = np.linspace(1.0, a50, 500)
    my_colors1 = np.array([r, g, b, a]).transpose()
    r = np.linspace(r50, 0.0, 500)
    g = np.linspace(g50, 0.0, 500)
    b = np.linspace(b50, 0.0, 500)
    a = np.linspace(a50, 1.0, 500)
    my_colors2 = np.array([r, g, b, a]).transpose()
    #my_colors = np.vstack((my_colors, cmrmap_r_colors[500:]))
    my_colors = np.vstack((my_colors1, my_colors2)) 

    my_cmap = ListedColormap(my_colors)

    out_file = os.path.join(out_dir, '%s.%s.heatmap.png' % (out_prefix, sv_id))
    plt.figure(figsize=(10, 10))
    ax = sns.heatmap(pd_ovl_2d_array, cmap=my_cmap, square=True, xticklabels=int(xsize/10), yticklabels=int(ysize/10) )
    ax.invert_yaxis()

    ax.axhline(y=0,     color='k',linewidth=2)
    ax.axhline(y=ysize, color='k',linewidth=2)
    ax.axvline(x=0,     color='k',linewidth=2)
    ax.axvline(x=xsize, color='k',linewidth=2)

    plt.axis([0, xsize, 0, ysize])
    plt.xlabel(xchr)
    plt.ylabel(ychr)
    plt.xticks(rotation='vertical') 
    plt.yticks(rotation='horizontal') 
    plt.title('Number of overlapping barcodes (%s, %s)' % (sv_id, svtype))
    plt.rcParams.update({'font.size': 12})
    plt.show()
    plt.savefig(out_file, dpi=200)
    plt.close('all')
    my_utils.myprint('saved figure: %s' % out_file)


    return

if __name__ == '__main__':
    main()
