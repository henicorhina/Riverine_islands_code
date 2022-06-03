#!/usr/bin/env python

"""
Created on Tue May 28 12:49:21 2019

@author: ojohnson
"""

import os
import argparse
import numpy as np

#os.chdir('/Volumes/Brumfield_Lab_Drive/River_islands/1_analysis/dxy_fst_final/results/conbic/')

def get_args():
	parser = argparse.ArgumentParser(
			description="""Program description""")
	parser.add_argument(
			"in_dir",
			type=str,
			help="""The input directory of results.
                       Must contain four files, the output of among_within_indiv_dists.py"""
		)
	parser.add_argument(
			"out_file",
			type=str,
			help="""The desired output file"""
		)
	return parser.parse_args()


def main():
    args = get_args()
    os.chdir(args.in_dir)
    #get number of columns
    x = np.genfromtxt('p_distances.txt',dtype='str', skip_header=1)
    cols = len(x[0]) - 1

    p_dist = np.loadtxt('p_distances.txt', usecols=range(1,cols), skiprows=1)
    p_dist[p_dist==0] = np.nan
    av_p_dist = np.mean(np.nanmean(p_dist, axis=1))
    print("p_distances: {}".format(av_p_dist))


    p_std = np.loadtxt('p_std.txt', usecols=range(1,cols), skiprows=1)
    p_std[p_std==0] = np.nan
    av_p_std = np.mean(np.nanmean(p_std, axis=1))
    print("p_std: {}".format(av_p_std))


    jc_dist = np.loadtxt('jc_distances.txt', usecols=range(1,cols), skiprows=1)
    jc_dist[jc_dist==0] = np.nan
    av_jc_dist = np.mean(np.nanmean(jc_dist, axis=1))
    print("jc_distances: {}".format(av_jc_dist))


    jc_std = np.loadtxt('jc_std.txt', usecols=range(1,cols), skiprows=1)
    jc_std[jc_std==0] = np.nan
    av_jc_std = np.mean(np.nanmean(jc_std, axis=1))
    print("jc_std: {}".format(av_jc_std))

    out = open(args.out_file, 'wb')
    out.write("P_distance\t{0}\n".format(av_p_dist))
    out.write("P_stdev\t{0}\n".format(av_p_std))
    out.write("Jc_distance\t{0}\n".format(av_jc_dist))
    out.write("Jc_stdev\t{0}\n".format(av_jc_std))
    out.close()

if __name__ == '__main__':
    main()
