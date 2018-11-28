#!/usr/bin/env python

import pandas as pd
import glob
import re

## read in all output csv files from isoSegmenter & noReads
isoSegmenter_chs = glob.glob('../ec2_results/isosegmenter/*.csv')
no_reads_chs = glob.glob("../ec2_results/no_reads/*.csv")

isoSegmenter_dfs = []
no_reads_dfs = []

## for all chrs set ch<number>_isoSegmenter equal to data frame containing the isoSegmenter output (eg ch01_isoSegmenter)
for isoSegmenter_ch in isoSegmenter_chs:
    ch = re.search('(ch[0-9]+)', isoSegmenter_ch).group(0)
    exec(ch + "_isoSegmenter" + " = pd.read_csv('../ec2_results/isosegmenter/SL2.31%s.fa.csv', sep=',')" % ch)

## save all isoSegmenter variables that were set in an array
for variable in dir():
    if variable.endswith('_isoSegmenter'):
        isoSegmenter_dfs.append(variable)

## for all chrs & all fastqs set data frames containing the noReads output equal to ch<numer>_SRR<number> eg ch01_SRR346617_no_reads
for no_reads_ch in no_reads_chs:
    ch_fastq = re.search('(ch[0-9]+_SRR[0-9]+)', no_reads_ch).group(0)
    exec(ch_fastq + "_no_reads" + " = pd.read_csv('../ec2_results/no_reads/no_reads_output_SL2.31%s.csv', sep=',')" % ch_fastq)

## save all no_reads variables that were set in an array
for variable in dir():
    if variable.endswith('_no_reads'):
        no_reads_dfs.append(variable)

## append no_reads_output to each isoSegmenter output file, as columns, for each FASTQ file
for isoSegmenter in isoSegmenter_dfs:
    counter = 0
    for no_reads in no_reads_dfs:
        if re.search('(ch[0-9]+)', isoSegmenter).group(0) == re.search('(ch[0-9]+)', no_reads).group(0):
            counter+=1
            vars()[isoSegmenter][('No of Reads %s' % counter)] = vars()[no_reads]

# example of what ^ evaluates to
# ch01_isoSegmenter['No of Reads 1'] = ch01_SRR346617_no_reads
# ch01_isoSegmenter['No of Reads 2'] = ch01_SRR346618_no_reads
# ch01_isoSegmenter['No of Reads 3'] = ch01_SRR346619_no_reads
# ch02_isoSegmenter['No of Reads 1'] = ch02_SRR346617_no_reads
# ch02_isoSegmenter['No of Reads 2'] = ch02_SRR346618_no_reads
# ch02_isoSegmenter['No of Reads 3'] = ch02_SRR346619_no_reads

#TODO: need to do everyhting below here for all files

## join all of the isoSegmenter_output files (which contain the no_reads_output)
for isoSegmenter in isoSegmenter_dfs:
    result = pd.concat([ch01_isoSegmenter, ch02_isoSegmenter], ignore_index=True)

## rename & reorder columns

## output csv
result.to_csv('hope_this_works.csv')
