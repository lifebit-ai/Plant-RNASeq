#!/usr/bin/env python

import os
import argparse

parser = argparse.ArgumentParser(description='Computing gene expressions using IsoXtractor', epilog="", formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-w', '--window_size', type=int, required=False, help="Set window size of isochores in bp for IsoSegmenter (default: 100,000).")
args = parser.parse_args()
if args.window_size == None:
	args.window_size = 100000

os.makedirs( "graphs")
#run isosegmenter on all chromosomes 
os.makedirs( "isochores")
for filename in os.listdir('chromosomes/'):
	print("-Running isoSegmenter on " + str(filename) )	
	command = r'isoSegmenter.py ' + ' -i chromosomes/'+str(filename) + ' --y_min 1 --y_max 100 -g graphs/'+str(filename)+'.jpg' + ' --window_size ' + str(args.window_size) + ' -o isochores/'+str(filename)+'.csv'	
	os.system(command)