import os
import sys
import random
import string
import argparse
import math


parser = argparse.ArgumentParser(description='Computing gene expressions using IsoXtractor', epilog="", formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-g', '--genome_dir', type=str, required=True, help="The directory name where the genome file exists")
parser.add_argument('-r', '--reads_dir', type=str, required=True, help="The directory name where the read files exist.")
parser.add_argument('-c', '--conditions', type=str, required=True, help="The conditions file identifying which read belongs to which condition.")
parser.add_argument('-w', '--window_size', type=int, required=False, help="Set window size of isochores in bp for IsoSegmenter (default: 100,000).")
parser.add_argument('-s', '--seed_errors', type=int, required=False, help="Maximum number of seed errors for REAL (default: 2).")
parser.add_argument('-e', '--total_errors', type=int, required=False, help="Total number of errors for REAL (default: 5).")
parser.add_argument('-l', '--seed_length', type=int, required=False, help="Length of seed in bp for REAL (default: 32).")
parser.add_argument('-t', '--threads', type=int, required=False, help="Number of threads to use (default: 1).")

args = parser.parse_args()

if args.window_size == None:
	args.window_size = 100000
if args.seed_errors == None:
	args.seed_errors = 2
if args.total_errors == None:
	args.total_errors = 5
if args.seed_length == None:
	args.seed_length = 32
if args.threads == None:
	args.threads = 1


if not os.path.exists(args.genome_dir+'/'):
	print("-Error: Genome folder " + args.genome_dir+ " does not exist")
	exit()
else: 
	path = args.genome_dir+'/'  # this is a path given by user of where genome folder is

if not os.path.exists(args.reads_dir+'/'):
	print("-Error: Reads folder " + args.reads_dir+ " does not exist")
	exit()
else:
	path_reads = args.reads_dir+'/' # path where reads are contained
	
if not os.path.isfile(args.conditions):
	print("-Error: Conditions file " + args.conditions+" does not exist")	
	exit()

#check genome ends in .fa
for filename in os.listdir(path):
	if( str(filename).endswith(".fa") == 0 ):
		print( "-Error: Genome fasta filename must end in .fa")
		exit()

split genome into chromosomes and store in folder chromosomes
os.makedirs( "chromosomes")
print("-Splitting genome into chromosomes")
for filename in os.listdir(path):
	text_file = open(os.path.join(path, filename), "r").read() #file path of genome provided by user, open filename inside directory
	seq = text_file.splitlines()
	for k in range(0 , len(seq)):
		if(seq[k].startswith(">")):
			file_in = str(seq[k][-2:])
			text_file2 = open("chromosomes/ch"+str(seq[k][-2:]), "w")
			text_file2.close()
		text_file3 = open("chromosomes/ch"+str(file_in), "a")
		text_file3.write(seq[k]+"\n")
		text_file3.close()
	

#run REAL on all reads
os.makedirs( "aligned")
for filename_genome in os.listdir(path):
	for filename_read in os.listdir(path_reads): 
		print("-Running REAL on " + str(filename_read) )
		command = r'./real -t ' + path + str(filename_genome) + ' -s ' +  str(args.seed_errors) + ' -e ' + str(args.total_errors) + ' -l ' + str(args.seed_length) + ' -T ' + str(args.total_errors) + ' -p ' + path_reads+filename_read + ' -o aligned/'+str(filename_read)+'.OUT'
		os.system(command)


os.makedirs( "graphs")
#run isosegmenter on all chromosomes 
os.makedirs( "isochores")
for filename in os.listdir('chromosomes/'):
	print("-Running isoSegmenter on " + str(filename) )	
	command = r'./isoSegmenter/scripts/isoSegmenter.py ' + ' -i chromosomes/'+str(filename) + ' --y_min 1 --y_max 100 -g graphs/'+str(filename)+'.jpg' + ' --window_size ' + str(args.window_size) + ' -o isochores/'+str(filename)+'.csv'	
	os.system(command)


