#!/usr/bin/env python

import os
import sys
import random
import string
import argparse
import math

#computes number of reads for each chromosome and isochore file

#print titles for number of reads file Chromosome
temp1 = open("temp1.csv", 'a')
temp1.write("Chromosome,Isochore Class,GC Level,Isochore Start,Isochore End,Isochore Size\n")

#filling in column by column

for isochore in os.listdir('isochores'):
	isochore_file = open(os.path.join('isochores', isochore ), "r").read() #output from isoSegmenter
	isochores = isochore_file.splitlines()
		
	for i in range(1, len(isochores)):
		iso = isochores[i].split(',')

		if ( iso[3] != "gap" ): 
			counter = 0

			temp1.write( str(isochore).split('.')[0][2:] +','+ iso[3]  +','+ iso[4]  +','+ iso[0] +','+ iso[1]  +','+ iso[2] +'\n' )


temp1.flush()
temp1.close()

for aligned in os.listdir('aligned'):
	aligned_file = open(os.path.join('aligned',aligned), "r").read() #output from REAL
	print("-Computing reads for " + str(aligned) )

	reads = aligned_file.splitlines()

	temp2 = open("temp2.csv", 'a')
	temp2.write(aligned.split('.')[0]+'\n' ), #header for reads

	for isochore in os.listdir('isochores'):

		print("-Computing reads for " + str(isochore) )	
		isochore_file = open(os.path.join('isochores', isochore ), "r").read() #output from isoSegmenter
		isochores = isochore_file.splitlines()

		for i in range(1, len(isochores)):
			iso = isochores[i].split(',')
			counter = 0
			
			if ( iso[3] != "gap" ): 
				for j in range(0, len(reads) ):
		
					re = reads[j].split()
				
					if( re[8].endswith( str(isochore).split('.')[0] ) ):
						if float(re[9]) <  float(iso[0]) and (float(re[9]) + float(re[6])  > float(iso[0])) and   ( ( float(re[9]) + float(re[6]) ) -  float(iso[0]) ) > ( float(re[9]) + float(re[6]) )/2:
							counter = counter +1 
						elif float(re[9]) > float(iso[0]) and ( ( float(re[9]) +  float(re[6]) ) < float(iso[1]) ):
							counter = counter +1
						elif float(re[9]) < float(iso[1]) and (float(re[9]) + float(re[6]) ) > iso[1] and float(iso[1]) - float(re[9]) > ( float(re[9]) + float(re[6]) )/2:
							counter=counter+1
		

				temp2.write( str(counter)+'\n')	
				temp2.flush()

	
	temp2.flush()
	tempOne = open("temp1.csv", 'r').read()	
	tempTwo = open("temp2.csv", 'r').read()
	temp3 = open("temp3.csv", 'w')

	lines = tempOne.splitlines()
	lines2 = tempTwo.splitlines()
	
	for i in range ( 0, len(lines) ):
		temp3.write( lines[i] + ',' + lines2[i] + '\n')
		temp3.flush()

	del(lines)
	del(lines2)
	os.remove("temp1.csv")
	os.remove("temp2.csv")
	os.rename("temp3.csv", "temp1.csv") 
	

os.rename("temp1.csv", "reads.csv") 
