#!/usr/bin/env python
import os
import sys
import re

isochores = []

input = open(sys.argv[1], 'r') # output from IsoSegmenter
inputList = input.read()

seq = inputList.splitlines()

for i in range(1, len(seq)):
	isochores.append ( seq[i].split(',') )

reads = []
input2 = open(sys.argv[2], 'r') # output from REAL
inputList2 = input2.read()
seq2 = inputList2.splitlines()
for i in range(1, len(seq2)):
	reads.append ( seq2[i].split() )

counter = 0

no = 0;
for i in range(0, len(isochores)):
	if ( isochores[i][3] != "gap" ):
		counter = 0
		for j in range(0, len(reads)):
			start = float(reads[j][9])
			end = float(reads[j][9]) + float(reads[j][6])

			if start > float(isochores[i][1]):
				break

			if end > float(isochores[i][0]):
				if float(reads[j][9]) <  float(isochores[i][0]) and (float(reads[j][9]) + float(reads[j][6])  > float(isochores[i][0])) and   ( ( float(reads[j][9]) + float(reads[j][6]) ) -  float(isochores[i][0]) ) > ( float(reads[j][9]) + float(reads[j][6]) )/2:
					counter = counter +1
				elif float(reads[j][9]) > float(isochores[i][0]) and ( ( float(reads[j][9]) +  float(reads[j][6]) ) < float(isochores[i][1]) ):
					counter = counter +1
				elif float(reads[j][9]) < float(isochores[i][1]) and (float(reads[j][9]) + float(reads[j][6]) )> float(isochores[i][1]) and float(isochores[i][1]) - float(reads[j][9]) > ( float(reads[j][9]) + float(reads[j][6]) )/2:
					counter=counter+1


		print( counter )
