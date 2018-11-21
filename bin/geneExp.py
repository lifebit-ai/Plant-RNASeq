#!/usr/bin/env python
import sys
import math

isochores = []
total_reads = []
isochore_family = []

conditions = int( sys.argv[2] )

input = open(sys.argv[1], 'r')
inputList = input.read()

seq = inputList.splitlines()


for i in range(0, len(seq)):
	isochores.append ( seq[i].split(',') )

arguments = len(sys.argv) - 1


#compute number of isochore families
for i in range(1, len(isochores) ):
	isochore_family.append ( isochores[i][2] ) #isochore family must be in column 2 of input

isochore_family = set( isochore_family )

iso_family = list( isochore_family )

#identify column ranges of each condition
con = 3
column_ranges = []

for j in range(0, conditions ):
	column_ranges.append( sys.argv[con].split(',') )
	con = con + 1

#compute total reads for each column
for j in range( int( column_ranges[0][0] ), int( column_ranges[ len( column_ranges ) - 1 ][1] )+1 ):
	total = 0
	for k in range( 1, len( isochores) ):
		if(  float(isochores[k][j]) > 10 ):
			total = total + float(isochores[k][j])
	total_reads.append(total)


expression = [0] * conditions #columns

for j in range(0, conditions ):
	expression[j] = [0] * len(isochores) #rows


#compute average for each for the conditions
for j in range(0, conditions):
	for l in range (1, len(isochores) ):
		avg = 0
		total_reads_counter = j*conditions
		for k in range( int(column_ranges[j][0]), int(column_ranges[j][1])+1 ):

			if(  float(isochores[l][k]) > 10 ):
				formula = float(isochores[l][k]) / ( total_reads[ total_reads_counter ] * float(isochores[l][5]) ) #lengths must be in column 5 of input

				#print( float(isochores[l][k]) , " ", total_reads[ total_reads_counter ], " ", float(isochores[l][5]) )
				if formula == 0:
					avg = avg + 0
				else: avg = avg + ( math.log( formula, 2.0 ) )
				total_reads_counter = total_reads_counter = total_reads_counter + 1

		expression[j][l-1] = avg/(int(column_ranges[j][1])- int(column_ranges[j][0])+1 )
		#print( expression[j][l-1] )

final_table = [0] * ( 1+len(iso_family) )

for j in range(0, 1+len(iso_family) ):
	final_table[j] = [0] * ( 3 * conditions + 2 )

#add titles to final table
final_table[0][0] = "isochore family"
c = 0
inp = 1
for j in range(0, 3*conditions + 1 ):
	cond = "Condition "
	cond+= str(inp)
	if c == 0: cond+= " avg"
	if c == 1: cond+= " std dev"
	if c == 2: cond+= " var"
	if c == 2:
		inp = inp + 1
		c = 0
	else: c = c + 1
	final_table[0][j+1] = cond

final_table[0][3 * conditions + 1 ] = " Count"

mean = [0] * len(iso_family)
for j in range(0, len(iso_family) ):
	mean[j] = [0] * conditions

#compute mean
for i in range( 0, len(iso_family) ):
	for j in range(0, conditions ):
		counter = 0
		total = 0
		for k in range( 1, len( isochores ) ):
			if isochores[k][2] == iso_family[i]:
				counter = counter + 1
				total = total + expression[j][k-1]

		mean[i][j] = total / counter


#compute final table
for i in range( 0, len(iso_family) ):
	final_table[i+1][0] = iso_family[i]
	cond = 0
	for j in range(0, conditions ):
		counter = 0
		standard_dev = 0
		#variance = 0
		total = 0
		for k in range( 1, len( isochores ) ):
			if isochores[k][2] == iso_family[i]:
				counter = counter + 1
				total = total + expression[j][k-1]
				standard_dev = standard_dev + ( ( expression[j][k-1] - mean[i][j] )** 2 )

		#print( " ")
		#print( total )
		final_table[i+1][cond + 1] = mean[i][j] #average
		final_table[i+1][cond + 2] = math.sqrt( standard_dev / counter )
		final_table[i+1][cond + 3] = standard_dev / counter
		final_table[i+1][3*conditions + 1] = counter #count
		cond = cond + 3

for i in range ( 1+len(iso_family) ):
	print( final_table[i] )
