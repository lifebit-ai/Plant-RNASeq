import os
import sys
import random
import string
import argparse
import math

isochores = []
total_reads = []
isochore_family = []

conditions_file = open(args.conditions)

cond_count = 0

lines = conditions_file.readlines()

for i in range( 0, len(lines) ):
	no_newline = lines[i].strip('\n')
	cond = no_newline.split(',')
	if( int(cond[1]) > cond_count ):
		cond_count = int(cond[1])


conditions =[[] for i in xrange(cond_count)]


for i in range( 0, len(lines) ):
	no_newline = lines[i].strip('\n')
	cond = no_newline.split(',')
	conditions[ int(cond[1])-1 ].append( cond[0] )


#lines is number of read files
conditions_file.close()

input = open("reads.csv")
inputList = input.read()

seq = inputList.splitlines()


for i in range(0, len(seq)):
	isochores.append ( seq[i].split(',') )


#compute number of isochore families
for i in range(1, len(isochores) ):
	isochore_family.append ( isochores[i][1] ) #isochore family must be in column 1 of input

isochore_family = set( isochore_family )

iso_family = list( isochore_family )

################################################################################################

#identify which read header is in which column

conditions_val=[[] for i in xrange(cond_count)]

for i in range(0, len(conditions) ):
	for j in range( 0, len(conditions[i] ) ):
		conditions_val[i].append( 0 )

head = seq[0]
header = head.split(',')[6:]

for j in range(0, cond_count):
	for k in range( 0, len(conditions[j]) ):
		for l in range( 0, len(header) ):
			if( str( conditions[j][k] ) == str( header[l] ) ):
				conditions_val[j][k] = l + 6


#compute total reads for each column
for j in range ( 6, 6+len(lines) ):
	total = 0
	for k in range( 1, len( isochores) ):
		total = total + float(isochores[k][j])
	total_reads.append(total)

###################################################################################
expression = [0] * cond_count #columns

for j in range(0, cond_count ):
	expression[j] = [0] * len(isochores) #rows


#compute average for each of the conditions
for j in range(0, cond_count):
	for l in range (1, len(isochores) ):
		avg = 0
		#total_reads_counter = j*cond_count
		for k in range( 0, len(conditions_val[j]) ):
			formula = float(isochores[l][conditions_val[j][k]]) / ( total_reads[int(conditions_val[j][k])-6] * float(isochores[l][5]) ) 
			#lengths must be in column 5 of input

			if formula == 0:
				avg = avg + 0
			else: avg = avg + ( math.log( formula, 2.0 ) ) # compute log of formula above

		expression[j][l-1] = avg/ len(conditions_val[j])


final_table = [0] * ( 1+len(iso_family) )

for j in range(0, 1+len(iso_family) ):
	final_table[j] = [0] * ( 3 * cond_count + 2 )

#add titles to final table
final_table[0][0] = "isochore family"
c = 0
inp = 1
for j in range(0, 3*cond_count + 1 ):
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

final_table[0][3 * cond_count + 1 ] = " Count"

mean = [0] * len(iso_family)
for j in range(0, len(iso_family) ):
	mean[j] = [0] * cond_count

#compute mean
for i in range( 0, len(iso_family) ):
	for j in range(0, cond_count ):
		counter = 0
		total = 0
		for k in range( 1, len( isochores ) ):
			if isochores[k][1] == iso_family[i]: #isochore family in column 1
				counter = counter + 1
				total = total + expression[j][k-1]

		mean[i][j] = total / counter

#########################table 2 ####################################################



print("-Computing average expression for isochores.")


input_read = open("reads.csv")
input_readlist = input_read.read()

all_reads = input_readlist.splitlines()

avg_exp = open("avg_exp.csv", 'w')

split = all_reads[0].split(',')

for j in range(0, 6 ):
	avg_exp.write( split[j] + ',' )
for k in range(0, cond_count):
	avg_exp.write( "condition " + str(k+1) + ',' )
avg_exp.write('\n')	


for i in range(1, len(isochores) ) :
	split = all_reads[i].split(',')
	for j in range(0, 6 ):
		avg_exp.write( split[j] + ',' )
	for k in range(0, cond_count):
		avg_exp.write( str(expression[k][i-1]) + ',' )
		
	avg_exp.write('\n')




###############################################################################

print("-Computing average expression for isochore classes.")	

#compute final table
for i in range( 0, len(iso_family) ):
	final_table[i+1][0] = iso_family[i]
	cond = 0
	for j in range(0, cond_count ):
		counter = 0
		standard_dev = 0
		variance = 0
		total = 0
		for k in range( 1, len( isochores ) ):
			if isochores[k][1] == iso_family[i]: #isochore family in column 1
				counter = counter + 1
				total = total + expression[j][k-1]
				standard_dev = standard_dev + ( ( expression[j][k-1] - mean[i][j] )** 2 )


		final_table[i+1][cond + 1] = mean[i][j] #average
		final_table[i+1][cond + 2] = math.sqrt( standard_dev / counter )
		final_table[i+1][cond + 3] = standard_dev / counter
		final_table[i+1][3*cond_count + 1] = counter #count
		cond = cond + 3


exp = open("expression.csv", 'w')
for i in range ( 1+len(iso_family) ):
	for j in range(0, 3 * cond_count + 2):
		exp.write( str(final_table[i][j])+',' )	
	exp.write("\n")
			

