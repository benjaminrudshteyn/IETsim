#! /usr/bin/python

import sys
import re

def com_to_crystbind(file):
	""" 
        Assuming Crystal is cubic or rectangular, for surfaces
	"""
	fin = open(file,'r')
	paper = []
	#
	# Need only coordinates, filter out all line with less than
	# four entries
	for line in fin.readlines():
		newline = line.split()
		if len(newline) > 3:
			paper.append(newline)

	fin.close()

	# remove any gaussian keyword commands if remaining
	find_1 = re.compile('#')
	find_2 = re.compile('%')
	i = 0
	for line in paper:
		if find_1.search(line[0]):
			paper.pop(i)
		elif find_2.search(line[0]):
			paper.pop(i)
		i = i + 1

	# output filename is name of com file with .com removed
	fout = open(file[:-4]+'.bind','w')

	# Write title of Calculation
	fout.write(file+'\n\n')

	# Find Crystal Basis Vectors
	crystal_a = 0.0
	crystal_b = 0.0
	crystal_c = 0.0
	for atom in paper:
		if atom[0] == 'Tv':
			crystal_a += float(atom[1])
			crystal_b += float(atom[2])
			crystal_c += float(atom[3])
	#print 'a = %f b = %f c = %f' % (crystal_a,crystal_b,crystal_c)		

	# Count number of atoms
	number = len(paper)
	line = 'Geometry Crystallographic\n%d  \n' % number
	fout.write(line)

	# Write out geometry
	i = 1
	for atom in paper:
		if atom[0] == 'Tv':
			newline = '%3d  %2s  %9.6f  %9.6f  %9.6f\n' % (i,'&',(float(atom[1])+float(paper[0][1]))/crystal_a,(float(atom[2])+float(paper[0][2]))/crystal_b,(float(atom[3])+float(paper[0][3]))/crystal_c)
		else:
			newline = '%3d  %2s  %9.6f  %9.6f  %9.6f\n' % (i,atom[0],float(atom[1])/crystal_a,float(atom[2])/crystal_b,float(atom[3])/crystal_c)
		fout.write(newline)
		i = i + 1


	# Count the number of electrons
	line = '\nCharge\n 0 \n'
	fout.write(line)

	# write out Lattice definition
	dimension = 0
	for atom in paper:	
		if atom[0] == 'Tv':
			dimension = dimension + 1
	line = '\nLattice\n%d\n' % dimension
	fout.write(line)
		
	# Choose number of overlaps per lattice vector
	line = '%d %d %d \n' % (1,1,1)
	fout.write(line)
	
	# Lattice vectors
	number = len(paper)
	line = '1 %3d\n1 %3d\n1 %3d\n' % (number-2,number-1,number)
	fout.write(line+'\n')
	# end Lattice definition

	# Crystal Specs
	line = '\nCrystal Spec\n%f  %f  %f\n90  90  90\n' % (crystal_a,crystal_b,crystal_c)
	fout.write(line)

	# Do basic Analysis
	line = '\nAverage Properties\n\nkpoints\n1\n0.000000  0.000000  0.000000 1\n\n'
	fout.write(line)

	# The End!
	fout.close()


def p_table():
	dict = {
            'H' : 1,
            'Li': 1,
            'C' : 4,
            'N' : 5,
            'O' : 6,
            'F' : 7,
            'Na': 1,
            'Al': 3,
            'Si': 4,
            'P' : 5,
            'S' : 6,
            'Cl': 7,
            'Ti': 4,
            'Mn': 7,
            'Zn':12
            }
	return dict

if __name__=='__main__':
	com_to_crystbind(sys.argv[1])	
