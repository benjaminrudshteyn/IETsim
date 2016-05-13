#!/usr/bin/env python
from sys import argv
import re





def get_coords(filename):
	paper = open(filename).readlines()

	find_coords = re.compile('[A-Z][a-u,w-z]?\s+-?\d+.\d+\s+-?\d+.\d+\s+-?\d+.\d+')
	coords = [line.split() for line in paper if find_coords.search(line)]
	coords = [[line[0].upper()]+[float(x) for x in line[1:]] for line in coords]

	return coords


def get_basis(filename):
	paper = open(filename).readlines()
	
	basis = [line.split() for line in paper]
	basis = [line for line in basis if len(line)>5]
	basis = [line for line in basis if line[0]!=';']	
	basis = [[line[0]]+[line[5]] for line in basis]


	b_set = {}
	for atom in basis:
		if atom[1] == 's':
			b_set[atom[0]] = 1
		elif atom[1] == 'p':
			b_set[atom[0]] = 4
		elif atom[1] == 'd':
			b_set[atom[0]] = 9
		elif atom[1] == 'f':
			b_set[atom[0]] = 23
		
	return b_set


def basis_numbers(atom_num,coords,basis):
	b_nums = []
	
	num_start = 0
	for i in range(atom_num-1):
		num_start += basis[coords[i][0]]
	
	for i in range(1,basis[coords[atom_num-1][0]]+1):
		b_nums.append(num_start+i)

	return b_nums
		



if __name__=='__main__':
	xyzname = argv[1]
	coords = get_coords(xyzname)

	#basisname = argv[2]	
	basisname  = argv[0].split('/')
	basisname[-1] = 'eht_parms.dat'
	basisname = '/'.join(basisname)
	basis  = get_basis(basisname)
	
	nums = []
	for i in argv[2:]:
		atom_num = int(i)
		nums += basis_numbers(atom_num,coords,basis)

	outline = ''
	print 'Absorbing'
	print '%d %4.3f' % (len(nums),1.0)
	for i in range(len(nums)):
		if i%10==9:
			outline += '%d   \\\n' % nums[i]
		elif i==len(nums)-1:
			outline += '%d\n' % nums[i]
		else:
			outline += '%d, ' % nums[i]
	print outline
	print argv[3]

