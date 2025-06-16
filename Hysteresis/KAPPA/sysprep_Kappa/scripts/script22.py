# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 11:50:32 2015

@author: karl
"""

import os
import glob
import matplotlib
matplotlib.use('agg') 
import matplotlib.pyplot as plt
import csv
import numpy


numbers = []

#files = glob.glob('cr_*.xvg')
files = glob.glob('*/density*xvg')
results = []
print files


# extract number of ions in the box
particle_numbers = []
for fname in files:
    #part_num = fname.split('_')[2][:-4]
    part_num = fname.split('_')[1][:-4]
    particle_numbers.append(int(part_num))

print particle_numbers

# iterate over files
for i, n in zip(files, particle_numbers):
    numbers = []
    with open(i) as f:
        header = ['distance']
        for line in f:
            if line.startswith(' '):
                coords = line.split()
                numbers.append(float(coords[1]))
                header.append(coords[0])
    numbers = numpy.array(numbers)
    results.append(numbers/sum(numbers)*n)


# sort by distances
matrix = numpy.array(results)
corect_idxs = numpy.argsort(particle_numbers)
matrix = matrix[corect_idxs]
particle_numbers = sorted(particle_numbers)

# write results
with open('density.csv', 'w') as f:
    wr = csv.writer(f)
    wr.writerow(header)
    for n, row in zip(particle_numbers, list(matrix)):
        wr.writerow( [n] + list(row))
        
# integrate

# determine region


# only works if first simulation has accurate representation of mono-layer
first_line  =matrix[0, :]
#print list(enumerate(first_line))
begin = None
end = 0

prev_number = 1e-5
for i, num in enumerate(first_line):
    if num > 1e-5:
        end = i
    if num > prev_number and not begin:
        begin = i - 1
    prev_number = num
    
#print begin, end
# integrate

integrals = []
for row in matrix:
    integrals.append(sum(row[begin:end]))

#calculate n(monolayer)
nmono=0
for i,integral in enumerate(integrals):
 if (i+1)<len(integrals):
  if integral>integrals[i+1]:
   nmono=sum(integrals[i:-1])/len(integrals[i:-1])
   break
 else:
  nmono=max(integrals)

#plot
xy=[0,2*nmono]
yy=[nmono,nmono]
plt.plot(xy, xy, "--k",xy,yy,"--k", particle_numbers, integrals, "bo")
plt.ylabel('Numerical counterion density of the 1st layer')
plt.xlabel('Total number of counterions')
#plt.show()
plt.xlim((0,2*nmono))
plt.ylim((0,2*nmono))
if nmono==max(integrals):
 plt.title('Monolayer representation not accurate', fontsize=12)
else:
 f = open("nmono.txt","w")
 f.write(str(nmono))
plt.savefig('density_graph.svg')
