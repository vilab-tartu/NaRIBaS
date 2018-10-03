import os
import fileinput
import shutil
import numpy
import math

### Get distances
dist = ["0.20", "0.25", "0.30", "0.35", "0.40", "0.45", "0.50", "0.55", "0.60", "0.65", "0.70", "0.75", "0.80", "0.85", "0.90", "0.95", "1.00", "1.05", "1.10", "1.15", "1.20"]

### Calculate
repn = [5]
avgf = []
stdf = []
avgg = []
stdg = []

o = open('output-pmf-fullstd-sc.dat', 'w')

for i in range(len(dist)):
  name = str(dist[i])
  force = []

  for j in range(len(repn)):
    for line in open('./'+name+'/'+str(repn[j])+'/NVT_pullf.xvg'):
      if not line.startswith(('#', '@')):                      #ignore lines starting with # and @
        columns = line.split('\t')
        if len(columns) >= 2:
          force.append(float(columns[1].rstrip('\n')))         #make an array of force values

  avgf.append(numpy.mean(force[0:]))
  stdf.append(numpy.std(force[0:]))
  force[:] = []

output = sorted(zip(dist,avgf,stdf), key=lambda dist: dist[0], reverse=True)

cumsum = 0
for f in range(len(output)-2):
  cumsum += (float(output[f][1]) + 4*float(output[f+1][1]) + float(output[f+2][1]))*(float(output[f][0])-float(output[f+1][0]))/6.0
  avgg.append(cumsum)
shift = numpy.mean(avgg[0:5])

for f in range(len(output)-2):
  errsum = (float(output[f][2]) + 4*float(output[f+1][2]) + float(output[f+2][2]))*(float(output[f][0])-float(output[f+1][0]))/6.0
  stdg.append(errsum)

for k in range(len(avgg)):
  o.write(str(output[k+1][0])+'\t'+str(avgg[k]+shift)+'\t'+str(stdg[k])+'\n')
