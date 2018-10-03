import os
import fileinput
import shutil
import math

### Get distances
dist = ["0.20", "0.25", "0.30", "0.35", "0.40", "0.45", "0.50", "0.55", "0.60", "0.65", "0.70", "0.75", "0.80", "0.85", "0.90", "0.95", "1.00", "1.05", "1.10", "1.15", "1.20"]

### Calculate
repn = [5]
avgf = []
stdf = []
avgg = []
stdg = []
data = []

for i in range(len(dist)):
  name = str(dist[i])
  mean = []

  for j in range(len(repn)):
    with open(name+'/'+str(repn[j])+'/density_cation.xvg') as file:
      for lines in file.readlines():
        line = lines.strip()
        if not line.startswith("#"):
          if not line.startswith("@"):  
            cols = line.split()
            data.append(cols)

  dis = [float(i) for i in zip(*data)[0]]
  den = [float(i) for i in zip(*data)[1]]

  #CALCULATE SUM OF ALL REPLICA VALUES

  for j in range(len(data)/len(dist)):
    for k in range(len(dist)-1):
      l = j+len(data)/len(dist)*(k+1)
      den[j] += den[l]

  #CALCULTE AVERAGE OF ALL REPLICA VALUES
  buffer_string=""
  for j in range(len(data)/len(dist)):
    buffer_string=buffer_string+str(dis[j])+" "+str( den[j]/float(len(dist)))+"\n"
f = open("density_cation.out", "w")
f.write(buffer_string)
f.close()
