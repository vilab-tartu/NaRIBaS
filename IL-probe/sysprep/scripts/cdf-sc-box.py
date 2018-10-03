import os
import fileinput
import load_cube
import math
import numpy
import shutil
import sys

### Get distances
dist = []
for folder in os.listdir('./'):
  if os.path.isfile(folder): pass
  else:
    if "." in folder:
      dist.append(str(folder))
dist.sort()

### Calculate
repn = [0, 1, 2]

for i in range(len(dist)):

  #COPY tpr and index  
  shutil.copyfile(str(dist[i])+'/'+str(repn[0])+'/NVT.tpr',str(dist[i])+'/NVT.tpr')
  shutil.copyfile(str(dist[i])+'/'+str(repn[0])+'/index.ndx',str(dist[i])+'/index.ndx')
  merge = 'trjcat -f '+str(dist[i])+'/0/NVT.xtc '+str(dist[i])+'/1/NVT.xtc '+str(dist[i])+'/2/NVT.xtc -o  '+str(dist[i])+'/merge.xtc'
  os.system(merge)

  #CREATE cube file

  conv = 'trjconv -s '+str(dist[i])+'/NVT.tpr -f '+str(dist[i])+'/merge.xtc -n '+str(dist[i])+'/index.ndx -o '+str(dist[i])+'/noPBC.xtc -pbc mol -ur compact -center << EOF\n6\n0\nEOF'
  os.system(conv)

  conv = 'trjconv -s '+str(dist[i])+'/NVT.tpr -f '+str(dist[i])+'/noPBC.xtc -n '+str(dist[i])+'/index.ndx -o '+str(dist[i])+'/fit.xtc -fit rot+trans << EOF\n0\n0\nEOF'
  os.system(conv)

  sdf = 'g_spatial -f '+str(dist[i])+'/fit.xtc -s '+str(dist[i])+'/NVT.tpr -n '+str(dist[i])+'/index.ndx -nab 20 << EOF\n3\n5\nEOF'
  mv = 'mv grid.cube '+str(dist[i])+'/sdf_BMI.cube'
  os.system(sdf)
  os.system(mv)
  
  sdf = 'g_spatial -f '+str(dist[i])+'/fit.xtc -s '+str(dist[i])+'/NVT.tpr -n '+str(dist[i])+'/index.ndx -nab 20 << EOF\n4\n5\nEOF'
  mv = 'mv grid.cube '+str(dist[i])+'/sdf_BF4.cube'
  os.system(sdf)
  os.system(mv)

  #PLOT rdf

  rdf = 'g_rdf -com -f '+str(dist[i])+'/fit.xtc -s '+str(dist[i])+'/NVT.tpr -n '+str(dist[i])+'/index.ndx -o '+str(dist[i])+'/rdf_3D_BMI.xvg << EOF\n5\n3\nEOF'
  os.system(rdf)

  rdf = 'g_rdf -com -f '+str(dist[i])+'/fit.xtc -s '+str(dist[i])+'/NVT.tpr -n '+str(dist[i])+'/index.ndx -o '+str(dist[i])+'/rdf_3D_BF4.xvg << EOF\n5\n4\nEOF'
  os.system(rdf)

  #PLOT cylindrical distribution function
  cube_an =load_cube.CUBE(str(dist[i])+'/sdf_BF4.cube')
  cube_ca =load_cube.CUBE(str(dist[i])+'/sdf_BMI.cube')
  f = open(str(dist[i])+'/charge_radial_2D_box.dat', 'w+')

  r_an = []
  r_ca = []
  d_an = []
  d_ca = []
  d_an_norm = []
  d_ca_norm = []
  d_rev = []

  n_an = float(sys.argv[1]) #$anion_num_name
  n_ca = float(sys.argv[2]) #$cation_num_name

  sum_an = cube_an.data.sum()
  sum_ca = cube_ca.data.sum()
  
  axis=40 #2nm
  z = 80
  
  for j in range(z):
    #anion
    for k in range(cube_an.NY):
      for l in range(cube_an.NX):
        r_an.append(math.sqrt((k-cube_an.NY/2)**2+(l-cube_an.NX/2)**2))
        d_an.append(cube_an.data[l][k][j+cube_an.NZ/2-z]*n_an/sum_an)
      rd_an = zip(r_an,d_an)
      rd_an.sort(key=lambda element:element[0])
    #cation
    for k in range(cube_ca.NY):
      for l in range(cube_ca.NX):
        r_ca.append(math.sqrt((k-cube_ca.NY/2)**2+(l-cube_ca.NX/2)**2))
        d_ca.append(cube_ca.data[l][k][j+cube_ca.NZ/2-z]*n_ca/sum_ca)
      rd_ca = zip(r_ca,d_ca)
      rd_ca.sort(key=lambda element:element[0])

    d_el = 0
    #anion
    R = 1
    for m in range(len(rd_an)):
      if rd_an[m][0] <= R:
        d_el += rd_an[m][1]/(R**2-(R-1)**2)
      else:
        d_an_norm.append(d_el)
        R += 1
        d_el = 0
    R = 1
    #cation
    R = 1
    for m in range(len(rd_ca)):
      if rd_ca[m][0] <= R:
        d_el += rd_ca[m][1]/(R**2-(R-1)**2)
      else:
        d_ca_norm.append(d_el)
        R += 1
        d_el = 0

    for n in range(axis/2):
      d_rev.append(d_ca_norm[n] - d_an_norm[n])
    d_rev[:] = d_rev[::-1]
    for n in range(axis/2-1):
      print >>f, d_rev[n],
    for n in range(axis/2):
      print >>f, d_ca_norm[n] - d_an_norm[n],
    print >>f, ''
    r_an[:] = []
    r_ca[:] = []
    d_an[:] = []
    d_ca[:] = []
    d_an_norm[:] = []
    d_ca_norm[:] = []
    d_rev[:] = []
  f.close()
  
  filename = 'charge_radial_2D_box.dat'  
  name = 'charge_radial_2D_box'            #get file name without extension
  wall = "{0:.5f}".format(float(dist[i])) #get the distance between the electrode and the wall from the folder name, e.g. 0.20 -> 0.20000
  if str(sys.argv[3]) is "PLi":
    diameter = 1.78
  elif str(sys.argv[3]) is "PKA":
    diameter = 3.04
  else:
    diameter = 3.04
  probe= "{0:.5f}".format(float(axis))
  for line in fileinput.FileInput('simple_contour_plot_box.m',inplace=1):
    line = line.replace('FILE_NAME', str(dist[i])+'/'+filename).replace('NAME', name).replace('MOLS','d='+str(dist[i])).replace('CENTER', str(probe)).replace('WALL', wall).replace('DIAMETER',str(diameter))
    print line,                   #replace variable in matlab-script

  matlab = 'matlab -nosplash -nodesktop -r simple_contour_plot_box'
  os.system(matlab)               #plot figure and save it

  for line in fileinput.FileInput('simple_contour_plot_box.m',inplace=1):
    line = line.replace(str(dist[i])+'/'+filename, 'FILE_NAME').replace(name, 'NAME').replace('d='+str(dist[i]), 'MOLS').replace(str(probe), 'CENTER').replace(wall, 'WALL').replace(str(diameter),'DIAMETER')
    print line,                   #revert the matlab-script to its original state

#CLEAN UP
rm = 'find . -name \*# -type f -delete'
os.system(rm)
rm = 'rm */noPBC.xtc */merge.xtc'
os.system(rm)
