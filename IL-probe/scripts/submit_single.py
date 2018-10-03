import os
import fileinput
import shutil

for i in range(24):
  r = str(i)
   
  ###### copy submission file

  for line in fileinput.input("template.sub", inplace=1):
    print line.replace("XXXXXX", 'NVT'+r),
  shutil.copyfile('template.sub', r+'/ch-r'+r+'.sub')
  for line in fileinput.input("template.sub", inplace=1):
    print line.replace('NVT'+r, "XXXXXX"),


  subm = 'qsub NVT'+r+'.sub'
  os.system(subm)
