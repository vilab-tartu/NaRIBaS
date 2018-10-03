import os
import fileinput
import shutil
files = ['.edr','.gro','.log','.tpr','.trr', '.xtc', '_pullf.xvg', '_pullx.xvg']
for i in range(len(files)):
  cmd5 = "sed -e 's/.tpr/"+files[i]+"/g' < copyback.txt > copyback.sh"
  os.system(cmd5)
  cmd6 = 'chmod 777 copyback.sh'
  os.system(cmd6)
  cmd7 = 'sh copyback.sh'
  os.system(cmd7)
