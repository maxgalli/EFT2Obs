import os
import subprocess
from collections import OrderedDict as od
import json

grep_width = subprocess.check_output(["grep 'Width :' logs/*gridpack.log"], shell=True)
grep_cross = subprocess.check_output(["grep 'Cross-section :' logs/*gridpack.log"], shell=True)
grep = grep_width+grep_cross
print(grep)

widths = od()
for line in grep.split("\n")[:-1]:
  decay_mode = line.split("/")[1].split("SMEFT")[0][:-1]
  width = float(line.split(":")[2].split("+")[0])
  print(decay_mode, width)
  widths[decay_mode] = width

with open("widths_cross_sections.json", "w") as f:
  json.dump(widths, f, indent=4)
