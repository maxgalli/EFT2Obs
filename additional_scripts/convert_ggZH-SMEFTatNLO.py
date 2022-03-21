import json
from collections import OrderedDict
import argparse
import numpy as np
import os

def getTerm(eqn, key):
  if key in eqn.keys(): return [eqn[key], eqn["u_"+key]]
  else: return [0,0]

conversion = {
 "cll1221": "cll1",
 "cpdc": "chdd",
 "cdp": "chbox",
 "c3pl1": "chl3",
 "c3pl2": "chl3",
 "cpqm": "chq1",
 "cpt": "cht",
 "ctp": "cthre",
 "cpd": "chd",
 "cpu": "chu"
}

to_drop = ["cpq3", "i"]

def convert_SMEFTatNLO_To_SMEFTsim(eqns):
  for stxs_bin in eqns.keys():
    params = filter(lambda x: x[:2] != "u_", eqns[stxs_bin].keys())
   
    converted_eqn = OrderedDict() 
    for param in sorted(params):
      drop_param = False
      for to_drop_param in to_drop:
        if to_drop_param in param:
          drop_param = True

      if not drop_param:
        term = getTerm(eqns[stxs_bin], param)
        for key, value in conversion.items():
          param = param.replace(key, value)
          if param == "B_chl3_chl3":
            param = "B_chl3_2"
  
        if param not in converted_eqn.keys():
          converted_eqn[param] = 0
          converted_eqn["u_"+param] = 0
      
        converted_eqn[param] += term[0]
        converted_eqn["u_"+param] = np.sqrt(converted_eqn["u_"+param]**2 + term[1]**2)
    eqns[stxs_bin] = converted_eqn
  return eqns

parser = argparse.ArgumentParser()
parser.add_argument('--input-dir', '-i', type=str, default="ConvertedEquations/")
parser.add_argument('--output', '-o', default="ConvertedEquations/ggZH_SMEFTatNLO_combined.json")
parser.add_argument('--postfix', type=str, default="")
args = parser.parse_args()

with open(os.path.join(args.input_dir, "ggZH_SMEFTatNLO_%s.json"%args.postfix), "r") as f:
  eqns = json.load(f)

converted_eqns = convert_SMEFTatNLO_To_SMEFTsim(eqns)

with open(args.output, "w") as f:
  json.dump(converted_eqns, f, indent=4)

