import json
from collections import OrderedDict
import argparse
import numpy as np

#helper function for reading json and converting keys to integers
def keysToInt(x):
  od = OrderedDict()
  for k, v in x:
    od[int(k)] = v
  return od

def jsonToNewDict(eft2obs_json, key):
  new_json_dict = OrderedDict()
  
  tags = [edges[0] for edges in eft2obs_json['edges']]

  for i, tag in enumerate(tags):
    if key is not None:
      if tag not in key.keys(): #if this bin shouldn't exist
        warnings.warn("Found a bin in EFT2Obs json that does not exist in the key. Expected (and probably isn't a problem) if >1 STXS stage 0 processes in a dataset, e.g. VH")
        continue

    bin_dict = OrderedDict()
    for param_info in eft2obs_json['bins'][i]:
      #linear terms
      if len(param_info) == 3:
        bin_dict["A_%s"%param_info[2]] = param_info[0]
        bin_dict["u_A_%s"%param_info[2]] = param_info[1]
      #quadratic terms
      elif param_info[2]==param_info[3]:
        bin_dict["B_%s_2"%param_info[2]] = param_info[0]
        bin_dict["u_B_%s_2"%param_info[2]] = param_info[1]
      #cross terms
      else:
        bin_dict["B_%s_%s"%(param_info[2], param_info[3])] = param_info[0]
        bin_dict["u_B_%s_%s"%(param_info[2], param_info[3])] = param_info[1]
    
    #if key available use it
    if key is not None:
      name = key[tag]
    else:
      name = tag
    new_json_dict[name] = bin_dict

  return new_json_dict

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
 "ctp": "cth",
 "cpd": "chd", 
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

"""
def addOtherParams(eqns):
  v = 0.24622
  params = OrderedDict()
  params["chbox"] = 1
  params["chdd"] = -0.25
  params["chl3"] = -1
  params["cll1"] = 0.5

  print(params.keys())
  n_par = len(params)

  stxs_bins = eqns.keys()

  for stxs_bin in stxs_bins:
    if eqns[stxs_bin].keys() == []: continue
    original_params = filter(lambda x: x[:2] == "A_", eqns[stxs_bin].keys())
    original_params = [param[2:] for param in original_params]

    for param in params.keys():
      eqns[stxs_bin]["A_%s"%param] = params[param]*2*v**2
      eqns[stxs_bin]["u_A_%s"%param] = 0

    for i, p1 in enumerate(params.keys()):
      for j, p2 in enumerate(params.keys()):
        if p1 == p2:
          eqns[stxs_bin]["B_%s_2"%p1] = params[p1]*params[p1]*v**4
          eqns[stxs_bin]["u_B_%s_2"%p1] = 0
        elif j > i:
          eqns[stxs_bin]["B_%s_%s"%(p1, p2)] = params[p1]*params[p2]*2*v**4
          eqns[stxs_bin]["u_B_%s_%s"%(p1, p2)] = 0
    
      for p3 in original_params:
        eqns[stxs_bin]["B_%s_%s"%(p1, p3)] = (params[p1]*2*v**2 * eqns[stxs_bin]["A_%s"%p3]) / 2
        eqns[stxs_bin]["u_B_%s_%s"%(p1, p3)] = (params[p1]*2*v**2 * eqns[stxs_bin]["u_A_%s"%p3]) / 2
  return eqns
"""

parser = argparse.ArgumentParser()
parser.add_argument('--output', '-o', default="ggZH-SMEFTatNLO_converted.json")
parser.add_argument('--key', help="Key to interpret bin numbers.")
args = parser.parse_args()

with open(args.key, "r") as f:
  key = json.loads(f.read(), object_pairs_hook=keysToInt)
key = {i:key[tag] for i, tag in enumerate(key.keys())}

with open("ggZH-SMEFTatNLO.json", "r") as f:
  eqns = jsonToNewDict(json.load(f), key)

converted_eqns = convert_SMEFTatNLO_To_SMEFTsim(eqns)
#combined_eqn = addOtherParams(combined_eqn)

with open(args.output, "w") as f:
  json.dump(converted_eqns, f, indent=4)

