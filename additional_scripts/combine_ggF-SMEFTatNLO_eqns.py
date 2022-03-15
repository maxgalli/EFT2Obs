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

def combineEqns(loop, tree, tree_loop_2, tree_loop_4):
  combined_eqn = OrderedDict()
  for stxs_bin in loop.keys():
    all_poss_keys = set(loop[stxs_bin].keys() + tree[stxs_bin].keys() + tree_loop_2[stxs_bin].keys() + tree_loop_4[stxs_bin].keys())
    params = filter(lambda x: x[:2] != "u_", all_poss_keys)

    combined_bin = OrderedDict()

    for param in params:
      loop_term = getTerm(loop[stxs_bin], param)
      tree_term = getTerm(tree[stxs_bin], param)
      tree_loop_2_term = getTerm(tree_loop_2[stxs_bin], param)
      tree_loop_4_term = getTerm(tree_loop_4[stxs_bin], param)

      if param == "B_cpg_2":
        tree_loop_4_term = [0, 0]

      combined_bin[param] = loop_term[0] + tree_term[0] + tree_loop_2_term[0] + tree_loop_4_term[0]
      combined_bin["u_"+param] = np.sqrt(loop_term[1]**2 + tree_term[1]**2 + tree_loop_2_term[1]**2 + tree_loop_4_term[1]**2)

    combined_eqn[stxs_bin] = combined_bin

  return combined_eqn

def getTerm(eqn, key):
  if key in eqn.keys(): return [eqn[key], eqn["u_"+key]]
  else: return [0,0]

def convert_SMEFTatNLO_To_SMEFTsim(eqns):
  a_s = 0.118
  g_s = 2*np.sqrt(a_s*np.pi)

  for stxs_bin in eqns.keys():
    params = filter(lambda x: x[:2] != "u_", eqns[stxs_bin].keys())
   
    converted_eqn = OrderedDict() 
    for param in sorted(params):
      term = getTerm(eqns[stxs_bin], param)
      param = param.replace("ctp", "cth")
      param = param.replace("cpg", "chg")
   
      if "ctg" in param:
        if "_2" not in param:
           term[0] /= -g_s
           term[1] /= g_s
        else:
           term[0] /= (-g_s)**2
           term[1] /= g_s**2
      converted_eqn[param] = term[0]
      converted_eqn["u_"+param] = term[1]
    eqns[stxs_bin] = converted_eqn
  return eqns

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

parser = argparse.ArgumentParser()
parser.add_argument('--output', '-o', default="combined_ggF-SMEFTatNLO.json")
parser.add_argument('--key', help="Key to interpret bin numbers.")
args = parser.parse_args()

with open(args.key, "r") as f:
  key = json.loads(f.read(), object_pairs_hook=keysToInt)
key = {i:key[tag] for i, tag in enumerate(key.keys())}

with open("ggF-SMEFTatNLO_loop.json", "r") as f:
  loop = jsonToNewDict(json.load(f), key)
with open("ggF-SMEFTatNLO_tree.json", "r") as f:
  tree = jsonToNewDict(json.load(f), key)
with open("ggF-SMEFTatNLO_tree_loop_2.json", "r") as f:
  tree_loop_2 = jsonToNewDict(json.load(f), key)
with open("ggF-SMEFTatNLO_tree_loop_4.json", "r") as f:
  tree_loop_4 = jsonToNewDict(json.load(f), key)

combined_eqn = combineEqns(loop, tree, tree_loop_2, tree_loop_4)
combined_eqn = convert_SMEFTatNLO_To_SMEFTsim(combined_eqn)
combined_eqn = addOtherParams(combined_eqn)

with open(args.output, "w") as f:
  json.dump(combined_eqn, f, indent=4)

