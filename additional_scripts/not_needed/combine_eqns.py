import json
from collections import OrderedDict as od
import sys
import math

CLEAN=False

def weightedSumBin(eqns, weights):
  combined_eqn = od()
  for param in params:
    A = 0
    u_A = 0 #sum of weighted error squared
    B = 0
    u_B = 0
    for i, eqn in enumerate(eqns):
      if "A_%s"%param in eqn.keys():
        A += weights[i] * eqn["A_%s"%param]
        u_A += (weights[i] * eqn["u_A_%s"%param])**2
      if "B_%s_2"%param in eqn.keys(): 
        B += weights[i] * eqn["B_%s_2"%param]
        u_B += (weights[i] * eqn["u_B_%s_2"%param])**2
    if A != 0: 
      combined_eqn["A_%s"%param] = A / sum(weights)
      combined_eqn["u_A_%s"%param] = math.sqrt(u_A) / sum(weights)
    if B != 0: 
      combined_eqn["B_%s_2"%param] = B / sum(weights)
      combined_eqn["u_B_%s_2"%param] = math.sqrt(u_B) / sum(weights)

  for i in range(len(params)):
    for j in range(i+1, len(params)):
      B = 0
      u_B = 0
      for k, eqn in enumerate(eqns):
        if "B_%s_%s"%(params[i],params[j]) in eqn.keys(): 
          B += weights[k] * eqn["B_%s_%s"%(params[i],params[j])]
          u_B += (weights[k] * eqn["u_B_%s_%s"%(params[i],params[j])])**2 
      if B != 0: 
        combined_eqn["B_%s_%s"%(params[i],params[j])] = B / sum(weights)
        combined_eqn["u_B_%s_%s"%(params[i],params[j])] = math.sqrt(u_B) / sum(weights)

  return combined_eqn

def jsonToText(eqns):
  string = ""
  for bin_name in eqns.keys():
    string += bin_name + ": 1 + "
    for term in eqns[bin_name].keys():
      if term[:2] != "u_":
        string += "+ %.2g %s "%(eqns[bin_name][term], term[2:])
    string += "\n"
  return string    

def cleanUp(new_json_dict):
  for tag in new_json_dict.keys():
    params = new_json_dict[tag].keys()
    params = filter(lambda x: x[0] != "u", params)
    for param in params:
      param_value = new_json_dict[tag][param]
      if abs(param_value) < 0.001:
        del new_json_dict[tag][param]
        del new_json_dict[tag]["u_"+param]
  return new_json_dict

prod_modes = ["bbH",  "ggH",  "qqH",  "tHq",  "tHW",  "ttH",  "WH_lep",  "ZH_lep"]

decay_modes = ['H_4j', 'H_4l', 'H_4tau', 'H_4v', 'H_aa', 'H_bb', 'H_cc', 'H_gg', 'H_jjll', 'H_jjlv', 'H_jjvv',
               'H_lltautau', 'H_llvv', 'H_lvlv', 'H_lvtauv', 'H_mumu', 'H_tautau', 'H_tauvtauv', 'H_za']

with open("config_christmas.json", "r") as f:
  config = json.load(f)
params = [each['name'] for each in config['parameters']]

eqns = od()

for mode in prod_modes+decay_modes:
  with open("converted_eqns/%s_SMEFT.json"%mode, "r") as f:
    eqns[mode] = json.load(f, object_pairs_hook=od)

with open("widths_cross_sections.json", "r") as f:
  weights = json.load(f, object_pairs_hook=od)

"""
For production, we can just add up all the json files except...
- tHq and tHW which must be summed up according to their cross sections
- gg->zh gets equations from qq->zh
"""

TH_FWDH_eqn = weightedSumBin([eqns["tHq"]["TH_FWDH"], eqns["tHW"]["TH_FWDH"]], [weights["tHq"], weights["tHW"]])
TH_eqn = weightedSumBin([eqns["tHq"]["TH"], eqns["tHW"]["TH"]], [weights["tHq"], weights["tHW"]])
del eqns["tHq"]
del eqns["tHW"]
eqns["tH"] = od()
eqns["tH"]["TH_FWDH"] = TH_FWDH_eqn
eqns["tH"]["TH"] = TH_eqn

eqns["ggZH"] = od()
for stxs_bin in eqns["ggH"].keys():
  if stxs_bin.split("_")[0] == "QQ2HLL":
    eqns["ggZH"]["GG2HLL_%s"%stxs_bin[7:]] = eqns["ZH_lep"][stxs_bin]

bin_proc_map = od()
bin_proc_map["GG2H"] = "ggH"
bin_proc_map["QQ2HQQ"] = "qqH"
bin_proc_map["QQ2HLNU"] = "WH_lep"
bin_proc_map["QQ2HLL"] = "ZH_lep"
bin_proc_map["GG2HLL"] = "ggZH"
bin_proc_map["TTH"] = "ttH"
bin_proc_map["BBH"] = "bbH"
bin_proc_map["TH"] = "tH"

prod_eqns = od()
for stxs_bin in eqns["ggH"].keys():
  if stxs_bin != "UNKNOWN":
    bin_tag = stxs_bin.split("_")[0]
    print("%s -> %s"%(stxs_bin, bin_proc_map[bin_tag]))
    prod_eqns[stxs_bin] = eqns[bin_proc_map[bin_tag]][stxs_bin]

if CLEAN: prod_eqns = cleanUp(prod_eqns)

with open("prod.txt", "w") as f:
  f.write(jsonToText(prod_eqns))

with open("prod.json", "w") as f:
  json.dump(prod_eqns, f, indent=4)

"""
For decay, we need eqns for ZZ, gamgam, bb, WW, tautau.
Also need to sum up all of them for the total width equation
"""

"""
print(json.dumps(eqns["H_aa"], indent=4))
print(json.dumps(eqns["H_4l"], indent=4))
print(weights["H_aa"])
print(weights["H_4l"])

combined_eqn = weightedSumBin([eqns["H_aa"]["GG2H_PTH_300_450"], eqns["H_4l"]["GG2H_PTH_300_450"]], [weights["H_aa"], weights["H_4l"]])
print(json.dumps(combined_eqn, indent=4))
"""

decay_eqns = od()

b = "GG2H_PTH_300_450" #the bin the decay events land in
decay_eqns["ZZ"] = eqns["H_4l"][b]
decay_eqns["gamgam"] = eqns["H_aa"][b]
decay_eqns["bb"] = eqns["H_bb"][b]
decay_eqns["WW"] = eqns["H_lvlv"][b]
decay_eqns["tautau"] = eqns["H_tautau"][b]


decay_eqns["tot"] = weightedSumBin([eqns[mode][b] for mode in decay_modes], [weights[mode] for mode in decay_modes])

if CLEAN: decay_eqns = cleanUp(decay_eqns)

with open("decay.txt", "w") as f:
  f.write(jsonToText(decay_eqns))

with open("decay.json", "w") as f:
  json.dump(decay_eqns, f, indent=4)
