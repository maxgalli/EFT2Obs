import json
from collections import OrderedDict
import argparse
import numpy as np
import os

def getTerm(eqn, key):
  if key in eqn.keys(): return [eqn[key], eqn["u_"+key]]
  else: return [0,0]

one_to_one_conversion = OrderedDict()
one_to_one_conversion["cll1221"] = "cll1"
one_to_one_conversion["cpdc"] = "chdd"
one_to_one_conversion["cdp"] = "chbox"
one_to_one_conversion["c3pl1"] = "chl3"
one_to_one_conversion["c3pl2"] = "chl3"
one_to_one_conversion["cpq3i"] = "chj3"
one_to_one_conversion["cpq3"] = "chq3"
one_to_one_conversion["cpu"] = "chu"
one_to_one_conversion["cpt"] = "cht"
one_to_one_conversion["cpd"] = "chd"
one_to_one_conversion["ctp"] = "cthre"
one_to_one_conversion["ctg"] = "ctgre"
one_to_one_conversion["cpg"] = "chg"

#print(one_to_one_conversion)

# def rotationConversion(eqns):
#   #rotation conversion: cpqmi = chj1-chj3, cpqm = chq1-chq3
#   for stxs_bin in eqns.keys():
#     for key in eqns[stxs_bin].keys():
#       if ("cpqmi" in key) or ("cpqm" in key):
#         value = eqns[stxs_bin][key]
#         del eqns[stxs_bin][key]

#         if "cpqmi" in key:
#           eqns[stxs_bin][key.replace("cpqmi", "chj1")] = value
#           if key[0] == "A":
#             eqns[stxs_bin][key.replace("cpqmi", "chj3")] = -value
#           else:
#             eqns[stxs_bin][key.replace("cpqmi", "chj3")] = value

#           if key[-1] == "B": #if key == B_c_2
#             eqns[stxs_bin]["B"]
#         else:
#           eqns[stxs_bin][key.replace("cpqm", "chq1")] = value
#           if key[0] == "A":
#             eqns[stxs_bin][key.replace("cpqm", "chq3")] = -value
#           else:
#             eqns[stxs_bin][key.replace("cpqm", "chq3")] = value
#   return eqns
#   #print(json.dumps(eqns["GG2HLL"], indent=4))

def rotationConversion(eqns):
  #rotation conversion: cpqmi = chj1-chj3, cpqm = chq1-chq3
  for stxs_bin in eqns.keys():
    eqns[stxs_bin]["A_chj1"] = eqns[stxs_bin]["A_cpqmi"]
    eqns[stxs_bin]["u_A_chj1"] = eqns[stxs_bin]["u_A_cpqmi"]
    eqns[stxs_bin]["A_chj3"] = -eqns[stxs_bin]["A_cpqmi"]
    eqns[stxs_bin]["u_A_chj3"] = eqns[stxs_bin]["u_A_cpqmi"]

    eqns[stxs_bin]["A_chq1"] = eqns[stxs_bin]["A_cpqm"]
    eqns[stxs_bin]["u_A_chq1"] = eqns[stxs_bin]["u_A_cpqm"]
    eqns[stxs_bin]["A_chq3"] = -eqns[stxs_bin]["A_cpqm"]
    eqns[stxs_bin]["u_A_chq3"] = eqns[stxs_bin]["u_A_cpqm"]


    eqns[stxs_bin]["B_chj1_2"] = eqns[stxs_bin]["B_cpqmi_2"]
    eqns[stxs_bin]["u_B_chj1_2"] = eqns[stxs_bin]["u_B_cpqmi_2"]
    eqns[stxs_bin]["B_chj3_2"] = eqns[stxs_bin]["B_cpqmi_2"]
    eqns[stxs_bin]["u_B_chj3_2"] = eqns[stxs_bin]["u_B_cpqmi_2"]

    eqns[stxs_bin]["B_chq1_2"] = eqns[stxs_bin]["B_cpqm_2"]
    eqns[stxs_bin]["u_B_chq1_2"] = eqns[stxs_bin]["u_B_cpqm_2"]
    eqns[stxs_bin]["B_chq3_2"] = eqns[stxs_bin]["B_cpqm_2"]
    eqns[stxs_bin]["u_B_chq3_2"] = eqns[stxs_bin]["u_B_cpqm_2"]


    eqns[stxs_bin]["B_chj1_chj3"] = -2*eqns[stxs_bin]["B_cpqmi_2"]
    eqns[stxs_bin]["u_B_chj1_chj3"] = 2*eqns[stxs_bin]["u_B_cpqmi_2"]
    eqns[stxs_bin]["B_chq3_chq1"] = -2*eqns[stxs_bin]["B_cpqm_2"]
    eqns[stxs_bin]["u_B_chq3_chq1"] = 2*eqns[stxs_bin]["u_B_cpqm_2"]


    for param in filter(lambda x: x[:2] != "u_", eqns[stxs_bin].keys()):
      if (param[0] == "B") and (param[-2:] != "_2"): #if cross term
        if "cpqmi" in param:
          eqns[stxs_bin][param.replace("cpqmi", "chj1")] = eqns[stxs_bin][param]
          eqns[stxs_bin]["u_"+param.replace("cpqmi", "chj1")] = eqns[stxs_bin]["u_"+param]
          eqns[stxs_bin][param.replace("cpqmi", "chj3")] = -eqns[stxs_bin][param]
          eqns[stxs_bin]["u_"+param.replace("cpqmi", "chj3")] = eqns[stxs_bin]["u_"+param]
        elif "cpqm" in param:
          eqns[stxs_bin][param.replace("cpqm", "chq1")] = eqns[stxs_bin][param]
          eqns[stxs_bin]["u_"+param.replace("cpqm", "chq1")] = eqns[stxs_bin]["u_"+param]
          eqns[stxs_bin][param.replace("cpqm", "chq3")] = -eqns[stxs_bin][param]
          eqns[stxs_bin]["u_"+param.replace("cpqm", "chq3")] = eqns[stxs_bin]["u_"+param]

    for key in eqns[stxs_bin].keys():
      if ("cpqmi" in key) or ("cpqm" in key):
        del eqns[stxs_bin][key]

  return eqns
  #print(json.dumps(eqns["GG2HLL"], indent=4))

def convert_SMEFTatNLO_To_SMEFTsim(eqns):
  eqns = rotationConversion(eqns)

  a_s = 0.1181
  g_s = 2*np.sqrt(a_s*np.pi)

  for stxs_bin in eqns.keys():
    print(stxs_bin)
    params = filter(lambda x: x[:2] != "u_", eqns[stxs_bin].keys())
   
    converted_eqn = OrderedDict() 
    for param in sorted(params):
      term = getTerm(eqns[stxs_bin], param)
      for key, value in one_to_one_conversion.items():
        param = param.replace(key, value)

      #print(param, param[-2:], param.split("_"))
      #if term is like B_chl3_chl3
      if (param[0] == "B") and (param[-2:] != "_2") and (param.split("_")[1] == param.split("_")[2]):
        param = "B_%s_2"%(param.split("_")[1])
   
      if "ctg" in param:
        if "_2" not in param:
           term[0] /= -g_s
           term[1] /= g_s
        else:
           term[0] /= (-g_s)**2
           term[1] /= g_s**2

      #if chl3 term already exsits, add contribution on top
      if param in converted_eqn.keys():
        #print(param)
        #if "chj3" in param: print(converted_eqn[param], term[0])
        converted_eqn[param] += term[0]
        converted_eqn["u_"+param] = np.sqrt(converted_eqn["u_"+param]**2 + term[1]**2)
      else:
        converted_eqn[param] = term[0]
        converted_eqn["u_"+param] = term[1]
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

#print(converted_eqns.keys())

with open(args.output, "w") as f:
  json.dump(converted_eqns, f, indent=4)

