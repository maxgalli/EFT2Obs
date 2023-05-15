import json
from collections import OrderedDict
import argparse
import numpy as np
import os

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

one_to_one_conversion = {
  "cll1221": "cll1",
  "cpdc": "chdd",
  "cdp": "chbox",
  "c3pl1": "chl3",
  "c3pl2": "chl3",
  "ctp": "cthre",
  "ctg": "ctgre",
  "cpg": "chg"
}

def convert_SMEFTatNLO_To_SMEFTsim(eqns):
  a_s = 0.1181
  g_s = 2*np.sqrt(a_s*np.pi)

  for stxs_bin in eqns.keys():
    params = filter(lambda x: x[:2] != "u_", eqns[stxs_bin].keys())
   
    converted_eqn = OrderedDict() 
    for param in sorted(params):
      term = getTerm(eqns[stxs_bin], param)
      for key, value in one_to_one_conversion.items():
        param = param.replace(key, value)

      if param == "B_chl3_chl3":
        param = "B_chl3_2"
   
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
        converted_eqn[param] += term[0]
        converted_eqn["u_"+param] = np.sqrt(converted_eqn["u_"+param]**2 + term[1]**2)
      else:
        converted_eqn[param] = term[0]
        converted_eqn["u_"+param] = term[1]
    eqns[stxs_bin] = converted_eqn
  return eqns

# def addOtherParams(eqns):
#   v = 0.24622
#   params = OrderedDict()
#   params["chbox"] = 1
#   params["chdd"] = -0.25
#   params["chl3"] = -1
#   params["cll1"] = 0.5

#   print(params.keys())
#   n_par = len(params)

#   stxs_bins = eqns.keys()

#   for stxs_bin in stxs_bins:
#     if eqns[stxs_bin].keys() == []: continue
#     original_params = filter(lambda x: x[:2] == "A_", eqns[stxs_bin].keys())
#     original_params = [param[2:] for param in original_params]

#     for param in params.keys():
#       eqns[stxs_bin]["A_%s"%param] = params[param]*2*v**2
#       eqns[stxs_bin]["u_A_%s"%param] = 0

#     for i, p1 in enumerate(params.keys()):
#       for j, p2 in enumerate(params.keys()):
#         if p1 == p2:
#           eqns[stxs_bin]["B_%s_2"%p1] = params[p1]*params[p1]*v**4
#           eqns[stxs_bin]["u_B_%s_2"%p1] = 0
#         elif j > i:
#           eqns[stxs_bin]["B_%s_%s"%(p1, p2)] = params[p1]*params[p2]*2*v**4
#           eqns[stxs_bin]["u_B_%s_%s"%(p1, p2)] = 0
    
#       for p3 in original_params:
#         eqns[stxs_bin]["B_%s_%s"%(p1, p3)] = (params[p1]*2*v**2 * eqns[stxs_bin]["A_%s"%p3]) / 2
#         eqns[stxs_bin]["u_B_%s_%s"%(p1, p3)] = (params[p1]*2*v**2 * eqns[stxs_bin]["u_A_%s"%p3]) / 2
#   return eqns

def contractEquations(loop, tree, tree_loop_2, tree_loop_4):
  """Make sure all equations have the same bins. Sometimes a bin is occupied in one equation but not another."""
  common_bin_names = set(loop.keys())
  for each in [tree, tree_loop_2, tree_loop_4]:
    common_bin_names = common_bin_names.intersection(each.keys())
  
  for each in [loop, tree, tree_loop_2, tree_loop_4]:
    for bin_name in each.keys():
      if bin_name not in common_bin_names:
        print("Deleting %s"%bin_name)
        del each[bin_name]

def main(args):
  input_dir = args.input_dir
  obs = args.observable
  output = "ConvertedEquations/ggH_SMEFTatNLO_{}_combined_{}.json"
  if obs == "pt_h":
    decays = ["Hgg", "HZZ", "HWW", "Htt", "Hbb", "HbbVBF", "HttBoosted"]
  elif obs == "njets":
    decays = ["Hgg", "HZZ", "HWW", "Htt"]
  elif obs == "pt_j0":
    decays = ["HZZ", "Htt"]
  elif obs == "deta_jj":
    decays = ["HZZ"]
  elif obs == "deltaphijj":
      decays = ["Hgg", "HZZ"]
  for dec in decays:
    with open(os.path.join(input_dir, "ggH_SMEFTatNLO_{}_loop_{}.json".format(dec, obs)), "r") as f:
      loop = json.load(f)
    with open(os.path.join(input_dir, "ggH_SMEFTatNLO_{}_tree_{}.json".format(dec, obs)), "r") as f:
      tree = json.load(f)
    with open(os.path.join(input_dir, "ggH_SMEFTatNLO_{}_tree_loop_2_{}.json".format(dec, obs)), "r") as f:
      tree_loop_2 = json.load(f)
    with open(os.path.join(input_dir, "ggH_SMEFTatNLO_{}_tree_loop_4_{}.json".format(dec, obs)), "r") as f:
      tree_loop_4 = json.load(f)

    contractEquations(loop, tree, tree_loop_2, tree_loop_4)
    combined_eqn = combineEqns(loop, tree, tree_loop_2, tree_loop_4)
    combined_eqn = convert_SMEFTatNLO_To_SMEFTsim(combined_eqn)
    #combined_eqn = addOtherParams(combined_eqn)

    with open(output.format(dec, obs), "w") as f:
      json.dump(combined_eqn, f, indent=4)

if __name__=="__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('--input-dir', '-i', type=str, default="ConvertedEquations/")
  parser.add_argument('--observable', type=str, default="pt_h")
  #parser.add_argument('--output', '-o', default="ConvertedEquations/ggH_SMEFTatNLO_combined.json")
  #parser.add_argument('--postfix', type=str, default="")
  args = parser.parse_args()

  main(args)



