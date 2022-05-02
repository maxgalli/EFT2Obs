"""
Multiply equations together and expand to linear or quadratic order.
"""

import json
import argparse
from collections import OrderedDict as od
import ast

def getEqnVal(eqn, bin_name, param_vals):
  non_empty_names = getNonEmptyBinNames(eqn)
  if len(non_empty_names) == 1:
    bin_eqn = eqn[non_empty_names[0]]
  else:
    bin_eqn = eqn[bin_name]

  sum = 1
  for term_name, coeff in bin_eqn.items():
    if term_name[0] == "u": continue
    params = term_name.split("_")[1:] #term name can be A_chw, B_chw_2, B_chw_chb
    if params[-1] == "2": 
      params[-1] = params[0]

    whole_term = coeff #calculate coeff * param1 * param2
    for param in params:
      if param in param_vals.keys():
        whole_term *= param_vals[param]
      else:
        #if param not given in param_vals, assume it is set to be zero
        whole_term = 0
    sum += whole_term

  return sum

def getNonEmptyBinNames(eqn):
  names = eqn.keys()
  #print(names)
  non_empty_names = []
  for name in names:
    if len(eqn[name]) > 0:
      non_empty_names.append(name)
  return non_empty_names

def assertBinNames(eqns):
  #make sure each eqn has either all the same bin names or just 1 bin
  bin_names = [getNonEmptyBinNames(eqn) for eqn in eqns]
  all_one_bin = True
  for names in bin_names:
    if len(names) > 1:
      all_one_bin = False
      names_to_check_list = names
      names_to_check = set(names)
      break
  if all_one_bin: return "all_one_bin"

  for names in bin_names:
    if len(names) == 1: continue
    if len(names_to_check.union(names)) != len(names_to_check):
      return False
  
  return names_to_check_list

def findAllParams(eqns):
  params = []

  for eqn in eqns:
    for bin_name, bin_eqn in eqn.items():
      for term in bin_eqn.keys():
        if term[0] == "u": continue
        new_params = term.split("_")[1:]
        if "2" in new_params: new_params.remove("2")
        
        for param in new_params:
          if param not in params:
            params.append(param)

  return params

def createGetVal(eqns, function_str):
  assert function_str.count("eq") == len(eqns)
  for i in range(len(eqns)):
    function_str = function_str.replace("eq%d "%i, "getEqnVal(eqns[%d], bin_name, param_vals)"%i)
  function_str = "lambda eqns, bin_name, param_vals: " + function_str

  getVal = eval(function_str)
  return getVal

def deriveEquation(eqns, function_str, linear_only, x=1e-4):
  result = assertBinNames(eqns)
  assert result != False

  if result == "all_one_bin":
    bin_names = ["inclusive"]
  else:
    bin_names = result

  #print(bin_names)

  all_params = findAllParams(eqns)
  new_eqns = {}

  getVal = createGetVal(eqns, function_str)

  for bin_name in bin_names:
    #get linear and quadratic (not cross)
    new_eqn = od()
    for param in all_params:
      mu1 = getVal(eqns, bin_name, {param: x/2})
      mu2 = getVal(eqns, bin_name, {param: x})

      new_eqn["A_%s"%param] = (4*mu1-mu2-3) / x
      if not linear_only:  new_eqn["B_%s_2"%param] = (2*(mu2-2*mu1+1)) / x**2
    
    if not linear_only:
      n_params = len(all_params)
      for i in range(n_params):
        for j in range(i+1, n_params):
          p1, p2 = all_params[i], all_params[j]

          mu1 = getVal(eqns, bin_name, {p1: x})
          mu2 = getVal(eqns, bin_name, {p2: x})
          mu3 = getVal(eqns, bin_name, {p1: x, p2: x})
          new_eqn["B_%s_%s"%(p1,p2)] = ((mu3-1) - (mu1-1 + mu2-1)) / x**2

    new_eqns[bin_name] = new_eqn

  return new_eqns

def main(equations_jsons, output, function_str, linear_only):
  eqns = []
  for eqn_path in equations_jsons:
    with open(eqn_path, "r") as f:
      eqns.append(json.load(f))

  for i in range(len(eqns)):
    item = eqns[i].items()[0][1]
    if type(item) == float:
      eqns[i] = {"inclusive": eqns[i]}

  new_eqns = deriveEquation(eqns, function_str, linear_only)

  with open(output, "w") as f:
    json.dump(new_eqns, f, indent=4)
  
  return eqns

if __name__=="__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('--equations-jsons', '-e', type=str, nargs="+", required=True)
  parser.add_argument('--output', '-o', type=str, required=True)
  parser.add_argument('--function', type=str, default=None)
  parser.add_argument('--linear-only', action="store_true")

  args = parser.parse_args()

  if args.function == None:
    args.function = "".join(["eq%d*"%i for i in range(len(args.equations_jsons))])[0:-1]
  print(args.function)

  eqns = main(args.equations_jsons, args.output, args.function, args.linear_only)

  #x = 2
  #print(getVal(eqns, "GG2HLL_PTV_0_75", {"chl3":x}))
  #print(1 + 0.009761624020791556*x +0.002957358458270093*x**2)
  #print(1 -0.24250108297146233*x + 0.015003116844696606*x**2)

  print(-0.24250108297146233 + 0.009761624020791556)
  print(0.015003116844696606 + 0.002957358458270093 + (-0.24250108297146233*0.009761624020791556))

  print(-0.24250108297146233 - 0.009761624020791556)
  print(0.015003116844696606 - 0.002957358458270093 - (-0.24250108297146233*0.009761624020791556))

  