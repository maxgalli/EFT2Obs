import sys
import os
import subprocess
import numpy as np
import tools
import argparse
from collections import OrderedDict as od
import json

def skimResults(args):
  #results_dir = os.path.join("directResults", args.process)

  try:
    grep_output = subprocess.check_output(["grep 'Width ' directResults/%s_*/*/*"%args.process], shell=True).split("\n")[:-1]
  except:
    grep_output = subprocess.check_output(["grep 'Cross-section ' directResults/%s_*/*/*"%args.process], shell=True).replace("Cross-section", "Width").replace("pb", "GeV").split("\n")[:-1]
  print("\n".join(grep_output))

  results = {}
  for result in grep_output:
    rw_num = int(result.split("rw_")[1].split("_")[0])
    width_statement = result.split("Width :")[1]

    width = float(width_statement.split("+-")[0])
    error = float(width_statement.split("+-")[1].split("GeV")[0])
    results[rw_num] = [width, error]

  cfg = tools.GetConfigFile(args.config)
  n_par = len(cfg['parameters'])
  max_rw = 2*n_par + n_par*(n_par-1)/2

  for i in range(max_rw+1):
    if i not in results.keys():
      results[i] = [0, -1]

  for i in range(max_rw+1):
    print(i, results[i])

  return results 

def getQuadXs(eqn, nominal_eqn, param, x, sm):
  if eqn["B_%s_2"%param] == 0:
    if "B_%s_2"%param in nominal_eqn.keys():
      xs = nominal_eqn["B_%s_2"%param]*x*x*sm[0]
      error = abs(np.sqrt( (nominal_eqn["u_B_%s_2"%param]/nominal_eqn["B_%s_2"%param])**2 + (sm[1]/sm[0])**2 ) * xs)
      return [xs, error]
    else:
      return [0, 0]
  else:
    xs = eqn["B_%s_2"%param]*x*x*sm[0]
    error = abs(np.sqrt( (eqn["u_B_%s_2"%param]/eqn["B_%s_2"%param])**2 + (sm[1]/sm[0])**2 ) * xs)
    return [xs, error]

def deriveEquations(results, args):
  cfg = tools.GetConfigFile(args.config)
  n_par = len(cfg['parameters'])
  x = cfg['parameter_defaults']['val']

  eqn = od()

  with open(args.nominal_eqn, "r") as f:
    nominal_eqn = getNonEmptyBin(json.load(f))

  sm = results[0]
  for i in range(n_par):
    rw1 = results[i*2+1]
    rw2 = results[i*2+2]

    A_err = 0
    B_err = 0

    A = (2*rw1[0])/(x*sm[0])
    if rw1[0] != 0: A_err = abs(np.sqrt(((2*rw1[1])/(2*rw1[0]))**2 + ((x*sm[1])/(x*sm[0]))**2) * A)
 
    B = rw2[0]/(x*x*sm[0])
    if rw2[0] != 0: B_err = abs(np.sqrt((rw2[1]/rw2[0])**2 + ((x*x*sm[1])/(x*x*sm[0]))**2) * B)

    param = cfg['parameters'][i]['name']
    eqn["A_%s"%param] = A
    eqn["u_A_%s"%param] = A_err
    eqn["B_%s_2"%param] = B
    eqn["u_B_%s_2"%param] = B_err

  #now cross terms
  k = n_par*2+1
  for i in range(n_par):
    for j in range(i+1, n_par):
      rw3 = results[k]
      pi = cfg['parameters'][i]['name']
      pj = cfg['parameters'][j]['name']

      print(pi, pj, rw3)

      if rw3[0] == 0: 
        eqn["B_%s_%s"%(pi,pj)] = 0
        eqn["u_B_%s_%s"%(pi,pj)] = 0
        k += 1
        continue
      
      xs_i = getQuadXs(eqn, nominal_eqn, pi, x, sm)
      xs_j = getQuadXs(eqn, nominal_eqn, pj, x, sm)

      new_xs = rw3[0] - xs_i[0] - xs_j[0]
      new_error = np.sqrt(rw3[1]**2 + xs_i[1]**2 + xs_j[1]**2)
      #new_xs = rw3[0] - nominal_eqn["B_%s_2"%pi]*x*x*sm[0] - nominal_eqn["B_%s_2"%pj]*x*x*sm[0]
      #new_error = np.sqrt(rw3[1]**2 + (nominal_eqn["u_B_%s_2"%pi]*x*x*sm[0])**2 + (nominal_eqn["u_B_%s_2"%pj]*x*x*sm[0])**2)
      #new_xs = rw3[0] - quad_xs[pi][0] - quad_xs[pj][0]
      #new_error = np.sqrt(rw3[1]**2 + quad_xs[pi][1]**2 + quad_xs[pj][1]**2)
      rw3 = [new_xs, new_error]

      B_ij_err = 0

      B_ij = rw3[0]/(x*x*sm[0])
      if rw3[0] != 0: B_ij_err = abs(np.sqrt((rw3[1]/rw3[0])**2 + ((x*x*sm[1])/(x*x*sm[0]))**2) * B_ij)
      
      eqn["B_%s_%s"%(pi,pj)] = B_ij
      eqn["u_B_%s_%s"%(pi,pj)] = B_ij_err  
      
      k += 1 
  return eqn

def genToyResults(n_par, nevents=10000, sm_width=1, x=0.01, smear=True):
   err = 1/np.sqrt(nevents)
   results = {0: [sm_width, sm_width*err]}
   actual_coeff = []

   for i in range(n_par):
     A = np.random.uniform(-0.1, 0.1)
     B = np.random.uniform(0, 10)
     actual_coeff.append([A, B])
     
     rw1 = (1+A*(x/2)+B*(x/2)**2)*sm_width
     if smear: rw1 = rw1 + np.random.normal(scale=rw1*err)
     rw2 = (1+A*(x)+B*(x)**2)*sm_width
     if smear: rw2 = rw2 + np.random.normal(scale=rw2*err)
     results[i*2+1] = [rw1, rw1*err]
     results[i*2+2] = [rw2, rw2*err]

   return results, actual_coeff

def getNonEmptyBin(eqns):
  for key in eqns.keys():
    if eqns[key] != od():
      return eqns[key]

def getTerm(eqn, key):
  if key in eqn.keys(): return [eqn[key], eqn["u_"+key]]
  else: return [0,0]

def combineEqns(direct_eqn, args):
  with open(args.nominal_eqn, "r") as f:
    nominal_eqn = getNonEmptyBin(json.load(f))

  combined_eqn = od()

  keys = set(direct_eqn.keys()).union(nominal_eqn.keys())

  for key in sorted(filter(lambda x: x[0:2] != "u_", keys)):
    direct_term = getTerm(direct_eqn, key)
    nominal_term = getTerm(nominal_eqn, key)

    if direct_term[0] == 0:
      combined_eqn[key], combined_eqn["u_"+key] = nominal_term[0], nominal_term[1]
    else:
      combined_eqn[key], combined_eqn["u_"+key] = direct_term        

  return combined_eqn

parser = argparse.ArgumentParser()
parser.add_argument('--process', '-p', default=None, help="Label of the process, must correspond to the dir name that was created in the MG dir")
parser.add_argument('--output', '-o', default='combined_eqn.json')
parser.add_argument('--config', '-c')
parser.add_argument('--nominal-eqn')
args = parser.parse_args()

if args.process == None: #if no direct results
  eqn = od()
else:
  results = skimResults(args)
  eqn = deriveEquations(results, args)
  for each in eqn:
    if eqn[each] != 0:
      print(each, eqn[each])

combined_eqn = combineEqns(eqn, args)

with open(args.output, "w") as f:
  json.dump({"inclusive":combined_eqn}, f, indent=4)