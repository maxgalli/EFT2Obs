import json
import os
import argparse
from collections import OrderedDict

def splitTH(converted_eqn, input_json):
  if "TH" in converted_eqn.keys():
    converted_eqn["THQ_FWDH"] = {}
    converted_eqn["THQ"] = {}
    converted_eqn["THW_FWDH"] = {}
    converted_eqn["THW"] = {}

    if "tHW" in input_json:
      converted_eqn["THW_FWDH"] = converted_eqn["TH_FWDH"]
      converted_eqn["THW"] = converted_eqn["TH"]
    elif "tHq" in input_json:
      converted_eqn["THQ_FWDH"] = converted_eqn["TH_FWDH"]
      converted_eqn["THQ"] = converted_eqn["TH"]

    del converted_eqn["TH_FWDH"]
    del converted_eqn["TH"]
  return converted_eqn

def getEqnPath(input_dir, proc, binning):
  if ("ggH" in proc) or ("ggZH" in proc):
    print(proc)
    file_name = "%s_%s_combined.json"%(proc, binning)
  else:
    file_name = "%s_%s.json"%(proc, binning)
  path = os.path.join(input_dir, file_name)
  return path

def loadEqns(procs, input_dir):
  eqns = []
  for binning in ["stage0", "stage1p2"]:
    for proc in procs:
      with open(getEqnPath(input_dir, proc, binning)) as f:
        eqn = json.load(f, object_pairs_hook=OrderedDict)
      eqn = splitTH(eqn, proc)
      eqns.append(eqn)
  return eqns

def getBinNames(eqns):
  #bin_names = eqns[0].keys()
  #for eqn in eqns:
  #  if set(bin_names).symmetric_difference(eqn.keys()) != set():
  #    raise Exception("Bin names should be same across equations")
  bin_names = set()
  for eqn in eqns:
    bin_names.update(eqn.keys())
  return sorted(bin_names)

def combineEqns(eqns):
  bin_names = getBinNames(eqns)
  print(bin_names)
  combined_eqn = OrderedDict()

  for bin_name in bin_names:
    print(bin_name)
    already_claimed = False
    for eqn in eqns:
      if bin_name in eqn.keys():
        if eqn[bin_name] != {}:
          if already_claimed and ("FWDH" not in bin_name) and (bin_name.count("_") > 1):
            print([key for key in eqn.keys() if eqn[key] != {}])
            raise Exception("Only one file should contain equation from same bin")
          else:
            print([key for key in eqn.keys() if eqn[key] != {}])
            combined_eqn[bin_name] = eqn[bin_name]
            already_claimed = True

  return combined_eqn

def main(input_dir, output):
  procs = ["ggH_SMEFTatNLO", "ggZH_SMEFTatNLO", "bbH_SMEFTsim_topU3l", "qqH_SMEFTsim_topU3l", "tHq_SMEFTsim_topU3l", "tHW_SMEFTsim_topU3l", "ttH_SMEFTsim_topU3l", "WH_lep_SMEFTsim_topU3l", "ZH_lep_SMEFTsim_topU3l"]

  eqns = loadEqns(procs, input_dir)
  combined_eqn = combineEqns(eqns)

  with open(output, "w") as f:
    json.dump(combined_eqn, f, indent=4)

if __name__=="__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('--input-dir', '-i', type=str, default="ConvertedEquations/")
  parser.add_argument('--output', '-o', type=str, default="prod.json")
  
  args = parser.parse_args()

  main(args.input_dir, args.output)