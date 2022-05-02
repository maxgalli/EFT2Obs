import json
import os
import argparse
from collections import OrderedDict
import combine_and_expand
import numpy as np

def getPath(proc, input_dir):
  H_4f_channels = ["llll", "llnunu", "llqq", "lnuqq", "nunuqq", "qqqq"]
  analytic_channels = ["aa", "Za"]
  if proc in H_4f_channels:
    path = os.path.join(input_dir, "H_%s_SMEFTsim_topU3l_inclusive_combined.json"%proc)
  elif proc in analytic_channels:
    path = os.path.join(input_dir, "H_%s_SMEFTsim_topU3l_inclusive_analytic.json"%proc)
  else:
    path = os.path.join(input_dir, "H_%s_SMEFTsim_topU3l_inclusive.json"%proc)
  return path

def loadEqn(proc, input_dir):  
  path = getPath(proc, input_dir)

  with open(path, "r") as f:
    eqn = json.load(f, object_pairs_hook=OrderedDict)
    assert len(eqn.keys()) == 1, eqn.keys()
  
  return eqn.items()[0][1]

# def produceTotalWidth(channels, input_dir):
#   with open("widths_cross_sections_xswg.json", "r") as f:
#     widths = json.load(f)

#   equation_jsons = [getPath(channel, input_dir) for channel in channels]

#   function_str = "("
#   for i, channel in enumerate(channels):
#     function_str += "%.10f*eq%d +"%(widths["H_%s"%channel], i)
#   function_str = function_str[:-1]
#   function_str += ")"

#   total_sm_width = sum([widths["H_%s"%channel] for channel in channels])
#   function_str += "/%.10f"%total_sm_width
#   print(equation_jsons)
#   print(function_str)

#   combine_and_expand.main(equation_jsons, "%s/H_total_SMEFTsim_topU3l_inclusive.json"%input_dir, function_str, False)

def produceTotalWidth(channels, input_dir):
  with open("widths_cross_sections_xswg.json", "r") as f:
    widths = json.load(f)
  tot_width_sm = sum([widths["H_%s"%channel] for channel in channels])

  eqns = {channel: loadEqn(channel, input_dir) for channel in channels}

  terms = set()
  for channel in channels:
    terms.update(filter(lambda x: x[0] != "u", eqns[channel].keys()))

  tot = OrderedDict()
  for term in terms:
    sum_term = 0
    sum_error_2 = 0
    for channel in channels:
      #print(channel)
      if term in eqns[channel].keys():
        sum_term += widths["H_%s"%channel] * eqns[channel][term]
        if "u_"+term in eqns[channel].keys(): sum_error_2 += (widths["H_%s"%channel] * eqns[channel]["u_"+term])**2
    tot[term] = sum_term / tot_width_sm
    tot["u_"+term] = np.sqrt(sum_error_2) / tot_width_sm

  with open(os.path.join(input_dir, "H_total_SMEFTsim_topU3l_inclusive.json"), "w") as f:
    json.dump({"inclusive": tot}, f, sort_keys=True, indent=4)
    #json.dump({"inclusive": tot}, f, sort_keys=True, indent=4)

def main(input_dir, output):
  total_width_channels =  ["aa", "bb", "cc", "gg", "llll", "llnunu", "llqq", "lnuqq", "nunuqq", "qqqq", "tautau", "Za"]
  all_channels = total_width_channels + ["mumu"] + ["Zgam"]

  produceTotalWidth(total_width_channels, input_dir)

  combined_eqn = OrderedDict()
  combined_eqn["ZZ"] = loadEqn("llll", input_dir)
  combined_eqn["gamgam"] = loadEqn("aa", input_dir)
  combined_eqn["bb"] = loadEqn("bb", input_dir)
  combined_eqn["WW"] = loadEqn("llnunu", input_dir)
  combined_eqn["tautau"] = loadEqn("tautau", input_dir)
  combined_eqn["mumu"] = loadEqn("mumu", input_dir)
  combined_eqn["Zgam"] = loadEqn("Za", input_dir)
  combined_eqn["cc"] = loadEqn("cc", input_dir)
  combined_eqn["gluglu"] = loadEqn("gg", input_dir)

  combined_eqn["tot"] = loadEqn("total", input_dir)

  with open(output, "w") as f:
    json.dump(combined_eqn, f, indent=4)

if __name__=="__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('--input-dir', '-i', type=str, default="ConvertedEquations/")
  parser.add_argument('--output', '-o', type=str, default="decay.json")
  
  args = parser.parse_args()

  main(args.input_dir, args.output)