import json
import os

with open("widths_cross_sections_xswg.json", "r") as f:
  widths = json.load(f)

channels = ["aa", "bb", "cc", "gg", "llll", "llnunu", "llqq", "lnuqq", "nunuqq", "qqqq", "tautau", "Za"]
channels = ["H_%s"%channel for channel in channels]

def getEquationJsons(channels):
  jsons = []
  for channel in channels:
    if os.path.exists("ConvertedEquations/%s_SMEFT_topU3l.json"%channel):
      jsons.append("ConvertedEquations/%s_SMEFT_topU3l.json"%channel)
    else:
      jsons.append("ConvertedEquations/%s_SMEFT.json"%channel)

  return jsons

equation_jsons = getEquationJsons(channels)

import combine_and_expand

function_str = "("
for i, channel in enumerate(channels):
  function_str += "%.10f*eq%d +"%(widths[channel], i)
function_str = function_str[:-1]
function_str += ")"

total_sm_width = sum([widths[channel] for channel in channels])
function_str += "/%.10f"%total_sm_width
print(equation_jsons)
print(function_str)

combine_and_expand.main(equation_jsons, "ConvertedEquations/H_total_SMEFT_topU3l.json", function_str, False)

