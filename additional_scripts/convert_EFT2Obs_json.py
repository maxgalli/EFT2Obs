import json
from collections import OrderedDict
import sys
import warnings

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

from optparse import OptionParser
parser = OptionParser(usage="%prog input.json output.json")
parser.add_option("--key", dest="key", 
                    help="Key to interpret bin numbers.")
parser.add_option("--relative-threshold", dest="relative_threshold", default=1e-16, type=float,
                  help="Remove terms where the coefficients are, as a fraction to the biggest coefficient, smaller than the threshold")

(options, args) = parser.parse_args()

if len(args) < 2:
  parser.print_help()
  sys.exit(1)

input_json = args[0]
output_json = args[1]

with open(input_json, "r") as f:
  eft2obs_json = json.loads(f.read())

if options.key is not None:
  with open(options.key, "r") as f:
    key = json.loads(f.read(), object_pairs_hook=keysToInt)
  key = {i:key[tag] for i, tag in enumerate(key.keys())}
  options.key = key

new_json_dict = jsonToNewDict(eft2obs_json, options.key)
#new_json_dict = splitTH(new_json_dict, input_json)

with open(output_json, "w") as f:
  f.write(json.dumps(new_json_dict, indent=4))

