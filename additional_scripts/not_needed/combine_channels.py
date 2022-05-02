import json
from collections import OrderedDict

#helper function for reading json and converting keys to integers
def keysToInt(x):
  od = OrderedDict()
  for k, v in x:
    od[int(k)] = v
  return od

channels = ["llll", "llnunu", "llqq", "lnuqq", "nunuqq", "qqqq"]

combined_eqn = OrderedDict()

for channel in channels:
  with open("H_%s_SMEFT_combined.json"%channel, "r") as f:
    combined_eqn[channel] = json.load(f, object_pairs_hook=OrderedDict)

with open("decay_H_4j.json", "w") as f:
  json.dump(combined_eqn, f, indent=4)


