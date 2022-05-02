import yoda
import sys
from collections import OrderedDict
import json

def isHist(aos, name):
    return (type(aos[name]) == yoda.core.Histo1D) and ("RAW" in name) and ("AUX" not in name) and ("rw" in name)

def main(input_file, output_file):
  aos = yoda.read(input_file, asdict=True)

  converted = OrderedDict() 

  names = list(filter(lambda x: isHist(aos, x), aos.keys()))
  names = sorted(names)
  if len(names) == 0: exit(1)

  hist_names = set()
  for name in names:
    hist_names.update([name.split("/")[-1].split("[")[0]])

  for hist_name in hist_names:
    sm_name = list(filter(lambda x: "%s[rw0000]"%hist_name in x, names))
    print(sm_name)
    assert len(sm_name) == 1
    sm_name = sm_name[0]

    active_bins = []
    active_idx  = []
    total_n_events = 0
    for i, bin in enumerate(aos[sm_name]):
      if bin.numEntries > 0:
        active_bins.append([bin.xEdges[0], bin.numEntries])
        active_idx.append(i)
        total_n_events += bin.numEntries
    converted["%s_active_bins"%hist_name] = active_bins

    for name in filter(lambda x: hist_name in x, names):
      #bin_info = [[aos[name].bins[i].sumW*total_n_events, aos[name].bins[i].sumW2*total_n_events] for i in active_idx] #undo scaling by 1/n_events
      bin_info = [[aos[name].bins[i].sumW, aos[name].bins[i].sumW2] for i in active_idx]
      if bin_info != [[0., 0.] for bin in active_bins]:  converted[name.split("/")[-1]] = bin_info
      #converted[name.split("/")[-1]] = bin_info

    # #If rw0000 xs is zero (can happen for ggh), take sm xs from original xs
    # if sm_name.split("/")[-1] not in converted.keys(): 
    #   original_name = sm_name.split("[")[0]
    #   print(original_name)
    #   bin_info = [[aos[original_name].bins[i].sumW, aos[original_name].bins[i].sumW2] for i in active_idx] #no need to undo 1/n_events scaling
    #   converted[sm_name.split("/")[-1]] = bin_info

  with open(output_file, "w") as f:
    dumps = json.dumps(converted)
    dumps = dumps.replace(' "', '\n"')
    f.write(dumps)

if __name__=="__main__":
  input_file = sys.argv[1]
  output_file = sys.argv[2]
  main(input_file, output_file)