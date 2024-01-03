import json
import sys
from collections import OrderedDict
import numpy as np
import threading

def tqdm(iterable):
  for i, item in enumerate(iterable):
    j = i + 1
    sys.stdout.write('\r'+"%d/%d"%(j, len(iterable)))
    sys.stdout.flush()
    if j == len(iterable): print("")
    yield item

def getHistNames(loaded_file):
  return [name.split("_active_bins")[0] for name in loaded_file.keys() if "_active_bins" in name]

def getRwNames(loaded_file, hist_name):
  return list(filter(lambda x: (hist_name in x) and ("{}_".format(hist_name) not in x) and ("[rw" in x), loaded_file.keys()))

def getUniqueRwNames(loaded_files, hist_name):
  all_names = set()
  for loaded_file in loaded_files:
    all_names.update(getRwNames(loaded_file, hist_name))
  return all_names

def getUniqueBinEdges(loaded_files, hist_name):
  all_edges = set()
  for loaded_file in loaded_files:
    all_edges.update(getBinEdges(loaded_file, hist_name))
  return sorted(list(all_edges))

def getActiveBins(loaded_file, hist_name):
  return loaded_file["%s_active_bins"%hist_name]

def getBinEdges(loaded_file, hist_name):
  return [bin_info[0] for bin_info in getActiveBins(loaded_file, hist_name)]

def expandFile(loaded_file, hist_name, rw_names, bin_edges):
  """If rw_name does not exist in file, create it with sumw=0"""
  rw_names_not_exist = set(rw_names).difference(loaded_file.keys())
  for rw_name in rw_names_not_exist:
    loaded_file[rw_name] = [[0.,0.] for i in getActiveBins(loaded_file, hist_name)]

  active_bins = getActiveBins(loaded_file, hist_name)
  for i in range(len(bin_edges)):
    if (len(active_bins)==0) or (active_bins[i%len(active_bins)][0] != bin_edges[i]): #if this bin not in loaded_file
      active_bins.insert(i, [bin_edges[i], 0.0])
      for rw_name in rw_names:
        loaded_file[rw_name].insert(i, [0.,0.])
  assert getBinEdges(loaded_file, hist_name) == bin_edges, "\n%s\n%s"%(active_bins, bin_edges)

def combineActiveBins(file1, file2, hist_name):
  active_bins_1 = getActiveBins(file1, hist_name)
  active_bins_2 = getActiveBins(file2, hist_name)
  for i in range(len(active_bins_1)):
    active_bins_1[i][1] += active_bins_2[i][1]
  return active_bins_1

def loadJson(input_file, loaded_files):
  with open(input_file, "r") as f:
    loaded_files.append(json.load(f, object_pairs_hook=OrderedDict))

def merge(output_file, input_files):
  print(">> Loading files")
  loaded_files = []
  # for input_file in tqdm(input_files):
  #   with open(input_file, "r") as f:
  #     loaded_files.append(json.load(f, object_pairs_hook=OrderedDict))
  threads = []
  for input_file in input_files:
   thread = threading.Thread(target=loadJson, args=(input_file, loaded_files))
   thread.start()
   threads.append(thread)
  for thread in tqdm(threads):
    thread.join()  

  hist_names = getHistNames(loaded_files[0])
  #check that are have the same hist_names
  for i, loaded_file in enumerate(loaded_files):
    assert set(hist_names) == set(getHistNames(loaded_file)), "\n%s\n%s"%(hist_names, getHistNames(loaded_file))

  #find all unique names/keys
  all_names = set()
  for loaded_file in loaded_files:
    all_names.update(loaded_file.keys())

  combined = loaded_files[0]

  for hist_name in hist_names:
  #for hist_name in ["HTXS_stage0", "HTXS_stage1_2_pTjet30"]:
    print(">> Combining %s histograms"%hist_name)
    rw_names = getUniqueRwNames(loaded_files, hist_name)
    bin_edges = getUniqueBinEdges(loaded_files, hist_name)
    print("> Expanding files")
    for loaded_file in tqdm(loaded_files):
      expandFile(loaded_file, hist_name, rw_names, bin_edges)

    print("> Summing histograms")
    for loaded_file in tqdm(loaded_files[1:]):
      combined["%s_active_bins"%hist_name] = combineActiveBins(combined, loaded_file, hist_name)
      for rw_name in rw_names:
        combined[rw_name] = (np.array(combined[rw_name]) + np.array(loaded_file[rw_name])).tolist()

  with open(output_file, "w") as f:
    dumps = json.dumps(combined, sort_keys=True)
    dumps = dumps.replace(' "', '\n"')
    f.write(dumps)

def main(output_file, input_files, n_max_files=50):
  if len(input_files) <= n_max_files:
    merge(output_file, input_files)
    return True
  else:
    print(">> %d files left\n"%len(input_files))
    merge(output_file+".temp", input_files[:n_max_files])
    input_files = input_files[n_max_files:]
    input_files.insert(0, output_file+".temp")
    main(output_file, input_files, n_max_files)

if __name__=="__main__":
  output_file = sys.argv[1]
  input_files = sys.argv[2:]
  main(output_file, input_files)
