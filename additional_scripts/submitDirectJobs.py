import os
import sys
import argparse
import tools
from collections import OrderedDict as od
import importlib

"""SCRIPT IS SMEFTsim SPECIFIC!!!"""

CONDOR_TEMPLATE="""executable            = additional_scripts/run_rw_point.sh
arguments             = $(MG_process) $(rw_num) $(nevents) $(ncores)

output                = directResults/%(process)s/condor_$(ClusterId)_rw_$(rw_num)_nevents_$(nevents).out
error                 = directResults/%(process)s/condor_$(ClusterId)_rw_$(rw_num)_nevents_$(nevents).err
log                   = directResults/%(process)s/condor_$(ClusterId)_nevents_$(nevents).log

# Send the job to Held state on failure.
on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)
# Periodically retry the jobs every 10 minutes, up to a maximum of 5 retries.
periodic_release =  (NumJobStarts < 3) && ((CurrentTime - EnteredCurrentStatus) > 600)

#MG_process=%(process)s
nevents=%(nevents)d
ncores=%(ncores)d
request_cpus=$(ncores)

+JobFlavour = "longlunch"
#+MaxRuntime = 120

queue rw_num, MG_process from (
  %(rw_nums)s
)
"""

MG_DIR = os.getenv("MG_DIR")
TMPDIR = os.getenv("TMPDIR")
sys.path.append(os.path.join(MG_DIR))
sys.path.append(os.path.join(MG_DIR, "models"))

def loadModel(args):
  with open(os.path.join("cards", args.process, "proc_card.dat"), "r") as f:
    model_name = f.readline().strip("\n").split("import model")[1].split("-")[0].strip(" ")
  
  print(">> Loading model: %s"%model_name)  
  model = importlib.import_module(model_name)
  return model

def getNPOrders(model, args):
  cfg = tools.GetConfigFile(args.config)
  param_id_to_name = {int(each['index']): each['name'] for each in cfg['parameters']}

  orders = list(filter(lambda x: "NP" in x, dir(model.coupling_orders)))
  
  #try to find matches between parameters and NP coupling order
  param_orders = {}
  for param_id in args.pars:
    for order in orders:
      if param_id_to_name[param_id].lower() == order.split("NP")[1].lower():
        param_orders[param_id] = order
        break
  return param_orders  

def setupProcesses(args):
  """
  Need to copy the process, change the NP orders and remake it with new orders.
  A copy is needed for:
   - sm only: NP=0
   - sm-bsm interference: NP<=1 NP^2==1
   - bsm^2: NP==1
   - for each param in --pars: bsm-bsm interference: NP<=1 NP^2==2 NPc[a]^2==1
  """
  copies = od()
  copies["%s_sm_2"%args.process] = "NP=0"
  copies["%s_sm_bsm"%args.process] = "NP<=1 NP^2==1"
  copies["%s_bsm_2"%args.process] = "NP==1"

  model = loadModel(args)
  param_orders = getNPOrders(model, args)
  
  for param in args.pars:
    copies["%s_bsm_bsm_%s"%(args.process, param)] = "NP<=1 NP^2==2 %s^2==1"%param_orders[param]
  
  print(">> Need copies of the process according to:")
  for copy in copies.keys():
    print(copy, copies[copy])

  if not args.dry_run:
    for copy in copies.keys():
      if not os.path.exists("proc_dir_%s.tar.gz"%copy):
        print(">> Process directory for %s does not already exists -> will make it now"%copy)
        os.system("rm -rf cards/%s"%copy)
        os.system("rm -rf %s/%s"%(MG_DIR, copy))

        print("> Copying %s to %s"%(args.process, copy))
        os.system("cp -r cards/%s cards/%s"%(args.process, copy))
        print("> Editing proc card")
        with open(os.path.join("cards", copy, "proc_card.dat"), "r") as f:
          proc_card = f.read()
        proc_card = proc_card.replace("generate", "set group_subprocesses True \ngenerate")
        proc_card = proc_card.replace("NP<=1", copies[copy]) #replace NP coupling order
        proc_card = proc_card.replace(args.process, copy) #replace output line
        with open(os.path.join("cards", copy, "proc_card.dat"), "w") as f:
          f.write(proc_card)
        print("> New proc_card:")
        print(proc_card)
        os.system("./scripts/setup_process.sh %s"%copy)
        pwd = os.getcwd()
        os.chdir(MG_DIR)
        os.system("tar -zcf proc_dir_%s.tar.gz %s/"%(copy, copy))
        os.system("mv proc_dir_%s.tar.gz ../"%copy)
        os.system("rm -r %s"%copy)
        os.chdir(pwd)
    
def getRwNums(args):
  """
  Find out which rw points (rw_nums) we want to submit jobs for. The relevant points are those where the parameters
  provided in args.pars are turned on.
  """
  with open(os.path.join("cards", args.process, "reweight_card.dat"), "r") as f:
    rw_points = f.read().split("launch")[1:]

  cfg = tools.GetConfigFile(args.config)
  n_par = len(cfg['parameters']) 

  rw_nums = [[0, "%s_sm_2"%args.process]]

  for rw_point in rw_points:
    rw_point = rw_point.split("\n")
    rw_num = int(rw_point[0].split("=rw")[1]) #get num from "launch --rwgt_name=rw0001" 

    for set_line in rw_point[1:-1]: #set_line example: set SMEFT 6 0.1
      sett, block, param, val = set_line.split(" ")

      if float(val) != 0:
        if int(param) in args.pars:
          if float(rw_num)/2 > n_par: #cross term
            rw_nums.append([rw_num, "%s_bsm_bsm_%s"%(args.process, param)])
          elif rw_num % 2 == 0: #quadratic
            rw_nums.append([rw_num, "%s_bsm_2"%args.process])
          else: #linear
            rw_nums.append([rw_num, "%s_sm_bsm"%args.process])
          break

  return rw_nums
    
def submitJobs(args):
  try:
    os.makedirs(os.path.join("directResults", args.process))
  except OSError:
    pass

  condor_settings = CONDOR_TEMPLATE % {
    "process": args.process,
    "nevents": args.nevents,
    "ncores": args.ncores,
    "rw_nums": "\n  ".join(["%s, %s"%(each[0], each[1]) for each in args.rw_nums])
  }

  with open(os.path.join(TMPDIR, "submit.sub"), "w") as f:
    f.write(condor_settings)

  if not args.dry_run:
    os.system('condor_submit %s'%(os.path.join(TMPDIR, "submit.sub")))
  else:
    print(condor_settings)

parser = argparse.ArgumentParser()
parser.add_argument('--process', '-p', default='zh-HEL', help="Label of the process, must correspond to the dir name that was created in the MG dir")
parser.add_argument('--config', '-c', default='config.json')
parser.add_argument('--pars', type=str, nargs='+', help="List of parameters by block, e.g.: block1:1,2,3 block2:4,5,6", default=None)
parser.add_argument('--nevents', type=int, default=10000)
parser.add_argument('--ncores', type=int, default=1)
parser.add_argument('--dry-run', action="store_true", default=False)

args = parser.parse_args()

params_dict = {} #{"BLOCK1": [1,3,5], "BLOCK2":[2,3,4]}
for block_params in args.pars:
  block, params = block_params.split(":")
  params = params.split(",")
  params_dict[block] = [int(param) for param in params]

assert params_dict.keys() == ["SMEFT"]
args.pars = params_dict["SMEFT"]

args.rw_nums = getRwNums(args)
for each in args.rw_nums:
  print(each)
setupProcesses(args)
submitJobs(args)
