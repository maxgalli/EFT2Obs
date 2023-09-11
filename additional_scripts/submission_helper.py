import json
import os
import argparse
import time

EFT2OBS_DIR = "/afs/cern.ch/work/g/gallim/MGStudies/EFT2Obs_MattSetup3"

jobFlavour_time = {
  "espresso": 20,
  "microcentury": 60,
  "longlunch": 2*60,
  "workday": 8*60,
  "tomorrow": 26*60
}

def getNEventsAndJobs(proc, info, jobFlavour):
  total_n = info[proc][0]
  
  n_per_job = int((jobFlavour_time[jobFlavour]*60) * findRate(proc, info) * (0.50))
  n_jobs = total_n // n_per_job + 1

  if n_jobs == 1:
    n_per_job = total_n

  return n_per_job, n_jobs

def findRate(proc, info):
  with open("timing_info.json", "r") as f:
    timing_info = json.load(f)
  return timing_info[proc]

def writeSubmissionInfo(directory, n_per_job, n_jobs):
  json_path = os.path.join(directory, "submission_info.json")
  if os.path.exists(json_path):
    with open(json_path, "r") as f:
      submission_info = json.load(f)
  else:
    submission_info = []

  submission_info.append([time.ctime(), n_jobs, n_per_job])
  with open(json_path, "w") as f:
    json.dump(submission_info, f, indent=4)

def submit(procs, trial_run, dry_run, local, jobFlavour, decay):
  with open("submission_info.json", "r") as f:
    info = json.load(f)

  if procs == ["all"]: procs = info.keys()

  for proc in procs:
    print("\n>> %s"%proc)
    if trial_run:
      n_per_job, n_jobs = 100, 1
      directory = "condor/trial_run/%s"%proc
    else:
      n_per_job, n_jobs = getNEventsAndJobs(proc, info, jobFlavour)
      #directory = "condor/pass1/%s"%proc
      #directory = "condor/pass2/%s"%proc
      #directory = "condor/pass3/%s"%proc
      directory = "condor/pass_acc_test/%s"%proc
    
    rivet = info[proc][1]
    
    command = "python scripts/launch_jobs.py --gridpack gridpack_%s.tar.gz -j %d -s 0 -e %d -p %s -o '${PWD}' --dir %s "%(proc, n_jobs, n_per_job, rivet, directory)
    if local: command += "--job-mode interactive --parallel 1 "
    else:     
      command += "--job-mode condor --task-name %s "%proc
      command += "--sub-opts '+JobFlavour = \"%s\"' "%jobFlavour

    if decay:
      command += "--rivet-ignore-beams"
    else:
      command += "--env HIGGSPRODMODE=%s"%info[proc][2]

    print("> %s"%command)
    if not dry_run: 
      print("> Submitting")
      os.system(command)
      os.system("ln -s %s %s"%(os.path.join(EFT2OBS_DIR, "gridpack_%s.tar.gz"%proc), os.path.join(directory, "gridpack_%s.tar.gz"%proc)))
      writeSubmissionInfo(directory, n_per_job, n_jobs)
  

def getTimeFromLog(log_file):
  with open(log_file, "r") as f:
    for line in f.readlines():
      if "Job executing on host" in line:
        start_date, start_time = line.split(" ")[2:4]
      elif "Job terminated." in line:
        end_date, end_time = line.split(" ")[2:4]

  print(start_date)

  t1 = [2022, start_date.split("/")[0], start_date.split("/")[1], start_time.split(":")[0], start_time.split(":")[1], start_time.split(":")[2], 0, 0, 0]
  t1 = [int(num) for num in t1]
  t2 = [2022, end_date.split("/")[0], end_date.split("/")[1], end_time.split(":")[0], end_time.split(":")[1], end_time.split(":")[2], 0, 0, 0]
  t2 = [int(num) for num in t2]

  return time.mktime(t2) - time.mktime(t1)

# def checkTrial():
#   os.chdir("condor/trial_run")
#   procs = os.listdir(".")

#   for proc in procs:
#     print("\n>> %s"%proc)
#     success = False
#     rivet_path = os.path.join(proc, "Rivet_0.yoda")
#     if os.path.exists(rivet_path):
#       with open(rivet_path, "r") as f:
#         content = f.read()
#       if len(content) > 0:
#         success = True

#     if success: 
#       print("> Trial \033[92m success \033[0m")
#       log_file = os.path.join(proc, list(filter(lambda f: ".log" in f, sorted(os.listdir(proc))))[-1] )
#       time = getTimeFromLog(log_file)
#       rate = 100 / time
#       print("> Rate: %.2f events/s"%rate)

#       with open("../../timing_info.json", "r") as f:
#         timing_info = json.load(f)
#       timing_info[proc] = rate
#       with open("../../timing_info.json", "w") as f:
#         json.dump(timing_info, f, indent=4)
#     else:       
#       print("> Trial \033[91m fail \033[0m")

def checkJobs(procs, trial_run, write_timing):
  if trial_run: start_dir = "condor/trial_run"
  #else:         start_dir = "condor/pass1"
  #else:         start_dir = "condor/pass2"
  #else:         start_dir = "condor/pass3"
  else:         start_dir = "condor/pass_acc_test"

  if trial_run and write_timing:
    if os.path.exists("timing_info.json"):
      with open("timing_info.json", "r") as f:
        timing_dict = json.load(f)
    else:
      timing_dict = {}

  #procs = os.listdir(start_dir)
  for proc in procs:
    successful_job = printJobStatus(os.path.join(start_dir, proc))

    if successful_job and trial_run and write_timing:
      log_file = os.path.join(start_dir, proc, list(filter(lambda f: ".log" in f, sorted(os.listdir(os.path.join(start_dir, proc)))))[-1])
      events_per_second = 100 / getTimeFromLog(log_file)
      print("> rate: %.2f events/s"%events_per_second)
      timing_dict[proc] = events_per_second

  if trial_run and write_timing:
    with open("timing_info.json", "w") as f:
      json.dump(timing_dict, f, indent=4)

def checkRivet(path):
  if os.path.exists(path):
    with open(path, "r") as f:
      content = f.read()
    if len(content) > 0:
      return True
  return False

def printJobStatus(directory):
  with open(os.path.join(directory, "submission_info.json"), "r") as f:
    submission_info = json.load(f)[-1]
  n_jobs, n_per_job = submission_info[1:]

  n_successful_jobs = 0
  #for i in range(n_jobs):
  #  rivet_path = os.path.join(directory, "Rivet_%d.yoda"%i)
  #  if checkRivet(rivet_path): n_successful_jobs += 1

  os.system("find %s -name 'Rivet_*' | wc -w > %s/checkJobs.txt"%(directory, directory))
  os.system("find %s -name 'Rivet_*' -empty | wc -w >> %s/checkJobs.txt"%(directory, directory))
  with open(os.path.join(directory, "checkJobs.txt"), "r") as f:
    n_files, n_empty = f.read().split("\n")[:2]
    n_successful_jobs = int(n_files) - int(n_empty)

  proc = directory.split("/")[-1]
  if n_successful_jobs == n_jobs: colour = '\033[92m' #green
  elif n_successful_jobs > 0:     colour = '\033[93m' #yellow
  else:                           colour = '\033[91m' #red
  print("\n>> %s %s"%(proc, colour))
  print("> %d/%d jobs completed \033[0m"%(n_successful_jobs, n_jobs))

  return n_successful_jobs > 0


def main(mode, procs, trial_run, dry_run, local, jobFlavour, write_timing, decay):
  if mode == "submit":
    submit(procs, trial_run, dry_run, local, jobFlavour, decay)
  elif mode == "checkJobs":
    checkJobs(procs, trial_run, write_timing)
  else:
    raise Exception("No such mode as %s"%mode)

if __name__=="__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('--mode', '-m', type=str, default="submit")
  parser.add_argument('--procs', '-p', type=str, nargs="+", default=["all"])
  parser.add_argument('--trial-run', '-t', action="store_true")
  parser.add_argument('--dry-run', '-d', action="store_true")
  parser.add_argument('--local', '-l', action="store_true")
  parser.add_argument('--jobFlavour', '-f', default="espresso")
  parser.add_argument('--write-timing', action="store_true")
  parser.add_argument('--decay', action="store_true", default=False)


  args = parser.parse_args()

  main(args.mode, args.procs, args.trial_run, args.dry_run, args.local, args.jobFlavour, args.write_timing, args.decay)
  
