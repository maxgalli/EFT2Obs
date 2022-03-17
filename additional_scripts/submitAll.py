import json
import os

DRY_RUN=False

def getNEventsAndJobs(proc, info):
  total_n = info[proc][0]
  
  n_per_job = (20*60) / findTime(proc, info)
  n_jobs = total_n // n_per_job + 1
  return n_per_job, n_jobs

def findTime(proc, info):
  with open("timing_info.json", "r") as f:
    timing_info = json.load(f)
  
  if proc in timing_info.keys():
    return timing_info[proc]
  else:

    rivet = info[proc][1]
    decay = rivet == "inclusive"

    command = "python scripts/run_gridpack.py --gridpack gridpack_%s.tar.gz -s 1 -e 10 -p %s -o time_test "%(proc, rivet)

    if decay:
      command += "--rivet-ignore-beams"
    else:
      os.environ["HIGGSPRODMODE"]=info[proc][2]

    #command = "{ time " + command + " > /tmp/run_gridpack_for_timing.log ; } 2> time.txt"
    command = "{ time " + command + " > /dev/null ; } 2> time.txt"

    print(command)
    exit_status = os.system(command)
    assert exit_status == 0

    with open("time.txt", "r") as f:
      for line in f:
        if line[:4] == "real":
          minutes_and_seconds = line.split("\t")[1]

    minutes = float(minutes_and_seconds.split("m")[0])
    seconds = float(minutes_and_seconds.split("m")[1].split("s")[0])
    time = minutes*60 + seconds

    timing_info[proc] = time / 10
    with open("timing_info.json", "w") as f:
      json.dump(timing_info, f, indent=4)

    return timing_info[proc]

with open("submission_info.json", "r") as f:
  info = json.load(f)

for proc in info.keys():
  print(proc)
  if proc != "bbH_SMEFTsim_topU3l": continue
  n_per_job, n_jobs = getNEventsAndJobs(proc, info)
  
  rivet = info[proc][1]
  decay = rivet == "inclusive"
  
  command = "python scripts/launch_jobs.py --gridpack gridpack_%s.tar.gz -j %d -s 1 -e %d -p %s -o '${PWD}' --job-mode condor --dir condor/%s --task-name %s "%(proc, n_jobs, n_per_job, rivet, proc, proc)
  
  if decay:
    command += "--rivet-ignore-beams"
  else:
    command += "--env HIGGSPRODMODE=%s"%info[proc][2]

  print(command)
  if not DRY_RUN: os.system(command)
