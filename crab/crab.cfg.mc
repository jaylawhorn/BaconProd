[CRAB]

jobtype = cmssw
scheduler = remoteGlidein
use_server = 0

[CMSSW]

datasetpath = XX-DATASET-XX

pset = makingBacon_MC.py

total_number_of_events = -1
events_per_job = XX-NUM-XX

output_file = Output.root

[USER]

return_data = 0
copy_data = 1
storage_element = XX-STORAGE-ELEMENT-XX
storage_path = XX-STORAGE-PATH-XX
user_remote_dir = XX-USER-REMOTE-DIR-XX

[GRID]

rb = CERN
se_black_list = T0,T1
#se_white_list =
#ce_black_list =
#ce_white_list =
