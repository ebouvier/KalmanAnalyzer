#! /usr/bin/env python

import os, copy, datetime, pwd, re

from optparse import OptionParser
parser = OptionParser()
parser.add_option("", "--run", action="store_true", dest="run", default=False, help="run crab")
(options, args) = parser.parse_args()

#pset_name = "\'kalmananalyzer_el_cfg.py\'"
pset_name = "\'d0forrivet_el_cfg.py\'"

datasets = [
    ["\'/ElectronHad/verdier-ElectronHad_Run2012A-22Jan2013_02Oct13-v3-7ed5d64fb39097b01209acde5484d3b2/USER\'", "ElectronHad_Run2012A-22Jan2013_02Oct13"],
    ["\'/SingleElectron/verdier-SingleElectron_Run2012B-22Jan2013_02Oct13-v3-7ed5d64fb39097b01209acde5484d3b2/USER\'", "SingleElectron_Run2012B-22Jan2013_02Oct13"],
    ["\'/SingleElectron/verdier-SingleElectron_Run2012C-22Jan2013_02Oct13-v3-7ed5d64fb39097b01209acde5484d3b2/USER\'", "SingleElectron_Run2012C-22Jan2013_02Oct13"],
    ["\'/SingleElectron/verdier-SingleElectron_Run2012D-22Jan2013_02Oct13-v3-7ed5d64fb39097b01209acde5484d3b2/USER\'", "SingleElectron_Run2012D-22Jan2013_02Oct13"],

    ]

# Get username address
#user_name = "\'/store/user/%s\'" % (pwd.getpwuid(os.getuid()).pw_name)
user_name = "\'/store/user/ebouvier\'" 

d = datetime.datetime.now().strftime("%d%b%y")

version = 1

print("Creating configs for crab. Today is %s, you are %s and it's version %d" % (d, user_name, version))
print("")

for dataset in datasets:

  dataset_name = dataset[1]
  dataset_path = dataset[0]

  task_name = ("\'Data_%s\'") % (dataset_name)
  output_file = "crab_Data_%s.py" % (dataset_name)
  output_dir = ("\'crab_tasks/%s\'") % (d)

  print("Creating config file for %s" % (dataset_path))
  print("\tName: %s" % dataset_name)
  print("")

  os.system("sed -e \"s#@datasetname@#%s#\" -e \"s#@taskname@#%s#g\" -e \"s#@outputdir@#%s#g\" -e \"s#@username@#%s#g\" -e \"s#@psetname@#%s#g\" crab.cfg.template.ipnl > %s" % (dataset_path, task_name, output_dir, user_name, pset_name, output_file))

  cmd = "crab submit %s" % (output_file)
  if options.run:
    os.system(cmd)
