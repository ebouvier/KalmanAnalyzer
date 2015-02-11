#! /usr/bin/env python

import os, copy, datetime, pwd, re

from optparse import OptionParser
parser = OptionParser()
parser.add_option("", "--run", action="store_true", dest="run", default=False, help="run crab")
parser.add_option("-t", "--type", dest="type", type="string", default=False, help="D0 or MuTag")
(options, args) = parser.parse_args()

if not options.type or not options.type.lower() == "d0" and not options.type.lower() == "mutag":
    parser.error("you must specify a valid type")

datasets = [
    ["\'/MuHad/verdier-MuHad_Run2012A-22Jan2013_02Oct13-v3-7ed5d64fb39097b01209acde5484d3b2/USER\'", "MuHad_Run2012A-22Jan2013_02Oct13"],
    ["\'/SingleMu/verdier-SingleMu_Run2012B-22Jan2013_02Oct13-v3-7ed5d64fb39097b01209acde5484d3b2/USER\'", "SingleMu_Run2012B-22Jan2013_02Oct13"],
    ["\'/SingleMu/verdier-SingleMu_Run2012C-22Jan2013_02Oct13-v3-7ed5d64fb39097b01209acde5484d3b2/USER\'", "SingleMu_Run2012C-22Jan2013_02Oct13"],
    ["\'/SingleMu/verdier-SingleMu_Run2012D-22Jan2013_02Oct13-v3-7ed5d64fb39097b01209acde5484d3b2/USER\'", "SingleMu_Run2012D-22Jan2013_02Oct13"],
    ]

# Get username address
#user_name = "\'/store/user/%s\'" % (pwd.getpwuid(os.getuid()).pw_name)
user_name = "\'/store/user/ebouvier\'" 

d = datetime.datetime.now().strftime("%d%b%y")

version = 1

print("Creating configs for crab. Today is %s, you are %s and it's version %d" % (d, user_name, version))
print("")
if not os.path.exists(d):
    os.mkdir(d)
  
if options.type.lower() == "mutag":
    pset_name = "\'mutagforrivet_mu_cfg.py\'"
  
    print("Considering tasks for MuTag based selection")
    print("")  
    if not os.path.exists(d+"/muTag"):
        os.mkdir(d+"/muTag")

    for dataset in datasets:

        dataset_name = dataset[1]
        dataset_path = dataset[0]
 
        task_name = ("\'Data_muTag_%s\'") % (dataset_name)
        output_file = "%s/muTag/crab_Data_muTag_%s.py" % (d, dataset_name)
        output_dir = ("\'crab_tasks/%s/muTag\'") % (d)

        print("\tCreating config file for %s" % (dataset_path))
        print("\t\tName: %s" % dataset_name)
        print("")

        os.system("sed -e \"s#@datasetname@#%s#\" -e \"s#@taskname@#%s#g\" -e \"s#@outputdir@#%s#g\" -e \"s#@username@#%s#g\" -e \"s#@psetname@#%s#g\" crab.cfg.template.ipnl > %s" % (dataset_path, task_name, output_dir, user_name, pset_name, output_file))

        cmd = "crab submit %s" % (output_file)
        if options.run:
            os.system(cmd)

elif options.type.lower() == "d0":    
    for disc in ["pT", "csv"]:
        pset_name = "\'d0forrivet_"+disc+"_mu_cfg.py\'"
  
        if disc.startswith("pT"):
            print("Considering tasks for pT based selection")
        else:  
            print("Considering tasks for CSV based selection")
        print("")  
        if not os.path.exists(d+"/"+disc):
            os.mkdir(d+"/"+disc)

        for dataset in datasets:

            dataset_name = dataset[1]
            dataset_path = dataset[0]
 
            task_name = ("\'Data_%s_%s\'") % (disc, dataset_name)
            output_file = "%s/%s/crab_Data_%s_%s.py" % (d, disc, disc, dataset_name)
            output_dir = ("\'crab_tasks/%s/%s\'") % (d, disc)
    
            print("\tCreating config file for %s" % (dataset_path))
            print("\t\tName: %s" % dataset_name)
            print("")

            os.system("sed -e \"s#@datasetname@#%s#\" -e \"s#@taskname@#%s#g\" -e \"s#@outputdir@#%s#g\" -e \"s#@username@#%s#g\" -e \"s#@psetname@#%s#g\" crab.cfg.template.ipnl > %s" % (dataset_path, task_name, output_dir, user_name, pset_name, output_file))
    
            cmd = "crab submit %s" % (output_file)
            if options.run:
                os.system(cmd)

