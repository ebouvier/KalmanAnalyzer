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
    ["\'/TTJets_FullLeptMGDecays_8TeV-madgraph/verdier-TTbar-dilept_05Oct13-v3-80a352c723b543c465f0aa98c4bbef52/USER\'", "TTJets_FullLeptMGDecays_05Oct13"],
    ["\'/TTJets_SemiLeptMGDecays_8TeV-madgraph/verdier-TTbar-semilept_05Oct13-v3-80a352c723b543c465f0aa98c4bbef52/USER\'", "TTJets_SemiLeptMGDecays_05Oct13"],
    ["\'/W1JetsToLNu_TuneZ2Star_8TeV-madgraph/verdier-W1JetsToLNu_10Oct13-v3-80a352c723b543c465f0aa98c4bbef52/USER\'", "W1JetsToLNu_10Oct13"],
    ["\'/W2JetsToLNu_TuneZ2Star_8TeV-madgraph/verdier-W2JetsToLNu_10Oct13-v3-80a352c723b543c465f0aa98c4bbef52/USER\'", "W2JetsToLNu_10Oct13"],
    ["\'/W3JetsToLNu_TuneZ2Star_8TeV-madgraph/verdier-W3JetsToLNu_10Oct13-v3-80a352c723b543c465f0aa98c4bbef52/USER\'", "W3JetsToLNu_10Oct13"],
    ["\'/W4JetsToLNu_TuneZ2Star_8TeV-madgraph/verdier-W4JetsToLNu_10Oct13-v3-80a352c723b543c465f0aa98c4bbef52/USER\'", "W4JetsToLNu_10Oct13"],
    ["\'/WW_TuneZ2star_8TeV_pythia6_tauola/verdier-WW-incl_08Oct13-v3-80a352c723b543c465f0aa98c4bbef52/USER\'", "WW-incl_08Oct13"],
    ["\'/WZ_TuneZ2star_8TeV_pythia6_tauola/verdier-WZ-incl_08Oct13-v3-80a352c723b543c465f0aa98c4bbef52/USER\'", "WZ-incl_08Oct13"],
    ["\'/ZZ_TuneZ2star_8TeV_pythia6_tauola/verdier-ZZ-incl_08Oct13-v3-80a352c723b543c465f0aa98c4bbef52/USER\'", "ZZ-incl_08Oct13"],
    ["\'/T_s-channel_TuneZ2star_8TeV-powheg-tauola/verdier-T_s-channel_09Oct13-v3-80a352c723b543c465f0aa98c4bbef52/USER\'", "T_s-channel_09Oct13"],
    ["\'/T_t-channel_TuneZ2star_8TeV-powheg-tauola/verdier-T_t-channel_09Oct13-v3-80a352c723b543c465f0aa98c4bbef52/USER\'", "T_t-channel_09Oct13"],
    ["\'/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/verdier-T_tW-channel_09Oct13-v3-80a352c723b543c465f0aa98c4bbef52/USER\'", "T_tW-channel_09Oct13"],
    ["\'/Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola/verdier-Tbar_s-channel_09Oct13-v3-80a352c723b543c465f0aa98c4bbef52/USER\'", "Tbar_s-channel_09Oct13"],
    ["\'/Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola/verdier-Tbar_t-channel_09Oct13-v3-80a352c723b543c465f0aa98c4bbef52/USER\'", "Tbar_t-channel_09Oct13"],
    ["\'/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/verdier-Tbar_tW-channel_09Oct13-v3-80a352c723b543c465f0aa98c4bbef52/USER\'", "Tbar_tW-channel_09Oct13"],
    ["\'/DY1JetsToLL_M-50_TuneZ2Star_8TeV-madgraph/verdier-DY1JetsToLL_M-50_09Oct13-v3-80a352c723b543c465f0aa98c4bbef52/USER\'", "DY1JetsToLL_M-50_09Oct13"],
    ["\'/DY2JetsToLL_M-50_TuneZ2Star_8TeV-madgraph/verdier-DY2JetsToLL_M-50_09Oct13-v3-80a352c723b543c465f0aa98c4bbef52/USER\'", "DY2JetsToLL_M-50_09Oct13"],
    ["\'/DY3JetsToLL_M-50_TuneZ2Star_8TeV-madgraph/verdier-DY3JetsToLL_M-50_09Oct13-v3-80a352c723b543c465f0aa98c4bbef52/USER\'", "DY3JetsToLL_M-50_09Oct13"],
    ["\'/DY4JetsToLL_M-50_TuneZ2Star_8TeV-madgraph/verdier-DY4JetsToLL_M-50_09Oct13-v3-80a352c723b543c465f0aa98c4bbef52/USER\'", "DY4JetsToLL_M-50_09Oct13"],

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
    pset_name = "\'mutagforrivet_cfg.py\'"
  
    print("Considering tasks for MuTag based selection")
    print("")  

    if not os.path.exists(d+"/muTag"):
        os.mkdir(d+"/muTag")

    for dataset in datasets:

        dataset_name = dataset[1]
        dataset_path = dataset[0]

        task_name = ("\'MC_muTag_%s\'") % (dataset_name)
        output_file = "%s/muTag/crab_MC_muTag_%s.py" % (d, dataset_name)
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
        pset_name = "\'d0forrivet_"+disc+"_cfg.py\'"
  
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

            task_name = ("\'MC_%s_%s\'") % (disc, dataset_name)
            output_file = "%s/%s/crab_MC_%s_%s.py" % (d, disc, disc, dataset_name)
            output_dir = ("\'crab_tasks/%s/%s\'") % (d, disc)

            print("\tCreating config file for %s" % (dataset_path))
            print("\t\tName: %s" % dataset_name)
            print("")
    
            os.system("sed -e \"s#@datasetname@#%s#\" -e \"s#@taskname@#%s#g\" -e \"s#@outputdir@#%s#g\" -e \"s#@username@#%s#g\" -e \"s#@psetname@#%s#g\" crab.cfg.template.ipnl > %s" % (dataset_path, task_name, output_dir, user_name, pset_name, output_file))
    
            cmd = "crab submit %s" % (output_file)
            if options.run:
                os.system(cmd)
