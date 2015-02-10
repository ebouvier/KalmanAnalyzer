#! /usr/bin/env python

import os, shutil, subprocess, sys
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-d", "--date", dest="date", type="string", help="date of crab submission")
parser.add_option("-t", "--type", dest="type", type="string", default=False, help="D0 or MuTag")
(options, args) = parser.parse_args()

if not options.type.lower() == "d0" and not options.type.lower() == "mutag":
    parser.error("you must specify a valid type")
  
if options.date is None or not os.path.isdir("crab_tasks/"+options.date):
    parser.error("you must specify a valid date")
if not os.path.exists("crab_results/"+options.date):
    os.mkdir("crab_results/"+options.date)

for disc in ["pT", "csv"]:
    if options.type.lower() == "d0":
        crabFolders = [name for name in os.listdir("crab_tasks/"+options.date+"/D0/"+disc) if os.path.isdir(os.path.join("crab_tasks/"+options.date+"/D0/"+disc, name)) and name.startswith("crab_MC_D0"+disc+"_TTJets_SemiLept")]
    elif options.type.lower() == "mutag":
        crabFolders = [name for name in os.listdir("crab_tasks/"+options.date+"/MuTag/"+disc) if os.path.isdir(os.path.join("crab_tasks/"+options.date+"/MuTag/"+disc, name)) and name.startswith("crab_MC_MuTag"+disc+"_TTJets_SemiLept")]

    if options.type.lower() == "d0":
        outputName = "crab_results/"+options.date+"/D0ForRivet_"+disc+"_TTJets_SemiLeptMGDecays.root" 
    elif options.type.lower() == "mutag":
        outputName = "crab_results/"+options.date+"/MuTagForRivet_"+disc+"_TTJets_SemiLeptMGDecays.root" 
    if os.path.exists(outputName):
        sys.exit("'%s' already exists." % outputName) 
    command = "hadd %s " % outputName
  
    for crabFolder in crabFolders:
        dataset = crabFolder.rstrip("/").replace("crab_MC_", "")
        print("Browsing %s" % dataset)
        if options.type.lower() == "d0":
            fullPath = "crab_tasks/%s/D0/%s/%s/results/" % (options.date, disc, crabFolder)
        elif options.type.lower() == "mutag":
            fullPath = "crab_tasks/%s/MuTag/%s/%s/results/" % (options.date, disc, crabFolder)
        for name in os.listdir(fullPath):
            if name.endswith(".root"):
                command = "%s%s%s " % (command, fullPath, name)

    os.system(command)

