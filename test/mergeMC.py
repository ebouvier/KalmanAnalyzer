#! /usr/bin/env python

import os, shutil, subprocess, sys
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-d", "--date", dest="date", type="string", help="date of crab submission")
parser.add_option("-t", "--type", dest="type", type="string", default=False, help="D0 or MuTag")
(options, args) = parser.parse_args()

if not options.type or not options.type.lower() == "d0" and not options.type.lower() == "mutag":
    parser.error("you must specify a valid type")

if options.date is None or not os.path.isdir("crab_tasks/"+options.date):
    parser.error("you must specify a valid date")
if not os.path.exists("crab_results/"+options.date):
    os.mkdir("crab_results/"+options.date)

if options.type.lower() == "mutag":
    crabFolders = [name for name in os.listdir("crab_tasks/"+options.date+"/muTag") if os.path.isdir(os.path.join("crab_tasks/"+options.date+"/muTag", name)) and name.startswith("crab_MC_muTag_")]

    for crabFolder in crabFolders:
        datasetElt = crabFolder.rstrip("/").replace("crab_MC_muTag_", "").split("_")
        dataset = "_".join(datasetElt[:-1])
        print("Browsing %s" % dataset)

        outputName = "crab_results/"+options.date+"/MuTagForRivet_"+dataset+".root" 
        if os.path.exists(outputName):
            print "'%s' already exists." % outputName 
            continue
        command = "hadd %s " % outputName

        fullPath = "crab_tasks/%s/muTag/%s/results/" % (options.date, crabFolder)
        for name in os.listdir(fullPath):
            if name.endswith(".root"):
                command = "%s%s%s " % (command, fullPath, name)

        os.system(command)

elif options.type.lower() == "d0":
    for disc in ["pT", "csv"]:
        crabFolders = [name for name in os.listdir("crab_tasks/"+options.date+"/"+disc) if os.path.isdir(os.path.join("crab_tasks/"+options.date+"/"+disc, name)) and name.startswith("crab_MC_"+disc+"_")]

        for crabFolder in crabFolders:
            datasetElt = crabFolder.rstrip("/").replace("crab_MC_"+disc+"_", "").split("_")
            dataset = "_".join(datasetElt[:-1])
            print("Browsing %s" % dataset)

            outputName = "crab_results/"+options.date+"/D0ForRivet_"+disc+"_"+dataset+".root" 
            if os.path.exists(outputName):
                print "'%s' already exists." % outputName
                continue
            command = "hadd %s " % outputName

            fullPath = "crab_tasks/%s/%s/%s/results/" % (options.date, disc, crabFolder)
            for name in os.listdir(fullPath):
                if name.endswith(".root"):
                    command = "%s%s%s " % (command, fullPath, name)

            os.system(command)

