#! /usr/bin/env python

import os, shutil, subprocess, sys
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-d", "--date", dest="date", type="string", help="date of crab submission")

(options, args) = parser.parse_args()

if options.date is None or not os.path.isdir("crab_tasks/"+options.date):
  parser.error("you must specify a valid date")

crabFolders = [name for name in os.listdir("crab_tasks/"+options.date) if os.path.isdir(os.path.join("crab_tasks/"+options.date, name)) and (name.startswith("crab_Data_Mu") or name.startswith("crab_Data_SingleMu"))]

outputName = "crab_results/"+options.date+"/"+"D0ForRivet_Mu_merged.root" 
if os.path.exists(outputName):
  sys.exit("'%s' already exists." % outputName) 
command = "hadd %s " % outputName

for crabFolder in crabFolders:
  dataset = crabFolder.rstrip("/").replace("crab_Data_", "")
  print("Browsing %s" % dataset)
  fullPath = "crab_tasks/%s/%s/results/" % (options.date, crabFolder)
  for name in os.listdir(fullPath):
      if name.endswith(".root"):
        command = "%s%s%s " % (command, fullPath, name)

if not os.path.exists("crab_results/"+options.date):
  os.mkdir("crab_results/"+options.date)

os.system(command)
