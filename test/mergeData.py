#! /usr/bin/env python

import os, shutil, subprocess, sys
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-d", "--date", dest="date", type="string", help="date of crab submission")

(options, args) = parser.parse_args()

if options.date is None or not os.path.isdir("crab_results/"+options.date):
  parser.error("you must specify a valid date")

os.chdir("crab_results/"+options.date)  

for disc in ["pT", "csv"]:
  if os.path.isfile("D0ForRivet_"+disc+"_El_merged.root") and os.path.isfile("D0ForRivet_"+disc+"_Mu_merged.root"):
    os.system("hadd D0ForRivet_"+disc+"_Data_merged.root D0ForRivet_"+disc+"_El_merged.root D0ForRivet_"+disc+"_Mu_merged.root");  
    print("D0ForRivet_"+disc+"_Data_merged.root has been produced")
  else:
    print("can't produce D0ForRivet_"+disc+"_Data_merged.root")

os.chdir("../..")

