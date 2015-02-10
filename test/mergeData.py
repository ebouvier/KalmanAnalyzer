#! /usr/bin/env python

import os, shutil, subprocess, sys
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-d", "--date", dest="date", type="string", help="date of crab submission")
parser.add_option("-t", "--type", dest="type", type="string", default=False, help="D0 or MuTag")
(options, args) = parser.parse_args()

if not options.type.lower() == "d0" and not options.type.lower() == "mutag":
    parser.error("you must specify a valid type")

if options.date is None or not os.path.isdir("crab_results/"+options.date):
    parser.error("you must specify a valid date")

os.chdir("crab_results/"+options.date+"/")  

if options.type.lower() == "d0":
    for disc in ["pT", "csv"]:
        if os.path.isfile("D0ForRivet_"+disc+"_El_merged.root") and os.path.isfile("D0ForRivet_"+disc+"_Mu_merged.root"):
            os.system("hadd D0ForRivet_"+disc+"_Data_merged.root D0ForRivet_"+disc+"_El_merged.root D0ForRivet_"+disc+"_Mu_merged.root");  
            print("D0ForRivet_"+disc+"_Data_merged.root has been produced")
        else:
            print("can't produce D0ForRivet_"+disc+"_Data_merged.root")

elif options.type.lower() == "mutag":
    for disc in ["pT", "csv"]:
        if os.path.isfile("MuTagForRivet_"+disc+"_El_merged.root") and os.path.isfile("MuTagForRivet_"+disc+"_Mu_merged.root"):
            os.system("hadd MuTagForRivet_"+disc+"_Data_merged.root D0ForRivet_"+disc+"_El_merged.root D0ForRivet_"+disc+"_Mu_merged.root");  
            print("MuTagForRivet_"+disc+"_Data_merged.root has been produced")
        else:
            print("can't produce MuTagForRivet_"+disc+"_Data_merged.root")

os.chdir("../..")

