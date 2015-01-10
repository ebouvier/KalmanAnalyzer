#! /usr/bin/env python

import os, shutil, subprocess, sys
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-d", "--date", dest="date", type="string", help="date of crab submission")

(options, args) = parser.parse_args()

if options.date is None or not os.path.isdir("crab_results/"+options.date):
  parser.error("you must specify a valid date")

fullPath = "crab_results/"+options.date
os.chdir(fullPath)
for rootFile in os.listdir("./"):
    if not rootFile.endswith("merged.root"):
        continue
    os.system("root2flat %s" % rootFile)
    fileName = rootFile.replace(".root", ".dat")
    os.system("cat MyRivetSelection/*.dat > %s" % fileName)
    os.system("flat2aida %s" % fileName)
    os.system("rm -rf MyRivetSelection %s" % fileName)
    print ("File %s has been created" % fileName.replace(".dat", ".aida"))
os.chdir("../..")
