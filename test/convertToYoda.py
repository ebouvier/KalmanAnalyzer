#! /usr/bin/env python

# should be called exclusively on lxplus , after sourcing 
# /afs/cern.ch/sw/lcg/external/MCGenerators_lcgcmt67c/rivet/2.2.0/x86_64-slc6-gcc47-opt/rivetenv-genser.sh

# ssh -XY ebouvier@lxplus.cern.ch
# cd /afs/cern.ch/work/e/ebouvier/CMSSW_5_3_18-TopMassSecVtx/src/
# source INIT
# cd Flat2Yoda
# scp -r bouvier@lyoserv.in2p3.fr:/gridgroup/cms/bouvier/CMSSW_5_3_18/src/UserCode/KalmanAnalyzer/test/crab_results/DATE .
# ./convertToYoda.py -d DATE
# scp DATE/*.yoda bouvier@lyoserv.in2p3.fr:/gridgroup/cms/bouvier/CMSSW_5_3_18/src/UserCode/KalmanAnalyzer/test/crab_results/DATE/

import os, shutil, subprocess, sys
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-d", "--date", dest="date", type="string", help="date of crab submission")

(options, args) = parser.parse_args()

if options.date is None or not os.path.isdir(options.date):
  parser.error("you must specify a valid date")

fullPath = options.date
os.chdir(fullPath)
for rootFile in os.listdir("./"):
    if not rootFile.endswith("merged.root"):
        continue
    if os.path.isdir(rootFile.replace(".root","")):
        os.system("rm -rf %s/*" % rootFile.replace(".root",""))
    else:
        os.mkdir(rootFile.replace(".root",""))
    os.chdir(rootFile.replace(".root",""))
    os.system("root2flat ../%s" % rootFile)
    os.system("flat2yoda MyRivetSelection/*.dat")
    os.chdir("../")
    os.system("cat %s/*.yoda > %s" % (rootFile.replace(".root",""),rootFile.replace(".root","_tmp.yoda")))
    print ("File %s has been created" % rootFile.replace(".root", "_tmp.yoda"))
for yodaFile in os.listdir("./"):
    if not yodaFile.endswith("merged_tmp.yoda"):
        continue
    os.system("sed \"s#Path=#Path=MyRivetSelection\/#\" %s > %s" % (yodaFile,yodaFile.replace("_tmp.yoda",".yoda")))
    os.remove(yodaFile)
    print ("Right paths have been set in %s" % yodaFile.replace("_tmp.yoda",".yoda"))
os.chdir("../")
