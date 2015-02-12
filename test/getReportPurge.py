#! /usr/bin/env python

import os, copy, datetime, pwd, re

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-d", "--date", dest="date", type="string", default=False, help="date when the task was submitted")
parser.add_option("-t", "--type", dest="type", type="string", default=False, help="csv, pT, or muTag")
parser.add_option("", "--status", action="store_true", dest="status", default=False, help="status of all tasks")
parser.add_option("", "--get", action="store_true", dest="get", default=False, help="getoutput of all tasks")
parser.add_option("", "--report", action="store_true", dest="report", default=False, help="report all tasks")
parser.add_option("", "--purge", action="store_true", dest="purge", default=False, help="purge all tasks")
parser.add_option("", "--kill", action="store_true", dest="kill", default=False, help="kill all tasks")
(options, args) = parser.parse_args()

if not options.date or not os.path.isdir(os.path.join("crab_tasks", options.date)):
    parser.error("you must specify a valid date")
rootName = os.path.join("crab_tasks", options.date)
if not options.type or not os.path.isdir(os.path.join(rootName, options.type)):
    parser.error("you must specify a valid type")
rootName = os.path.join(rootName, options.type)
if not options.status and not options.get and not options.report and not options.purge and not options.kill:
    parser.error("you must specify a valid action")
    
crabFolders = [name for name in os.listdir(rootName) if os.path.isdir(os.path.join(rootName, name))]

for crabFolder in crabFolders:

    folderName = os.path.join(rootName, crabFolder)

    if options.status:
        cmd = "crab status " + folderName
    if options.get:
        cmd = "crab getoutput " + folderName
    if options.report:
        cmd = "crab report " + folderName
    if options.purge:
        cmd = "crab purge --cache " + folderName
    if options.kill:
        cmd = "crab kill " + folderName
        
    os.system(cmd)

