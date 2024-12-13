#!/usr/bin/env python3

import re
import subprocess
import sys

jobid = sys.argv[1]

try:
    output = str(subprocess.check_output("sacct -j %s --format State --noheader" % jobid, shell=True))
except subprocess.CalledProcessError as e:
    print("FAILED")
    exit(0)

state = output.strip().split()[0]

running_status=["PENDING", "CONFIGURING", "COMPLETING", "RUNNING", "SUSPENDED"]
if "COMPLETED" in state:
    print("success")
elif any(s in state for s in running_status):
    print("running")
else:
    print("failed") 