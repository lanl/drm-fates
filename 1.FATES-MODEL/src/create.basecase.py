"""
# Creates a base case using a shell script
# Rutuja Chitra-Tarak
# July 25, 2021
"""
import os
import subprocess
runroot = os.environ["RUN_ROOT"]
archiveroot = os.environ["ARCHIVE_ROOT"]
din_loc_rootf = os.environ["DIN_LOC_ROOTF"]
subprocess.call(['sh', './src/create.basecase.sh',runroot, archiveroot, din_loc_rootf])

