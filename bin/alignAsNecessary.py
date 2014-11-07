#!/usr/bin/env python

import sys
from subprocess import Popen, PIPE, STDOUT

# Usage: alignAsNecessary.py fasta_file hmm output_stockholm

#read stdin
std = sys.stdin.read()

#if stdin empty, do nothing
# otherwise run cmd to fxtract/hmmalign
if std != '':
  command = 'fxtract -X -H -f /dev/stdin %s | hmmalign %s /dev/stdin > %s' % (
  sys.argv[1],
  sys.argv[2],
  sys.argv[3])

  pr = Popen(["/bin/bash", "-c", command], stdin=PIPE, stderr=STDOUT, stdout=PIPE)
  stdout = pr.communicate(input=std)[0]
