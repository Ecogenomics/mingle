#!/usr/bin/env python

###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = "Ben Woodcroft"
__copyright__ = "Copyright 2015"
__credits__ = ["Ben Woodcroft"]
__license__ = "GPL3"
__maintainer__ = "Ben Woodcroft"
__email__ = ""
__status__ = "Development"

import sys
from subprocess import Popen, PIPE, STDOUT

# Usage: alignAsNecessary.py fasta_file hmm output_stockholm

# read stdin
std = sys.stdin.read()

# if stdin empty, do nothing
# otherwise run cmd to fxtract/hmmalign
if std != '':
  command = 'fxtract -X -H -f /dev/stdin %s | hmmalign %s /dev/stdin > %s' % (
  sys.argv[1],
  sys.argv[2],
  sys.argv[3])

  pr = Popen(["/bin/bash", "-c", command], stdin=PIPE, stderr=STDOUT, stdout=PIPE)
  stdout = pr.communicate(input=std)[0]
