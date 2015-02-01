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

import Bio
import re

# Slight modification of the biopython stockholm format parser to
# also parse in global GC annoations rather than ignoring them.
# Also the interface is different - this class should be called
# directly, not through AlignIO.
#
# Copyright 2006-2013 by Peter Cock.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
class GraftMStockholmIterator:
    def readAlignedSequences(self, io):
        line = io.readline()
        if not line.strip() == '# STOCKHOLM 1.0':
            raise ValueError("Did not find STOCKHOLM header")

        # Note: If this file follows the PFAM conventions, there should be
        # a line containing the number of sequences, e.g. "#=GF SQ 67"
        # We do not check for this - per we should, and verify that
        # if present it agrees with our parsing.
        seqs = {}
        ids = []
        gs = {}
        gr = {}
        gf = {}
        gc = {}
        passed_end_alignment = False
        while True:
            line = io.readline()
            if not line:
                break  # end of file
            line = line.strip()  # remove trailing \n
            if line == '# STOCKHOLM 1.0':
                self._header = line
                break
            elif line == "//":
                # The "//" line indicates the end of the alignment.
                # There may still be more meta-data
                passed_end_alignment = True
            elif line == "":
                # blank line, ignore
                pass
            elif line[0] != "#":
                # Sequence
                # Format: "<seqname> <sequence>"
                assert not passed_end_alignment
                parts = [x.strip() for x in line.split(" ", 1)]
                if len(parts) != 2:
                    # This might be someone attempting to store a zero length sequence?
                    raise ValueError("Could not split line into identifier "
                                      + "and sequence:\n" + line)
                id, seq = parts
                if id not in ids:
                    ids.append(id)
                seqs.setdefault(id, '')
                seqs[id] += seq
            elif len(line) >= 5:
                # Comment line or meta-data
                if line[:5] == "#=GF ":
                    # Generic per-File annotation, free text
                    # Format: #=GF <feature> <free text>
                    feature, text = line[5:].strip().split(None, 1)
                    # Each feature key could be used more than once,
                    # so store the entries as a list of strings.
                    if feature not in gf:
                        gf[feature] = [text]
                    else:
                        gf[feature].append(text)
                elif line[:5] == '#=GC ':
                    # Generic per-Column annotation, exactly 1 char per column
                    # Format: "#=GC <feature> <exactly 1 char per column>"
                    pass
                elif line[:5] == '#=GS ':
                    # Generic per-Sequence annotation, free text
                    # Format: "#=GS <seqname> <feature> <free text>"
                    id, feature, text = line[5:].strip().split(None, 2)
                    # if id not in ids:
                    #    ids.append(id)
                    if id not in gs:
                        gs[id] = {}
                    if feature not in gs[id]:
                        gs[id][feature] = [text]
                    else:
                        gs[id][feature].append(text)
                elif line[:5] == "#=GR ":
                    # Generic per-Sequence AND per-Column markup
                    # Format: "#=GR <seqname> <feature> <exactly 1 char per column>"
                    id, feature, text = line[5:].strip().split(None, 2)
                    # if id not in ids:
                    #    ids.append(id)
                    if id not in gr:
                        gr[id] = {}
                    if feature not in gr[id]:
                        gr[id][feature] = ""
                    gr[id][feature] += text.strip()  # append to any previous entry
                    # TODO - Should we check the length matches the alignment length?
                    #       For iterlaced sequences the GR data can be split over
                    #       multiple lines
            # Next line...

        assert len(seqs) <= len(ids)
        # assert len(gs)   <= len(ids)
        # assert len(gr)   <= len(ids)

        self.ids = ids
        self.sequences = seqs
        self.seq_annotation = gs
        self.seq_col_annotation = gr
        self.gc = gc
        for k, v in seqs.items():
            seqs[k] = re.sub(r'[a-z\.\*]', '', v, 0)

        return seqs
