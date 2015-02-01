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

import re


# Working with tax strings e.g. "k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pseudomonadales; f__Moraxellaceae; g__Acinetobacter; s__;":
class TaxonomyString:
    def __init__(self, string):
        self.tax_string_split = re.sub(';$', '', string).split("; ")
        self.original_name = string

    def iterate(self):
        reg = re.compile('^.__;{0,1}$')
        for rank_i, taxon in enumerate(self.tax_string_split):
            # get out if there is no further info
            if reg.match(taxon):
                break
            else:
                yield taxon

    def full_name(self):
        current = ''
        for rank_i, taxon in enumerate(self.iterate()):
            if rank_i == 0:
                current = taxon
            else:
                current = "%s; %s" % (current, taxon)
        return current

    def names(self):
        names = []
        for i in self.iterate():
            names.append(i)
        return names

    def num_levels(self):
        return len(self.tax_string_split)
