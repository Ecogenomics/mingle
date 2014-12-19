#!/usr/bin/env python

import re

# Working with tax strings e.g. "k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pseudomonadales; f__Moraxellaceae; g__Acinetobacter; s__;":
class TaxonomyString:
    def __init__(self, string):
        self.tax_string_split = re.sub(';$','',string).split("; ")
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
