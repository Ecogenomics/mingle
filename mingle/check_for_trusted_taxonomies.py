#!/usr/bin/env python

# Given an unpdated file of taxonomies, check the tax string has
# taxonomy consistent with trusted taxonomy therein.

import fileinput

class Check_trusted_tax:
    # Define the list of trusted taxonomies.
    def load_trusted_tax(self):
        tax = {}
        
        for taxonomy in open('/srv/home/s4293811/gits/git_mingle/mingle/mingle/trusted_taxonomies.txt'):
            split = taxonomy.split(',')
            old_tax = split[0]
            new_tax = split[1]
            tax[old_tax] = new_tax
        
        return tax
    
    def check_tax(self, tax_file, trusted_tax):
                
        for old_tax, new_tax in trusted_tax.iteritems():
            
            for line in tax_file:
                
                if tax.startswith(old_tax):
                    tax.replace(old_tax, new_tax)