#!/usr/bin/env python

import argparse
import logging
import sys
import os
import subprocess
import re
from Bio import SeqIO
from subprocess import Popen, PIPE, STDOUT
import StringIO

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','lib'))
from phil_format_database_parser import PhilFormatDatabaseParser
from taxonomy_string import TaxonomyString


# Inputs:
# greengenes file from genome tree database => provides taxonomy
# HMM => the HMM to create a reference package out of
# path to proteome files - one for each genome => sequences to be put into the reference package
parser = argparse.ArgumentParser(description='''Creates a fasta file of each HMM's hit, where each is hit is annotated by its taxonomy''')
parser.add_argument('--aligned_fasta', help = 'aligned file output from createReference.py', required = True)
parser.add_argument('--greengenes', help = '.greengenes file to define taxonomy', required = True)
parser.add_argument('--proteomes_list_file', help = 'path to file containing faa files of each genome\'s proteome (of the form [path]<ace_ID>[stuff].faa e.g. /srv/home/ben/A00000001.fna.faa', required = True)
parser.add_argument('--hit_fasta_file', help = 'output taxonomy-annotated, unaligned FASTA file of hits to this file', required = True)

options = parser.parse_args()

# module load taxtastic hmmer fxtract fasttree
logging.basicConfig(level=logging.DEBUG)

# Read taxonomy into hash
logging.info('Reading taxonomy information from .greengenes file..')
ace_id_to_taxonomy = {}
for entry in PhilFormatDatabaseParser().each(open(options.greengenes)):
  ace_id = entry['db_name']
  try:
    taxonomy = entry['genome_tree_tax_string']
    ace_id_to_taxonomy[ace_id] = TaxonomyString(taxonomy).full_name()
  except KeyError:
    logging.warn("taxonomy information not found in Phil format file for ID: %s, skipping" % ace_id)


# Read file of reference faa file paths into hash of ACE_ID => faa file
ace_id_to_proteome_file = {}
with open(options.proteomes_list_file) as f:
  for pro in f:
    pro = pro.strip()

    ace_id = re.match(r'^([A-Z]\d+)',os.path.basename(pro)).group(1)
    ace_id_to_proteome_file[ace_id] = pro
if len(ace_id_to_proteome_file) < 1: raise Exception("Error: no proteome files found, cannot continue")
logging.info("Read in %s proteome files for processing e.g. %s => %s" % (len(ace_id_to_proteome_file), ace_id_to_proteome_file.keys()[0], ace_id_to_proteome_file[ace_id_to_proteome_file.keys()[0]]))



def extract_and_output(ace_id, ace_id_to_taxonomy, ace_id_to_proteome_file, current_protein_ids, output_fh):
  try:
    taxonomy = ace_id_to_taxonomy[ace_id]
  except KeyError:
    logging.warn("Taxonomy not found for ACE ID %s, skipping" % ace_id)
    taxonomy = 'unknown'

  proteome_file = ace_id_to_proteome_file[ace_id]
  std = "\n".join(current_protein_ids)
  command = 'fxtract -X -H -f /dev/stdin %s' % proteome_file

  pr = Popen(["/bin/bash", "-c", command], stdin=PIPE, stderr=STDOUT, stdout=PIPE)
  stdout = pr.communicate(input=std)[0]

  num_written = 0
  #import code; code.interact(local=locals())
  for record in SeqIO.parse(StringIO.StringIO(stdout), 'fasta'):
    output_fh.write('>%s_%s %s\n' % (ace_id, record.id, taxonomy))
    output_fh.write('%s\n' % record.seq)
    num_written += 1

  if num_written != len(current_protein_ids):
    import code; code.interact(local=locals())
    raise Exception("Unexpected number of proteins retrieved from fasta file!")






# Read aligned fasta file one line at a time,
with open(options.hit_fasta_file,'w') as output_fh:
  with open(options.aligned_fasta) as f:
    current_ace_id = None
    current_protein_ids = []
    for record in SeqIO.parse(f, 'fasta'):
      # parse ID into ACE_ID and protein ID
      # e.g. A00000001_1_NC_008009_1581
      # means ACE ID "A00000001" and protein ID "NC_008009_1581"
      matches = re.match(r'^([A-Z01-9]+)_\d+_(.+)$', record.id)
      if matches is None:
        raise Exception("Unparsable ID")
      #      import code; code.interact(local=locals())
      ace_id = matches.group(1)
      protein_id = matches.group(2)
      if current_ace_id is not None and ace_id != current_ace_id:
        # if the ACE ID has changed, then fxtract the proteins from the last ACE ID's faa file, adding the taxonomy to the ID line and print the now annotated fasta file
        extract_and_output(current_ace_id, ace_id_to_taxonomy, ace_id_to_proteome_file, current_protein_ids, output_fh)
        current_ace_id = ace_id
        current_protein_ids = [protein_id]
      else:
        current_ace_id = ace_id
        current_protein_ids.append(protein_id)

    if current_ace_id is None:
      raise Exception("No sequences parsed from aligned fasta file, it seems")
    else:
      # Parse the last one
      extract_and_output(ace_id, ace_id_to_taxonomy, ace_id_to_proteome_file, current_protein_ids, output_fh)





