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

__author__ = "Donovan Parks"
__copyright__ = "Copyright 2014"
__credits__ = ["Donovan Parks"]
__license__ = "GPL3"
__maintainer__ = "Donovan Parks"
__email__ = "donovan.parks@gmail.com"
__status__ = "Development"

import os
import sys
import logging

from mingle.fasttree import FastTreeRunner
from mingle.seq_io import SeqIO


class Reduce():
    """Workflow for inferring a tree over a reduced set of genes."""

    def __init__(self):
        """Initialization."""
        self.logger = logging.getLogger()

    def read_ids(self, gene_id_file):
        """Read gene id file.

        Read file with gene ids to retain. Each id
        should be put on a separate line.

        Parameters
        ----------
        gene_id_file : str
            File with gene ids.

        Returns
        -------
        set
          Gene or genome ids to retain.
        """

        genes_to_retain = set()
        for line in open(gene_id_file):
            genes_to_retain.add(line.strip())

        return genes_to_retain

    def run(self, msa_file, gene_id_file, output_tree, cpus):
        """Infer a tree over a reduced set of genes.

        Filter an existing multiple sequence alignment to
        a specified set of gene ids, and infer
        tree over this reduced set of sequences.

        Parameters
        ----------
        msa_file : str
            Fasta file containing multiple sequence alignment.
        gene_ids : str
            File with gene ids to retain in tree.
        output_tree: str
            Output file containing reduced tree.
        cpus : int
            Number of cpus to use during homology search.
        """

        if not os.path.exists(msa_file):
            self.logger.error('Missing multiple sequence alignment file: %s' % msa_file)
            sys.exit()

        if not os.path.exists(gene_id_file):
            self.logger.error('Missing gene id file: %s' % gene_id_file)
            sys.exit()

        # generate msa with reduced sequences
        self.logger.info('Extracting sequences to retain.')
        genes_to_retain = self.read_ids(gene_id_file)

        seq_io = SeqIO()
        seqs = seq_io.read_fasta(msa_file)

        reduced_seqs = {}
        for seq_id, seq in seqs.iteritems():
            if seq_id in genes_to_retain:
                reduced_seqs[seq_id] = seq

        reduced_msa_file = msa_file[0:msa_file.rfind('.')]
        reduced_msa_file += '.reduced.' + msa_file[msa_file.rfind('.') + 1:]
        seq_io.write_fasta(reduced_seqs, reduced_msa_file)

        self.logger.info('  Retained %d sequences.' % len(reduced_seqs))

        # infer tree
        self.logger.info('Inferring gene tree.')
        ft = FastTreeRunner(multithreaded=(cpus > 1))
        tree_log = output_tree + '.log'
        ft.run(reduced_msa_file, 'wag', output_tree, tree_log)
