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
__copyright__ = "Copyright 2015"
__credits__ = ["Donovan Parks"]
__license__ = "GPL3"
__maintainer__ = "Donovan Parks"
__email__ = "donovan.parks@gmail.com"
__status__ = "Development"

import os
import logging
import random
import tempfile
import shutil

from mingle.seq_io import SeqIO


class Bootstrap():
    """Calculate bootstrap support."""

    def __init__(self):
        """Initialization."""
        self.logger = logging.getLogger()

    def bootstrap(self, seqs, output_file):
        """Bootstrap multiple sequence alignment.

        Perform random sampling with replacement of
        columns within a multiple sequence alignment.

        Parameters
        ----------
        seqs : dict[seq_id] -> seq
            Sequences indexed by sequence id.
        output_file : str
            Name of file to write bootstrapped sequence alignment.
        """
        alignment_len = len(seqs[seqs.keys()[0]])
        cols = [random.randint(0, alignment_len - 1) for _ in xrange(alignment_len)]

        fout = open(output_file, 'w')
        for seqId, seq in seqs.iteritems():
            fout.write('>' + seqId + '\n')
            for col in cols:
                fout.write(seq[col])
            fout.write('\n')

        fout.close()

    def run(self, input_tree, msa_file, output_tree, num_replicates, cpus):
        """Calculate bootstraps.

        Calculate support for tree using  the non-parametric
        bootstrap methods.

        Parameters
        ----------
        input_tree : str
            Tree requiring bootstrap support values.
        msa_file : str
            Multiple sequence alignment used to infer input tree (fasta format).
        output_tree : float
            Output tree with bootstrap values.
        num_replicates : str
            Number of bootstrap replicates to perform.
        cpus : int
            Number of cpus to use.
        """

        seq_io = SeqIO()
        seqs = seq_io.read_fasta(msa_file)

        # create bootstrap trees
        self.logger.info('Creating bootstrap alignments.')
        tmp_dir = tempfile.mkdtemp()
        tree_list_file = os.path.join(tmp_dir, 'bootstrap_trees.txt')
        tree_list_out = open(tree_list_file, 'w')

        bootstrap_tree_files = []
        for i in xrange(0, num_replicates):
            bootstrap_tree_file = os.path.join(tmp_dir, 'bootstrap_tree.' + str(i) + '.tre')
            bootstrap_tree_files.append(bootstrap_tree_file)

            aln_file = os.path.join(tmp_dir, 'bootstrap_aln.' + str(i) + '.faa')
            self.bootstrap(seqs, aln_file)

            cmd = 'FastTree -quiet -nosupport -wag %s > %s 2> /dev/null\n' % (aln_file, bootstrap_tree_file)
            tree_list_out.write(cmd)

        tree_list_out.close()

        self.logger.info('Inferring bootstrap trees.')
        os.system('cat ' + tree_list_file + ' | parallel --max-procs ' + str(cpus))

        # create single file with bootstrap trees
        bootstrap_file = os.path.join(tmp_dir, 'bootstrap_trees.tre')
        bootstrap_out = open(bootstrap_file, 'w')
        for tree in bootstrap_tree_files:
            for line in open(tree):
                bootstrap_out.write(line)
        bootstrap_out.close()

        # determine bootstrap support for original tree
        self.logger.info('Determining bootstrap support for original tree.')
        os.system('CompareToBootstrap.pl -tree ' + input_tree + ' -boot ' + bootstrap_file + ' > ' + output_tree)

        # clean up temporary files
        shutil.rmtree(tmp_dir)
