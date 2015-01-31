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

import sys
import gzip
import logging


class SeqIO():
    """Methods for reading and writing sequence files."""

    def __init__(self):
        """Initialization."""
        self.logger = logging.getLogger()

    def read_fasta(self, fasta_file):
        """Read sequences from fasta file.

        Parameters
        ----------
        fasta_file : str
            Name of fasta file to read.

        Returns
        -------
        dict
            Dictionary of sequences indexed by sequence id.
        """

        try:
            if fasta_file.endswith('.gz'):
                openFile = gzip.open
            else:
                openFile = open

            seqs = {}
            for line in openFile(fasta_file):
                # skip blank lines
                if not line.strip():
                    continue

                if line[0] == '>':
                    seqId = line[1:].split(None, 1)[0]
                    seqs[seqId] = []
                else:
                    seqs[seqId].append(line[0:-1])

            for seqId, seq in seqs.iteritems():
                seqs[seqId] = ''.join(seq)
        except:
            logger = logging.getLogger()
            logger.error("  [Error] Failed to process sequence file: " + fasta_file)
            sys.exit()

        return seqs

    def write_fasta(self, seqs, output_file):
        """Write sequences to fasta file.

        Parameters
        ----------
        seqs : dict
            Dictionary of sequences indexed by sequence id.
        output_file : str
            Name of fasta file to produce.
        """

        if output_file.endswith('.gz'):
            fout = gzip.open(output_file, 'wb')
        else:
            fout = open(output_file, 'w')

        for seq_id, seq in seqs.iteritems():
            fout.write('>' + seq_id + '\n')
            fout.write(seq + '\n')
        fout.close()

    def extract_seqs(self, fasta_file, seqs_to_extract):
        """Extract specific sequences from fasta file.

        Parameters
        ----------
        fasta_file : str
            Fasta file containing sequences.
        seqs_to_extract : set
            Ids of sequences to extract.

        Returns
        -------
        dict
            Dictionary of sequences indexed by sequence id.
        """

        seqs = {}

        for line in open(fasta_file):
            if line[0] == '>':
                seq_id = line[1:].partition(' ')[0]

                seq_of_interest = False
                if seq_id in seqs_to_extract:
                    seqs[seq_id] = []
                    seq_of_interest = True
            elif seq_of_interest:
                seqs[seq_id].append(line[0:-1])

        for seq_id, seq in seqs.iteritems():
            seqs[seq_id] = ''.join(seq)

        return seqs
