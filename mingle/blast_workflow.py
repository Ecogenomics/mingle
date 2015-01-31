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
from collections import defaultdict

from mingle.blast import BlastRunner, BlastParser
from mingle.muscle import MuscleRunner
from mingle.fasttree import FastTreeRunner
from mingle.arb_parser import ArbParser
from mingle.seq_io import SeqIO


class BlastWorkflow():
    """Blast-based workflow for building a gene tree."""

    def __init__(self, output_dir):
        """Initialization."""

        self.logger = logging.getLogger()

        self.protein_seqs = '/srv/whitlam/bio/db/gtdb/prodigal/gtdb.genes.faa'
        self.protein_blast_db = '/srv/whitlam/bio/db/gtdb/blast/prot/gtdb.genes.faa'

        self.nucleotide_seqs = '/srv/whitlam/bio/db/gtdb/prodigal/gtdb.genes.fna'
        self.nucleotide_blast_db = '/srv/whitlam/bio/db/gtdb/blast/prot/gtdb.genes.fna'

        self.blast_output = os.path.join(output_dir, "blast_out.tsv")

        self.arb_greengenes_db = os.path.join(output_dir, 'gene_tree.greengenes')

        self.homolog_output = os.path.join(output_dir, "homologs.faa")
        self.msa_output = os.path.join(output_dir, "homologs.aligned.faa")
        self.msa_log = os.path.join(output_dir, "homologs.aligned.log")

        self.tree_output = os.path.join(output_dir, "gene_tree.tree")
        self.tree_log = os.path.join(output_dir, "gene_tree.log")
        self.tree_output_log = os.path.join(output_dir, "gene_tree.out")

    def modify_greengene_hashes(self, metadata, seqs):
        """Extract and modify GreenGene hashes for homologous genes.

        Provides sensible names for homologous genes which make the
        source genome clear and handles instances where a genome
        contains multiple homologous genes. A list of dictionaries
        describing genome metadata for each homologous genes is
        produced along with a corresponding dictionary of the renamed
        homologous sequences.

        Parameters
        ----------
        metadata : dict[genome_id] -> metadata dictionary
            Metadata for genomes.
        seqs: dict[seq_id] -> seq
            Homologous sequences indexed by sequence id.

        Returns
        -------
        list of dict
           Metadata for each homologous gene.
        dict
            Homologous sequences indexed with new sequence ids.
        """

        """"""
        metadata_for_homologs = []
        genes_in_genome = defaultdict(int)
        new_seqs = {}
        for seq_id, seq in seqs.iteritems():
            # get has for genome
            genome_id = seq_id.split('_')[0]

            # create hash for gene
            new_genome_id = genome_id

            count = genes_in_genome[genome_id]
            if count != 0:
                new_genome_id += '_' + str(count + 1)

            try:
                new_hash = metadata[genome_id].copy()
            except:
                self.logger.warning('Missing metadata information for genome: %s' % genome_id)
                new_hash = {}

            new_hash['db_name'] = new_genome_id
            new_hash['prokMSA_id'] = new_genome_id
            new_hash['name'] = new_genome_id
            new_hash['acc'] = seq_id
            new_hash['ACE_genome_id'] = genome_id
            new_hash['gene_id'] = seq_id[seq_id.find('_') + 1:]
            new_hash['gene_annotation'] = ''

            metadata_for_homologs.append(new_hash)

            genes_in_genome[genome_id] += 1

            # save sequence with new name
            new_seqs[new_genome_id] = seq

        return metadata_for_homologs, new_seqs

    def run(self, query_seqs, evalue, per_identity, per_aln_len, metadata, cpus):
        """Infer a gene tree for homologous genes identified by blast.

        Complete workflow for inferring a gene tree from homologous sequences
        to a set of query sequences. Homologous sequences are identified by blast
        search and a set of user-defined parameters.

        Parameters
        ----------
        query_seqs : str
            Fasta file containing query sequences.
        evalue : float
            E-value threshold used to define homologous gene.
        per_identit : float
            Percent identity threshold used to define a homologous gene.
        per_aln_len : float
            Alignment length threshold used to define a homologous gene.
        metadata : dict[genome_id] -> metadata dictionary
            Metadata for genomes.
        cpus : int
            Number of cpus to use during homology search.
        """

        if not os.path.exists(query_seqs):
            self.logger.error('Missing query sequences file: %s' % query_seqs)
            sys.exit()

        seq_io = SeqIO()

        # identify homologous genes using blast
        self.logger.info('Identifying homologous genes using BLAST.')
        br = BlastRunner()
        br.blastp(query_seqs, self.protein_blast_db, evalue, cpus, self.blast_output)

        bp = BlastParser()
        homologs = bp.identify_homologs(self.blast_output, evalue, per_identity, per_aln_len)
        self.logger.info('  Identified %d homologous genes.' % len(homologs))

        # extract homologous sequences
        self.logger.info('Extracting homologous sequences.')
        seqs = seq_io.extract_seqs(self.protein_seqs, homologs)

        # create GreenGenes style file for ARB
        self.logger.info('Creating GreenGenes-style file for ARB.')
        metadata_for_homologs, new_seqs = self.modify_greengene_hashes(metadata, seqs)
        seq_io.write_fasta(new_seqs, self.homolog_output)

        # infer multiple sequence alignment
        self.logger.info('Inferring multiple sequence alignment.')
        mr = MuscleRunner()
        mr.run(self.homolog_output, self.msa_output, self.msa_log)
        aligned_seqs = seq_io.read_fasta(self.msa_output)

        # finish constructing out GreenGenes style ARB file
        for metadata in metadata_for_homologs:
            gene_id = metadata['db_name']
            metadata['aligned_seq'] = aligned_seqs[gene_id]

        f = open(self.arb_greengenes_db, 'w')
        arb_parser = ArbParser()
        arb_parser.write(metadata_for_homologs, f)
        f.close()

        # infer tree
        self.logger.info('Inferring gene tree.')
        ft = FastTreeRunner(multithreaded=(cpus > 1))
        ft.run(self.msa_output, 'wag', self.tree_output, self.tree_log, self.tree_output_log)
