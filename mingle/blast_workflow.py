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

import biolib.seq_io as seq_io
import biolib.seq_tk as seq_tk
from biolib.common import concatenate_files
from biolib.taxonomy import Taxonomy
from biolib.external.blast import Blast
from biolib.external.muscle import Muscle
from biolib.external.fasttree import FastTree
from biolib.external.execute import check_dependencies

from mingle.arb_parser import ArbParser


class BlastWorkflow():
    """Blast-based workflow for building a gene tree."""

    def __init__(self, cpus):
        """Initialization.

        Parameters
        ----------
        cpus : int
            Number of cpus to use during homology search.
        """

        check_dependencies(['muscle', 'FastTreeMP', 'blastp', 't2t'])

        self.logger = logging.getLogger()

        self.cpus = cpus

    def extract_homologs(self, homologs, db_file, output_file):
        """Extract homologs sequences from database file.

        Parameters
        ----------
        homologs : iterable
            Unique identifiers of sequences to extract
        db_file : str
            Fasta file with sequences.
        output_file : str
            File to write homologs.
        """

        if type(homologs) is not set:
            homologs = set(homologs)

        fout = open(output_file, 'w')
        for seq_id, seq, annotation in seq_io.read_fasta_seq(db_file, keep_annotation=True):
            if seq_id in homologs:
                fout.write('>' + seq_id + ' ' + annotation + '\n')
                fout.write(seq + '\n')
        fout.close()

    def create_arb_metadata(self, homologs, msa_output, taxonomy, output_file):
        """Create metadata file suitable for import into ARB.

        Parameters
        ----------
        homologs : d[seq_id] -> namedtuple of BlastHit information
            BLAST results for identified homologs.
        msa_output : str
            Fasta file with aligned homologs.
        taxonomy : d[genome_id] -> list of taxa
            Taxonomic information for genomes.
        output_file : str
            File to write metadata information.
        """

        arb_metadata_list = []
        for seq_id, seq, annotation in seq_io.read_seq(msa_output, keep_annotation=True):
            if '~' in seq_id:
                scaffold_gene_id = seq_id[:seq_id.find('~')]
                genome_id = seq_id[seq_id.find('~') + 1:].split()[0]
            else:
                scaffold_gene_id = seq_id
                genome_id = ''

            arb_metadata = {}
            arb_metadata['db_name'] = seq_id
            arb_metadata['img_genome_id'] = genome_id
            arb_metadata['img_scaffold_id'] = scaffold_gene_id[0:scaffold_gene_id.rfind('_')]
            arb_metadata['img_scaffold_gene_id'] = scaffold_gene_id
            arb_metadata['img_tax_string'] = ';'.join(taxonomy.get(genome_id, ''))
            arb_metadata['aligned_seq'] = seq

            hit_info = homologs.get(seq_id, None)
            if hit_info:
                arb_metadata['blast_evalue'] = '%.1g' % hit_info.evalue
                arb_metadata['blast_bitscore'] = '%.1f' % hit_info.bitscore
                arb_metadata['blast_perc_identity'] = '%.1f' % hit_info.perc_identity
                arb_metadata['blast_subject_perc_alignment_len'] = '%.1f' % hit_info.subject_perc_aln_len
                arb_metadata['blast_query_perc_alignment_len'] = '%.1f' % hit_info.query_perc_aln_len
                arb_metadata['blast_query_id'] = hit_info.query_id

            if annotation:
                annotation_split = annotation.split('[')
                if len(annotation_split) == 3:
                    # assume format is <annotation> [<genome name>] [<IMG gene id>]
                    gene_annotation, organism_name, gene_id = annotation_split
                    organism_name = organism_name.replace(']', '')
                    gene_id = gene_id.replace(']', '').replace('IMG Gene ID: ', '')
                elif len(annotation_split) == 2:
                    # format is essentially unknown, but the most likely issue
                    # is that the gene itself just doesn't have an annotation
                    gene_annotation = ''
                    organism_name, gene_id = annotation_split
                    organism_name = organism_name.replace(']', '')
                    gene_id = gene_id.replace(']', '').replace('IMG Gene ID: ', '')
                else:
                    # no idea what the format is, so just save the annotation
                    gene_annotation = annotation
                    organism_name = ''
                    gene_id = ''

                arb_metadata['img_gene_annotation'] = gene_annotation
                arb_metadata['organism'] = organism_name
                arb_metadata['full_name'] = organism_name
                arb_metadata['img_gene_id'] = gene_id

            arb_metadata_list.append(arb_metadata)

        fout = open(output_file, 'w')
        arb_parser = ArbParser()
        arb_parser.write(arb_metadata_list, fout)
        fout.close()

    def run(self, query_proteins,
            db_file, taxonomy_file,
            evalue, per_identity, per_aln_len, max_matches, blast_mode,
            min_per_taxa, min_per_bp,
            output_dir):
        """Infer a gene tree for homologs genes identified by blast.

        Workflow for inferring a gene tree from sequences identified as being
        homologs to a set of query proteins. Homologs are identified using BLASTP
        and a set of user-defined parameters.

        Parameters
        ----------
        query_proteins : str
            Fasta file containing query proteins.
        db_file : str
            BLAST database of reference proteins.
        taxonomy_file : str
            Taxonomic assignment of each reference genomes.
        evalue : float
            E-value threshold used to define homolog.
        per_identity : float
            Percent identity threshold used to define a homolog.
        per_aln_len : float
            Alignment length threshold used to define a homolog.
        max_matches : int
            Maximum matches per query protein.
        metadata : dict[genome_id] -> metadata dictionary
            Metadata for genomes.
        blast_mode : str
            Type of blast to perform.
        min_per_taxa : float
            Minimum percentage of taxa required to retain a leading and trailing columns.
        min_per_bp : float
            Minimum percentage of base pairs required to keep trimmed sequence.
        output_dir : str
            Directory to store results.
        """

        if not os.path.exists(query_proteins):
            self.logger.error('Missing query file: %s' % query_proteins)
            sys.exit()

        if not os.path.exists(taxonomy_file):
            self.logger.error('Missing taxonomy file: %s' % taxonomy_file)
            sys.exit()

        if not os.path.exists(db_file):
            self.logger.error('Missing database file: %s' % db_file)
            sys.exit()

        # read taxonomy file
        self.logger.info('Reading taxonomy file.')
        taxonomy = Taxonomy().read(taxonomy_file)

        # identify homologs using BLASTP
        self.logger.info('Identifying homologs using %s.' % blast_mode)
        blast = Blast(self.cpus)
        blast_output = os.path.join(output_dir, 'blastp.hits.tsv')
        blast.blastp(query_proteins, db_file, blast_output, evalue, max_matches, output_fmt='custom', task=blast_mode)
        homologs = blast.identify_homologs(blast_output, evalue, per_identity, per_aln_len, output_fmt='custom')
        self.logger.info('Identified %d homologs.' % len(homologs))

        # extract homologs
        self.logger.info('Extracting homologs.')
        db_homologs_tmp = os.path.join(output_dir, 'homologs_db.tmp')
        self.extract_homologs(homologs.keys(), db_file, db_homologs_tmp)
        homolog_ouput = os.path.join(output_dir, 'homologs.faa')
        concatenate_files([query_proteins, db_homologs_tmp], homolog_ouput)
        os.remove(db_homologs_tmp)

        # infer multiple sequence alignment
        self.logger.info('Inferring multiple sequence alignment.')
        muscle = Muscle()
        msa_output = os.path.join(output_dir, 'homologs.aligned.faa')
        msa_log = os.path.join(output_dir, 'muscle.log')
        muscle.run(homolog_ouput, msa_output, msa_log)

        # trim multiple sequence alignment
        self.logger.info('Trimming leading and trailing columns of alignment.')
        seqs = seq_io.read_fasta(msa_output, keep_annotation=True)
        trimmed_seqs, pruned_seqs = seq_tk.trim_seqs(seqs, min_per_taxa / 100.0, min_per_bp / 100.0)
        trimmed_msa_output = os.path.join(output_dir, 'homologs.trimmed.aligned.faa')
        seq_io.write_fasta(trimmed_seqs, trimmed_msa_output)
        self.logger.info('Trimming sequences from %d bp to %d bp.' % (len(seqs.values()[0]), len(trimmed_seqs.values()[0])))
        self.logger.info('%d of %d taxa were deemed to be too short and removed.' % (len(pruned_seqs), len(seqs)))
        os.remove(msa_output)

        # infer tree
        self.logger.info('Inferring gene tree.')
        fasttree = FastTree(multithreaded=(self.cpus > 1))

        tree_output = os.path.join(output_dir, 'homologs.tree')
        tree_log = os.path.join(output_dir, 'homologs.tree.log')
        tree_output_log = os.path.join(output_dir, 'fasttree.log')
        fasttree.run(trimmed_msa_output, 'prot', 'wag', tree_output, tree_log, tree_output_log)

        # create tax2tree consensus map and decorate tree
        self.logger.info('Decorating internal tree nodes with tax2tree.')
        output_taxonomy_file = os.path.join(output_dir, 'taxonomy.tsv')
        fout = open(output_taxonomy_file, 'w')
        for homolog_id in homologs.keys():
            genome_id = homolog_id[homolog_id.find('~') + 1:].split()[0]
            fout.write(homolog_id + '\t' + ';'.join(taxonomy[genome_id]) + '\n')
        fout.close()

        t2t_tree = os.path.join(output_dir, 'homologs.tax2tree.tree')
        os.system('t2t decorate -m %s -t %s -o %s' % (output_taxonomy_file, tree_output, t2t_tree))

        # create ARB metadata file
        self.logger.info('Creating ARB metadata file.')
        arb_metadata_file = os.path.join(output_dir, 'arb.metadata.txt')
        self.create_arb_metadata(homologs, trimmed_msa_output, taxonomy, arb_metadata_file)
