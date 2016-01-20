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

import biolib.seq_io as seq_io
import biolib.seq_tk as seq_tk
from biolib.common import concatenate_files
from biolib.taxonomy import Taxonomy
from biolib.external.blast import Blast
from biolib.external.muscle import Muscle
from biolib.external.mafft import Mafft
from biolib.external.fasttree import FastTree
from biolib.external.execute import check_dependencies

from mingle.arb_parser import ArbParser
from mingle.common import validate_seq_ids

import dendropy


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

    def extract_homologs_and_context(self, homologs, db_file, output_file):
        """Extract homologs sequences from database file, and local gene context.

        This function extract sequences information for each
        homolog and writes this to file for downstream processing.
        In addition, it determines the local gene context for each
        gene. Specifically, it saves the annotations for the
        3 genes prior to and after a given gene.

        This function assumes the database is sorted according
        to the order genes are identified on each contig.

        Parameters
        ----------
        homologs : iterable
            Unique identifiers of sequences to extract
        db_file : str
            Fasta file with sequences.
        output_file : str
            File to write homologs.

        Returns
        -------
        dict
            d[seq_id] -> list of annotations for pre-context genes
        dict
            d[seq_id] -> list of annotations for post-context genes
        """

        gene_precontext = {}
        gene_postcontext = {}

        if len(homologs) == 0:
            return gene_precontext, gene_postcontext

        if type(homologs) is not set:
            homologs = set(homologs)

        fout = open(output_file, 'w')
        local_context = [('unknown_x~unknown', None)] * 3
        post_context_counter = {}
        for seq_id, seq, annotation in seq_io.read_fasta_seq(db_file, keep_annotation=True):
            if seq_id in homologs:
                fout.write('>' + seq_id + ' ' + annotation + '\n')
                fout.write(seq + '\n')

                gene_precontext[seq_id] = list(local_context)
                post_context_counter[seq_id] = 3

            # record 3 precontext genes
            local_context[0] = local_context[1]
            local_context[1] = local_context[2]
            local_context[2] = (seq_id, annotation)

            # record 3 postcontext genes
            if len(post_context_counter):
                key_to_remove = None
                for seq_id, count in post_context_counter.iteritems():
                    count -= 1
                    if count == -1:
                        gene_postcontext[seq_id] = list(local_context)
                        key_to_remove = seq_id
                    else:
                        post_context_counter[seq_id] = count

                if key_to_remove:
                    post_context_counter.pop(key_to_remove)

        fout.close()

        # filter gene context to contain only genes on the same scaffold
        gene_precontext = self._filter_gene_context(gene_precontext)
        gene_postcontext = self._filter_gene_context(gene_postcontext)

        return gene_precontext, gene_postcontext

    def _filter_gene_context(self, gene_context):
        """Filter gene context to contain only genes on the same scaffold.

        This function assumes sequence identifies have the following format:
            <scaffold_id>_<gene number>~<genome_id> [organism name] [IMG gene id]

        Parameters
        ----------
        gene_context : d[seq_id] -> [(seq_id, annotation), ..., (seq_id, annotation)]
            Gene context.

        Returns
        -------
        dict: d[seq_id] -> [annotation, ..., annotation]
            Filtered to contain only annotations from the same scaffold.
        """

        filtered_gene_context = {}
        for seq_id, context in gene_context.iteritems():
            gene_id = seq_id.split('~')[0]
            scaffold_id = gene_id[0:gene_id.rfind('_')]

            filtered_context = []
            for local_seq_id, annotation in context:
                local_gene_id = local_seq_id.split('~')[0]
                local_scaffold_id = local_gene_id[0:local_gene_id.rfind('_')]

                # strip organism name and IMG gene id
                annotation = annotation[0:annotation.rfind('[')]
                annotation = annotation[0:annotation.rfind('[')].strip()

                if scaffold_id == local_scaffold_id:
                    filtered_context.append(annotation)

            filtered_gene_context[seq_id] = filtered_context

        return filtered_gene_context

    def create_arb_metadata(self,
                            homologs, msa_output, taxonomy,
                            metadata,
                            gene_precontext, gene_postcontext,
                            output_file):
        """Create metadata file suitable for import into ARB.

        Parameters
        ----------
        homologs : d[seq_id] -> namedtuple of BlastHit information
            BLAST results for identified homologs.
        msa_output : str
            Fasta file with aligned homologs.
        taxonomy : d[genome_id] -> list of taxa
            Taxonomic information for genomes.
        metadata : d[key] - string
            Additional metadata to write to ARB file.
        gene_precontext : d[seq_id] -> list of annotations for pre-context genes
            Annotation for genes preceding a gene.
        gene_postcontext: d[seq_id] -> list of annotations for post-context genes
            Annotation for genes following a gene.
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
            arb_metadata['gtdb_tax_string'] = ';'.join(taxonomy.get(genome_id, ''))
            arb_metadata['aligned_seq'] = seq

            for k, v in metadata.iteritems():
                arb_metadata[k] = v

            arb_metadata['gene_precontext'] = ' -> '.join(gene_precontext.get(seq_id, []))
            arb_metadata['gene_postcontext'] = ' <- '.join(gene_postcontext.get(seq_id, []))

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
            db_file, custom_db_file,
            taxonomy_file, custom_taxonomy_file,
            evalue, per_identity, per_aln_len, max_matches, blast_mode,
            min_per_taxa, min_per_bp, restrict_taxon,
            msa_program,
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
        custom_db_file : str
            Custom database of proteins.
        taxonomy_file : str
            Taxonomic assignment of each reference genomes.
        custom_taxonomy_file : str
            Taxonomic assignment of genomes in custom database.
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
            Minimum percentage of taxa required to retain a column.
        min_per_bp : float
            Minimum percentage of base pairs required to keep trimmed sequence.
        restrict_taxon : str
            Restrict alignment to specific taxonomic group (e.g., k__Archaea).
        msa_program : str
            Program to use for multiple sequence alignment ['mafft', 'muscle'].
        output_dir : str
            Directory to store results.
        """

        assert(msa_program in ['mafft', 'muscle'])

        if not os.path.exists(query_proteins):
            self.logger.error('Missing query file: %s' % query_proteins)
            sys.exit()

        if not os.path.exists(taxonomy_file):
            self.logger.error('Missing taxonomy file: %s' % taxonomy_file)
            sys.exit()

        if not os.path.exists(db_file):
            self.logger.error('Missing database file: %s' % db_file)
            sys.exit()

        # validate query sequence names for use with mingle
        validate_seq_ids(query_proteins)

        # read taxonomy file
        self.logger.info('Reading taxonomy file.')
        taxonomy_tmp = Taxonomy().read(taxonomy_file)

        # *** [HACK] Fix genome ids to allow use of GTDBlite taxonomy files
        taxonomy = {}
        for genome_id, t in taxonomy_tmp.iteritems():
            genome_id = genome_id.replace('IMG_', '')
            taxonomy[genome_id] = t

        if custom_taxonomy_file:
            custom_taxonomy = Taxonomy().read(custom_taxonomy_file)
            taxonomy.update(custom_taxonomy)

        # identify homologs using BLASTP
        self.logger.info('Identifying homologs using %s.' % blast_mode)
        blast = Blast(self.cpus)
        blast_output = os.path.join(output_dir, 'blastp.reference_hits.tsv')
        blast.blastp(query_proteins, db_file, blast_output, evalue, max_matches, output_fmt='custom', task=blast_mode)
        homologs = blast.identify_homologs(blast_output, evalue, per_identity, per_aln_len)
        self.logger.info('Identified %d homologs in reference database.' % len(homologs))

        custom_homologs = None
        if custom_db_file:
            custom_blast_output = os.path.join(output_dir, 'blastp.custom_hits.tsv')
            blast.blastp(query_proteins, custom_db_file, custom_blast_output, evalue, max_matches, output_fmt='custom', task=blast_mode)
            custom_homologs = blast.identify_homologs(custom_blast_output, evalue, per_identity, per_aln_len)
            self.logger.info('Identified %d homologs in custom database.' % len(custom_homologs))

        # restrict homologs to specific taxonomic group
        if restrict_taxon:
            self.logger.info('Restricting homologs to %s.' % restrict_taxon)
            restricted_homologs = {}
            for query_id, hit in homologs.iteritems():
                genome_id = hit.subject_id[hit.subject_id.rfind('~') + 1:]
                if restrict_taxon in taxonomy[genome_id]:
                    restricted_homologs[query_id] = hit

            self.logger.info('%d of %d homologs are from the specified group.' % (len(restricted_homologs), len(homologs)))
            homologs = restricted_homologs

        if len(homologs) == 0:
            self.logger.error('Too few homologs were identified. Gene tree cannot be inferred.')
            sys.exit()

        # extract homologs
        self.logger.info('Extracting homologs and determining local gene context.')
        db_homologs_tmp = os.path.join(output_dir, 'homologs_db.tmp')
        gene_precontext, gene_postcontext = self.extract_homologs_and_context(homologs.keys(), db_file, db_homologs_tmp)

        homolog_ouput = os.path.join(output_dir, 'homologs.faa')
        if custom_homologs:
            custom_db_homologs_tmp = os.path.join(output_dir, 'custom_homologs_db.tmp')
            custom_gene_precontext, custom_gene_postcontext = self.extract_homologs_and_context(custom_homologs.keys(), custom_db_file, custom_db_homologs_tmp)
            gene_precontext.update(custom_gene_precontext)
            gene_postcontext.update(custom_gene_postcontext)
            homologs.update(custom_homologs)
            concatenate_files([query_proteins, db_homologs_tmp, custom_db_homologs_tmp], homolog_ouput)
            os.remove(custom_db_homologs_tmp)
        else:
            concatenate_files([query_proteins, db_homologs_tmp], homolog_ouput)

        os.remove(db_homologs_tmp)

        # infer multiple sequence alignment
        self.logger.info('Inferring multiple sequence alignment with %s.' % msa_program)

        if msa_program == 'mafft':
            mafft = Mafft(self.cpus)
            msa_output = os.path.join(output_dir, 'homologs.aligned.faa')
            msa_log = os.path.join(output_dir, 'mafft.log')
            mafft.run(homolog_ouput, msa_output, msa_log)
        elif msa_program == 'muscle':
            muscle = Muscle()
            msa_output = os.path.join(output_dir, 'homologs.aligned.faa')
            msa_log = os.path.join(output_dir, 'muscle.log')
            muscle.run(homolog_ouput, msa_output, msa_log)

        # trim multiple sequence alignment
        self.logger.info('Trimming poorly represented columns from alignment.')
        seqs = seq_io.read_fasta(msa_output, keep_annotation=True)

        for seq_id, seq in seqs.iteritems():
            if '>' in seq or '(' in seq:
                print seq_id

        trimmed_seqs, pruned_seqs = seq_tk.trim_seqs(seqs, min_per_taxa / 100.0, min_per_bp / 100.0)
        self.logger.info('%d of %d taxa were deemed to be too short and removed.' % (len(pruned_seqs), len(seqs)))

        if len(pruned_seqs) > 0:
            prune_seqs_out = os.path.join(output_dir, 'filtered_seqs.too_short.txt')
            self.logger.info('Pruned sequences written to %s.' % prune_seqs_out)
            seq_io.write_fasta(pruned_seqs, prune_seqs_out)

        if len(pruned_seqs) == len(seqs):
            self.logger.error('Too many sequences were pruned. Gene tree cannot be inferred.')
            sys.exit()

        trimmed_msa_output = os.path.join(output_dir, 'homologs.trimmed.aligned.faa')
        seq_io.write_fasta(trimmed_seqs, trimmed_msa_output)

        self.logger.info('Trimming alignment from %d bp to %d bp.' % (len(seqs.values()[0]), len(trimmed_seqs.values()[0])))

        # infer tree
        self.logger.info('Inferring gene tree.')
        fasttree = FastTree(multithreaded=(self.cpus > 1))

        tree_unrooted_output = os.path.join(output_dir, 'homologs.unrooted.tree')
        tree_log = os.path.join(output_dir, 'homologs.tree.log')
        tree_output_log = os.path.join(output_dir, 'fasttree.log')
        fasttree.run(trimmed_msa_output, 'prot', 'wag', tree_unrooted_output, tree_log, tree_output_log)

        # root tree at midpoint
        self.logger.info('Rooting tree at midpoint.')
        tree = dendropy.Tree.get_from_path(tree_unrooted_output, schema='newick', rooting="force-rooted", preserve_underscores=True)
        if len(trimmed_seqs) > 2:
            tree.reroot_at_midpoint(update_bipartitions=False)
        tree_output = os.path.join(output_dir, 'homologs.rooted.tree')
        tree.write_to_path(tree_output, schema='newick', suppress_rooting=True, unquoted_underscores=True)

        # create tax2tree consensus map and decorate tree
        self.logger.info('Decorating internal tree nodes with tax2tree.')
        output_taxonomy_file = os.path.join(output_dir, 'taxonomy.tsv')
        fout = open(output_taxonomy_file, 'w')
        for homolog_id in homologs.keys():
            genome_id = homolog_id[homolog_id.find('~') + 1:].split()[0]
            t = taxonomy.get(genome_id, None)
            if t:
                fout.write(homolog_id + '\t' + ';'.join(t) + '\n')
        fout.close()

        t2t_tree = os.path.join(output_dir, 'homologs.tax2tree.tree')
        os.system('t2t decorate -m %s -t %s -o %s' % (output_taxonomy_file, tree_output, t2t_tree))

        # setup metadata for ARB file
        src_dir = os.path.dirname(os.path.realpath(__file__))
        version_file = open(os.path.join(src_dir, 'VERSION'))

        metadata = {}
        metadata['mingle_version'] = version_file.read().strip()
        metadata['mingle_query_proteins'] = query_proteins
        metadata['mingle_db_file'] = db_file
        metadata['mingle_taxonomy_file'] = taxonomy_file
        metadata['mingle_blast_evalue'] = str(evalue)
        metadata['mingle_blast_per_identity'] = str(per_identity)
        metadata['mingle_blast_per_aln_len'] = str(per_aln_len)
        metadata['mingle_blast_max_matches'] = str(max_matches)
        metadata['mingle_blast_mode'] = blast_mode

        metadata['mingle_msa_min_per_taxa'] = str(min_per_taxa)
        metadata['mingle_msa_min_per_bp'] = str(min_per_bp)
        metadata['mingle_msa_program'] = msa_program

        # create ARB metadata file
        self.logger.info('Creating ARB metadata file.')
        arb_metadata_file = os.path.join(output_dir, 'arb.metadata.txt')
        self.create_arb_metadata(homologs, trimmed_msa_output, taxonomy,
                                 metadata,
                                 gene_precontext, gene_postcontext,
                                 arb_metadata_file)
