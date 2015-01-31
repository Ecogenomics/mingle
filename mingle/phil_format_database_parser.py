#!/usr/bin/env python

# Parser for this format of flat file:
#
# BEGIN
# db_name=A00000001
# organism=Korebacter versatilis Ellin345
# genome_tree_description=None
# prokMSA_id=A00000001
# owner=root
# genome_tree_tax_string=k__Bacteria; p__Acidobacteria; c__Acidobacteriia; o__Acidobacteriales; f__; g__; s__;
# greengenes_tax_string=k__Bacteria;p__Acidobacteria;c__;o__;f__Koribacteraceae;g__CandidatusKoribacter;s__CandidatusKoribacterversatilis
# blast_hits_16s=1144/1144 (100%) -- 157743 k__Bacteria; p__Acidobacteria; c__Acidobacteriia; o__Acidobacteriales; f__Koribacteraceae; g__Candidatus Koribacter; s__versatilis
# img_tax_string=k__Bacteria; p__Acidobacteria; c__unclassified; o__unclassified; f__unclassified; g__Candidatus Koribacter; s__versatilis;
# img_tax_tax2tree_corrected=
# checkm_completeness=1.0
# checkm_contamination=0.0
# core_list_status=public
# warning=
# aligned_seq=-KRTHKCGELRAADANKNVVLMGWVNRRRDLGGLIFIDLRDRTGITQIVFDNSSELQAKAGDLRSEYCIAVIGTVAKREANTVNKNLPTGEIEVVAKEMRLFNDSKVLPFSIANSNVNEEVRLKYRYLDLRRPEMQANVQMRHDVTFAIRNYLASQNFLEVETPIMTRSTPEGARDYLVPSRVHPGEFYALPQSPQIFKQ
# END
#
# BEGIN
# ...
#
# The particular keys are not necessarily as above, just the = BEGIN/END and general layout
class PhilFormatDatabaseParser:
    def __init__(self):
        self.output_order = [
                        'prokMSA_id',
                        'db_name',
                        'name',
                        'acc',
                        'organism',
                        'genome_tree_description',
                        'owner',
                        'genome_tree_tax_string',
                        'greengenes_tax_string',
                        'blast_hits_16s',
                        'img_tax_string',
                        'img_tax_tax2tree_corrected',
                        'checkm_completeness',
                        'checkm_contamination',
                        'core_list_status',
                        'warning',
                        'aligned_seq',
                        ]

    def each(self, file_handle):
        # state machine - are we inside or outside a BEGIN/END block
        state = 'outside'
        current_hash = {}

        for line_number, line in enumerate(file_handle):
            line = line.strip()
            # print "*%s*" % line
            if len(line) == 0:
                continue

            if line == 'BEGIN':
                if state == 'outside':
                    state = 'inside'
                    current_hash = {}
                else:
                    raise Exception("Badly formatted file type 1 on line %s, cannot continue, errored out at this line: %s" % (line_number + 1, line))

            elif line == 'END':
                if state == 'inside':
                    yield current_hash
                    state = 'outside'
                else:
                    raise Exception("Badly formatted file type 2 on line %s, cannot continue, errored out at this line: %s" % (line_number + 1, line))

            elif state == 'inside':
                splits = line.split('=', 1)
                if len(splits) == 2:
                    current_hash[splits[0]] = splits[1]
                else:
                    raise Exception("Badly formatted file type 3 on line %s, cannot continue, errored out at this line: %s" % (line_number + 1, line))

            else:
                raise Exception("Badly formatted file type 4 on line %s, cannot continue, errored out at this line: %s" % (line_number + 1, line))

        if state != 'outside':
            raise Exception("Badly formatted file type 5 on line %s, cannot continue, errored out at the end of the file (ended while expecting an END)" % (line_number + 1))

    def write(self, hashes, io):
        is_first = True
        for dahash in hashes:
            if is_first:
                is_first = False
            else:
                io.write('\n')

            io.write('BEGIN\n')
            for key in self.output_order:
                try:
                    io.write('='.join([key, dahash[key]]))
                    io.write('\n')
                except KeyError:
                    if key == 'warning':
                        io.write("warning=\n")  # needed for arb parsing
            io.write('END\n')
