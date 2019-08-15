#!/usr/bin/env python3

import nose
import functools
import os
from Bio import SeqIO

from BioExt.io import BamIO

from BioExt.args import (
    add_alphabet,
    add_reference,
    add_scorematrix
    )

from BioExt.uds import _align_par
from BioExt.misc import compute_cigar, gapless, translate_ambiguous

from BioExt.scorematrices import (
    DNAScoreMatrix,
    FrequenciesError,
    ProteinScoreMatrix,
    BLOSUM62
    )


def setup():
    ''' Define sequence reference and records '''


@nose.with_setup(setup=setup)
def test_align():
    ''' Ensure that sequence that ends with a '-' will not cause an error '''

    dir_path = os.path.dirname(os.path.realpath(__file__))

    ## Load reference sequence
    seqpath = os.path.join(dir_path, "./rsrc/SHORT.FASTA")
    output_file = os.path.join(dir_path, "./rsrc/SHORT.FASTA.test.bam")

    records = SeqIO.parse(seqpath, 'fasta')

    reference = gapless(next(records))

    def allseqs(records):
        yield compute_cigar(reference, reference)
        for record in records:
            print(record)
            yield record

    def output(records):
        BamIO.write(
            allseqs(records),
            output_file,
            reference
            )

    _align_par(
        reference,
        records,
        BLOSUM62.load(),
        True,
        False,
        None,
        None,
        output,
        False
        )


    # Read output file
    BamIO.sort(output_file)

@nose.with_setup(setup=setup)
def test_cigar():
   pass

@nose.with_setup(setup=setup)
def test_translate_ambiguous():
    ''' Ensure ambiguous sequences are translated '''

    dir_path = os.path.dirname(os.path.realpath(__file__))

    ## Load reference sequence
    seqpath = os.path.join(dir_path, "./rsrc/SHORT.FASTA")
    output_file = os.path.join(dir_path, "./rsrc/SHORT.FASTA.test.bam")

    records = SeqIO.parse(seqpath, 'fasta')

    translated = translate_ambiguous(next(records))


