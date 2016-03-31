#!/usr/bin/env python3

import nose
import functools
from Bio import SeqIO
from BioExt.uds import align_to_refseq
from BioExt.scorematrices import BLOSUM62
from BioExt.references import hxb2

def setup():
    ''' Define sequence reference and records '''


@nose.with_setup(setup=setup)
def test_align_to_refseq_suffix_pad():
    ''' Ensure that sequence that ends with a '-' will not cause an error '''

    # Load reference sequence
    refseq = hxb2.prrt.load()
    seqpath = "./TEST.FASTA"

    # Load sequences
    with open(seqpath) as fh:
        seqrecords = [record for record in SeqIO.parse(fh, "fasta")]

    if len (seqrecords) == 1:
        refseq = seqrecords[0].format ('fasta')
        return {'ref': refseq, 'alignment': refseq, 'seqs': seqrecords}

    sm = BLOSUM62.load()

    all([len(seqrecord) == len(seqrecords[0]) for seqrecord in seqrecords])

    ### find the longest sequence
    msa, discarded = align_to_refseq(
        refseq,
        seqrecords,
        score_matrix=sm,
        codon=True,
        expected_identity=0.6,
        keep_insertions=False
    )

    assert msa[3].seq == seqrecords[3].seq
