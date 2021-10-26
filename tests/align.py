#!/usr/bin/env python3

import nose

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from BioExt.align import Aligner
from BioExt.scorematrices import BLOSUM62

aln = Aligner(BLOSUM62.load())

class test_Aligner():
    def test_align_self(self):
        """Check alignment of reference against itself"""
        assert aln('GCTAGA', 'GCTAGA') == (4.5, 'GCTAGA', 'GCTAGA')

    def test_align_self_case(self):
        """Check case is irrelevant for self-alignment"""
        assert aln('GCTAGA', 'GCTAGA') == aln('GCTAGA', 'GcTaGa')

    def test_align_self_seq_to_str(self):
        """Check alignment of Seq instance against seq-identical str"""
        ref = 'GCTAGA'
        record = Seq(ref)
        assert aln(ref, record) == (4.5, ref, record)

    def test_align_self_seqrecord_to_str(self):
        """Check alignment of SeqRecord instance against seq-identical str"""
        ref = 'GCTAGA'
        record = SeqRecord(Seq(ref))
        score, ref_aln, query_aln = aln(ref, record)
        assert (score, ref_aln, str(query_aln.seq)) == (4.5, ref, ref)
