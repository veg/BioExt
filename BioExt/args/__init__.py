
from __future__ import division, print_function


__all__ = [
    'add_alphabet',
    'add_reference',
    'add_scorematrix'
    ]


def add_alphabet(parser, *args):
    kwargs = dict(
        metavar='ALPHABET',
        choices=('amino', 'dna', 'codon'),
        default='codon',
        help='perform an alignment using one of {{{0}}} [default=codon]'.format(
            ', '.join(('amino', 'dna', 'codon'))
            )
        )

    parser.add_argument(*args, **kwargs)

    return parser


def add_reference(parser, *args):
    from argparse import ArgumentTypeError
    from Bio import SeqIO
    from BioExt.references import hxb2, nl4_3, cov2

    references = {
        'HXB2_env': hxb2.env,
        'HXB2_gag': hxb2.gag,
        'HXB2_int': hxb2.int,
        'HXB2_nef': hxb2.nef,
        'HXB2_pol': hxb2.pol,
        'HXB2_pr': hxb2.pr,
        'HXB2_prrt': hxb2.prrt,
        'HXB2_rev': hxb2.rev,
        'HXB2_rt': hxb2.rt,
        'HXB2_tat': hxb2.tat,
        'HXB2_vif': hxb2.vif,
        'HXB2_vpr': hxb2.vpr,
        'HXB2_vpu': hxb2.vpu,
        'NL4-3_prrt': nl4_3.prrt,
        'CoV2-3C': cov2.threeC,
        'CoV2-E': cov2.E,
        'CoV2-endornase': cov2.endornase,
        'CoV2-exonuclease': cov2.exonuclease,
        'CoV2-helicase': cov2.helicase,
        'CoV2-leader': cov2.leader,
        'CoV2-methyltransferase': cov2.methyltransferase,
        'CoV2-M': cov2.M,
        'CoV2-N': cov2.N,
        'CoV2-nsp10': cov2.nsp10,
        'CoV2-nsp2': cov2.nsp2,
        'CoV2-nsp3': cov2.nsp3,
        'CoV2-nsp4': cov2.nsp4,
        'CoV2-nsp6': cov2.nsp6,
        'CoV2-nsp7': cov2.nsp7,
        'CoV2-nsp8': cov2.nsp8,
        'CoV2-nsp9': cov2.nsp9,
        'CoV2-ORF10': cov2.ORF10,
        'CoV2-ORF1a': cov2.ORF1a,
        'CoV2-ORF1b': cov2.ORF1b,
        'CoV2-ORF3a': cov2.ORF3a,
        'CoV2-ORF5': cov2.ORF5,
        'CoV2-ORF6': cov2.ORF6,
        'CoV2-ORF7a': cov2.ORF7a,
        'CoV2-ORF7b': cov2.ORF7b,
        'CoV2-ORF8': cov2.ORF8,
        'CoV2-RdRp': cov2.RdRp,
        'CoV2-S': cov2.S
        }

    def reference(string):
        if string in references:
            return references[string].load()
        try:
            with open(string) as handle:
                ref = next(SeqIO.parse(handle, 'fasta'))
            return ref
        except:
            msg = "'{0}' does not exist or is not a valid FASTA file".format(string)
            raise ArgumentTypeError(msg)

    kwargs = dict(
        metavar='REFERENCE',
        type=reference,
        help='REFERENCE FASTA file or {{{0}}}'.format(', '.join(references.keys()))
        )

    parser.add_argument(*args, **kwargs)

    return parser


def add_scorematrix(parser, *args):
    import BioExt.scorematrices
    from BioExt.scorematrices import (
        DNAScoreMatrix,
        LazyScoreMatrix,
        ProteinScoreMatrix
        )

    score_matrices = {}
    for obj in dir(BioExt.scorematrices):
        val = getattr(BioExt.scorematrices, obj)
        if isinstance(val, (DNAScoreMatrix, LazyScoreMatrix, ProteinScoreMatrix)):
            score_matrices[str(val)] = val

    kwargs = dict(
        metavar='SCOREMATRIX',
        type=lambda s: score_matrices.get(s, s),
        choices=sorted(score_matrices.values(), key=str),
        default=score_matrices['BLOSUM62'],
        help='parameterize using one of {{{0}}} [default=BLOSUM62]'.format(
            ', '.join(score_matrices.keys())
            )
        )

    parser.add_argument(*args, **kwargs)

    return parser
