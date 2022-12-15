import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio.pairwise2 import format_alignment
from Bio.Align import substitution_matrices
import numpy as np
from Bio import Align
import sys
aligner = Align.PairwiseAligner()
aligner.mode = 'global'
aligner.open_gap_score = -10
aligner.extend_gap_score = -1
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

# Read in the DNA and protein sequences from the input
dna_sequence = SeqIO.read(sys.argv[1], "fasta")
protein_sequence = SeqIO.read(sys.argv[2], "fasta")
blosum = substitution_matrices.load("BLOSUM62")


# Translate the DNA sequence into its corresponding protein sequence
translated_dna = dna_sequence.translate(table="Standard")
protein_sequence = protein_sequence.seq
translated_dna = translated_dna.seq

alignments = pairwise2.align.localds(translated_dna, protein_sequence, blosum, -10, -1)
globalalignments = pairwise2.align.globalds(translated_dna, protein_sequence, blosum, -10, -1)

# aligner1 = aligner.align(translated_dna, protein_sequence)

# for a in aligner1:
#     print(a)
#     print("Score = %.1f" % a.score)

for a in alignments:
    print(format_alignment(*a))

print("global:")
for a in globalalignments:
    print(format_alignment(*a))



