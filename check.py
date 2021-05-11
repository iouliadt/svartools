from Bio import SeqIO
import svartools as svt


# Function to create check format table
def checkFormat():
    sequences = []
    seq_IDs = []

    # * Parse fasta
    fastaFile = SeqIO.parse(svt.args.alignment, "fasta")
    for record in fastaFile:
        sequences.append(record.seq)
        seq_IDs.append(record.id)

    # Create checkFormat table and write to file
    if svt.args.outfile:
        fileName = svt.args.outfile
    else:
        fileName = svt.args.alignment

    # Calculate length, count Ns, gaps, codons, IUPAC ambiguous chars
    with open("{}_checkFormat.txt".format(fileName), "w") as f:
        f.write("Sequence\tLength\tCodon_cnt\tN_cnt\tGap_cnt\tAmbi_cnt\n")
        for i in range(len(sequences)):
            IDs = seq_IDs[i]
            length = len(sequences[i])
            Ns = sequences[i].count('N')
            gaps = sequences[i].count('-')
            codons = length/3 - gaps
            # IUPAC ambiguities: W,S,M,K,R,Y,B,D,H,V
            ambiguities = sequences[i].count('W') + sequences[i].count('S')
            + sequences[i].count('M') + sequences[i].count('K')
            + sequences[i].count('R') + sequences[i].count('Y')
            + sequences[i].count('B') + sequences[i].count('D')
            + sequences[i].count('H') + sequences[i].count('V')
            print(IDs, "\t", length, "\t", codons, "\t",
                  Ns, "\t", gaps, "\t", ambiguities, file=f,)
