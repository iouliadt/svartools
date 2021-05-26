# extract : slice genes/subregions out of a multifasta and write to (new) multifasta file.
# check : check sequence format: count length, Ns, gaps etc.
# count : create count/frequency tables for nucleotides and amino acids.

# Usage:
# python svartools.py extract -gb sequence.gb -aln cog_2020-06-26_alignment_with_ref.fasta [-c 21562:21700] [-o fileName]
# python svartools.py extract -gb sequence.gb -aln cog_2020-06-26_alignment_with_ref.fasta -g S -o [S.fasta]
# python svartools.py extract -gb sequence.gb -aln cog_2020-06-26_alignment_with_ref.fasta -all [-o ]

# python svartools.py check -aln S.fasta [-o fileName]

# python svartools.py count -aa -aln M.fasta -gaps -r NC_045512.2 [-o fileName]
# python svartools.py count -nuc -aln S.fasta -r NC_045512.2[-o fileName]
# python svartools.py count -nuc -aln Egaps.fasta -r NC_045512.2 -g

import parseArgs as pa


def main():

    # svartools extract - custom coords
    if pa.args.command == 'extract':
        import extract
        if pa.args.coordinates:
            extract.extract_customCoords()
        # svartools extract
        else:
            extract.extract_GBCoords()

    # svartools check
    elif pa.args.command == 'check':
        import check
        check.checkFormat()

    # svartools count - create frequency/entropy tables
    elif pa.args.command == 'count' and pa.args.gaps is None:
        import count
        count

    # svartools count - including gaps
    elif pa.args.command == 'count' and pa.args.gaps:
        import countGaps
        countGaps


if __name__ == "__main__":
    main()
