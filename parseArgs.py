# extract : slice genes/subregions out of a multifasta and write to (new) multifasta file.
# check : check sequence format: count length, Ns, gaps etc.
# count : create count/frequency tables for nucleotides and amino acids.

# Usage:
# python svartools.py extract -gb sequence.gb -aln cog_2020-06-26_alignment_with_ref.fasta [-c 21562:21700] [-o fileName]
# python svartools.py extract -gb sequence.gb -aln cog_2020-06-26_alignment_with_ref.fasta -g S -o [S.fasta]
# python svartools.py extract -gb sequence.gb -aln cog_2020-06-26_alignment_with_ref.fasta -all [-o ]

# python svartools.py check -aln S.fasta [-o fileName]

# python svartools.py count -aa -aln S.fasta -r NC_045512.2 [-o fileName]
# python svartools.py count -nuc -aln S.fasta -r NC_045512.2[-o fileName]
# python svartools.py count -nuc -aln Egaps.fasta -r NC_045512.2 -g


import argparse
import sys

# Parse command line arguments
parser = argparse.ArgumentParser(
    description='Extract regions, check sequence format or create '
    'nucleotide/amino acid count tables from a multifasta file.')

files = argparse.ArgumentParser(add_help=False)
files.add_argument('-aln', '--alignment',
                   help='Multifasta alignment file')
files.add_argument('-o', '--outfile',
                   help='Optional output file name')

subparsers = parser.add_subparsers(help='sub-command help', dest='command')

# Create subparsers for each positional argument

# extract parser
parser_extract = subparsers.add_parser('extract',
                                       help='Extract regions '
                                       'from multifasta file.',
                                       parents=[files])
parser_extract.add_argument('-gb', '--genbank',
                            help='Input genbank file')
parser_extract.add_argument('-g', '--gene',
                            help='The gene(s) to slice')
# parser_extract.add_argument('-all', '--allgenes',
#                            help='All gene names')
parser_extract.add_argument('-c', '--coordinates',
                            help='Option to extract a region of the gene')

# check parser
parser_check = subparsers.add_parser('check', description='Check sequence format',
                                     help='Create check format table',
                                     parents=[files])

# count parser
parser_count = subparsers.add_parser('count',
                                     help='Create nucleotide or amino acid '
                                     'frequency tables from a multifasta file',
                                     parents=[files])
parser_count.add_argument('-r', '--reference',
                          help='Specify the identifier of the reference to use')
parser_count.add_argument('-gaps', '--gaps', action='store_true',
                          help='Option to incude gaps')
parser_count.add_argument('-nuc', '--nucleotides', action='store_true',
                          help='Create nucleotide count table')
parser_count.add_argument('-aa', '--aminoacids', action='store_true',
                          help='Create amino acid frequency table')


# Display help message when svartools is called with wrong arguments
args = parser.parse_args()

if not vars(args) or len(sys.argv) <= 2:
    parser.print_help()
    parser.exit(1)

input_file = open(args.alignment)
