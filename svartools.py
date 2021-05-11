# extract : slice genes/subregions out of a multifasta and write to (new) multifasta file.
# check : check sequence format: count length, Ns, gaps etc.
# count : create count/frequency tables, determine the entropy, dN dS and polarity of mutations.

# Usage:
# python svartools.py extract -gb sequence.gb -aln cog_2020-06-26_alignment_with_ref.fasta [-c 21562:21700] [-o fileName]
# python svartools.py extract -gb sequence.gb -aln cog_2020-06-26_alignment_with_ref.fasta -g S -o [S.fasta]
# python svartools.py extract -gb sequence.gb -aln cog_2020-06-26_alignment_with_ref.fasta -all [-o ]

# python svartools.py check -aln S.fasta [-o fileName]

# python svartools.py count -aa -aln S.fasta -r NC_045512.2 [-o fileName]
# python svartools.py count -nuc -aln S.fasta -r NC_045512.2[-o fileName]


import argparse
import sys
import count
import extract
import check


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
parser_extract.add_argument('-c', '--coordinates',
                            help='Option to extract a region of the gene')

# check parser
parser_check = subparsers.add_parser('check', description='Check sequence format',
                                     help='Check sequence format by calculating length,'
                                     ' codons, ambiguous characters and gaps.',
                                     parents=[files])

# count parser
parser_count = subparsers.add_parser('count',
                                     help='Determine the frequency of particular mutations'
                                     ' relative to a reference genome,'
                                     ' and calculate the entropy, dN dS and polarity of mutations.',
                                     parents=[files])
parser_count.add_argument('-r', '--reference',
                          help='Specify the identifier of the reference to use')
parser_count.add_argument('-nuc', '-nucleotides', action='store_true',
                          help='Create nucleotide count table')
parser_count.add_argument('-aa', '-aminoacids', action='store_true',
                          help='Create amino acid frequency table')


# Display help message when svartools is called with wrong arguments
args = parser.parse_args()
if not vars(args) or len(sys.argv) <= 2:
    parser.print_help()
    parser.exit(1)

# svartools extract - custom coords
if args.command == 'extract' and args.coordinates:
    extract.extract_customCoords()


# svartools extract
elif args.command == 'extract':
    extract.extract_GBCoords()


# svartools check
elif args.command == 'check':
    check.checkFormat()


# svartools count - create frequency/entropy tables
elif args.command == 'count':
    count
