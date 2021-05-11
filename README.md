# svartools

#### svartools is a set of commands that process multiFASTA files. There are three main positional arguments: extract, check and count. 

#### Usage:
```
positional arguments:

    extract             Extract regions from multifasta file.
    check               Check sequence format by calculating length, codons, ambiguous characters and gaps.
    count               Determine the frequency of particular mutations relative to a reference genome, 
                        and calculate the entropy, dN dS and polarity of mutations.


extract

usage: svartools.py extract [-h] [-aln ALIGNMENT] [-o OUTFILE] [-gb GENBANK] [-g GENE] [-c COORDINATES]

optional arguments:
  -h, --help            show this help message and exit
  -aln ALIGNMENT, --alignment ALIGNMENT
                        Multifasta alignment file
  -o OUTFILE, --outfile OUTFILE
                        Optional output file name
  -gb GENBANK, --genbank GENBANK
                        Input genbank file
  -g GENE, --gene GENE  The gene(s) to slice
  -c COORDINATES, --coordinates COORDINATES
                        Option to extract a region of the gene


check

usage: svartools.py check [-h] [-aln ALIGNMENT] [-o OUTFILE]

Check sequence format

optional arguments:
  -h, --help            show this help message and exit
  -aln ALIGNMENT, --alignment ALIGNMENT
                        Multifasta alignment file
  -o OUTFILE, --outfile OUTFILE
                        Optional output file name


count

usage: svartools.py count [-h] [-aln ALIGNMENT] [-o OUTFILE] [-r REFERENCE] [-nuc] [-aa]

optional arguments:
  -h, --help            show this help message and exit
  -aln ALIGNMENT, --alignment ALIGNMENT
                        Multifasta alignment file
  -o OUTFILE, --outfile OUTFILE
                        Optional output file name
  -r REFERENCE, --reference REFERENCE
                        Specify the identifier of the reference to use
  -nuc, -nucleotides    Create nucleotide count table
  -aa, -aminoacids      Create amino acid frequency table
(base) 192:svartools juliatsatsani$ python svartools.py check -h
usage: svartools.py check [-h] [-aln ALIGNMENT] [-o OUTFILE]
```

#### Output:

- **extract** 

The output will be a multifasta file with sliced genes, as indicated by the coordinates given in -g. There is also the option for providing custom coordinates using -c. \
Example: 

\>NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome \
ATGGGCTATATAAACGTTTTCGCTTTTCCGTTTACGATATATAGTCTACTCTTGTGCAGA \
ATGAATTCTCGTAACTACATAGCACAAGTAGATGTAGTTAACTTTAATCTCACATAG \
\>Scotland/EDB129/2020 \
ATGGGCTATATAAACGTTTTCGCTTTTCCGTTTACGATATATAGTCTACTCTTGTGCAGA \
ATGAATTCTCGTAACTACATAGCACAAGTAGATGTAGTTAACTTTAATCTCACATAG \
\>England/CAMB-7961F/2020 \
ATGGGCTATATAAACGTTTTCGCTTTTCCGTTTACGATATATAGTCTACTCTTGTGCAGA \
ATGAATTCTCGTAACTACATAGCACAAGTAGATGTAGTTAACTTTAATCTCACATAG 


- **check** 

The output will be a csv file with sequence name, length, codon, N character, gap, and UPAC ambiguous character count columns.\
Example:

| Sequence |	Length  | Codon_cnt|	N_cnt |	Gap_cnt |	Ambi_cnt |
| --- | --- | --- | --- | ---| --- |
| NC_045512.2 |117| 39.0 |0 |	0 |	0 |
| Scotland/EDB129/2020 |117 |39.0 |	0 |	0	| 0 |
| England/CAMB-7961F/2020 | 	117 |	 39.0 |	0 |	0	| 0 |
| England/CAMB-71F1C/2020 |	117	| 39.0 |	0	| 0 |	0 |

- **count** -nuc 

The output will be a csv file with columns that correspond to the reference nucleotide at each site, synonymous/ nonsynonymous nucleotide count, enntropy and polarity of the corresponding aminoacid, so that each line shows how many times a single nucleotide mutation has been found at a site. This command also generates a summary nucfreq table showing the total number of sites where particular substitutions are observed, and a table with the total number of mutations that are in non-synonymous or synonymous codons. \
Example:

- nucfreq table

| NucSite | aa_position | RefNuc | Asyn | Anonsyn | Tsyn | Tnonsyn | Gsyn | Gnonsyn | Csyn | Cnonsyn | entropy | Ref_polarity |
| --- | --- | --- | --- | ---| --- | --- | --- | --- | --- | --- | ---| --- |
| 1 | 1 | A |	3 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | nonpolar  |
| 2 | 1 | T |	0 | 0 | 3 | 0 | 0 | 0 | 0 | 0 | 0 | nonpolar  |
| 3	| 1 | G |	0	| 0	| 0	| 0	| 3 | 0 | 0 | 0 | 0 | nonpolar  |

- nucfreq sum table

| RefNuc | dbNuc |  NbSitesWithChange |
| --- | --- | --- |
| A | C | 1125  |
| A | T | 1125  |
| A | G | 1125  |

- dNdS sum table

| RefNuc | dbNuc |  NonSyn  | Syn |
| --- | --- | --- |--- |
| A | G | 2 | 0 |
| T | K | 1 | 0 |

- **count** -aa 

The output will be a csv file with the mutation, reference amino acid, reference aa polarity, position of the aa, query/alternative aa, alternative polarity, count (how many times the same mutation has been found in the multifasta) and frequency columns. An entropy table is produced as well.  \
Example:

- mutfreq table

| Mutation |	Ref_aa	| Ref_polarity	| aa_position	| Query_aa | Alt_polarity	| frequency	| count |
| --- | --- | --- | --- | ---| --- | --- | --- |
| M1L	| M	| nonpolar	|1	| L	| nonpolar	| 3,82E+10	| 1 |
| M1X	| M	| nonpolar	|1	| X	| ambiguous	| 0.0012979576255010498	| 34 |
| G2X	| G	| nonpolar	|2	| X	| ambiguous	| 0.0012216071769421645	| 32 |


- entropy table

| aa_position | Ref_aa  | entropy | entropy(base 2) | nb_nonsyn | nb_syn  | dN/dS|
| --- | --- | --- | --- | ---| --- | --- |
| 614 | D | 0.09010335735736986 | 0.38997500048077083 | 2 | 0 | inf|
| 1214  | W | 0.1220680320742344  |  0.5283208335737187 | 1 | 0 | inf|
