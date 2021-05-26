from Bio.SeqFeature import FeatureLocation, CompoundLocation
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
import parseArgs as pa


# Slice out fasta using coordinates from gb file - create a file for each set of coords
def extract_GBCoords():

    coordinates = {}
    # * Parse genbank file
    gbFile = SeqIO.parse(pa.args.genbank, "genbank")
    for record in gbFile:
        # Find eference name
        gbRef = record.id  # NC_045512
        for f in record.features:
            # Find coding regions and gene/id to populate dict
            if f.type == "CDS" and "gene" in f.qualifiers:
                i = 1
                # If key already exists in the dict,
                if f.qualifiers["gene"][0] in coordinates.keys():
                    # add "_1"
                    f.qualifiers["gene"][0] += '_' + str(i)
                    i += 1
                coordinates[f.qualifiers["gene"][0]] = f.location

    # * Open alignment file
    input_file = pa.input_file
    # Write to dictionary
    fasta_dict = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"))

    # Check if the genbank ID exists in the alignment
    for id, seq in fasta_dict.items():
        if gbRef == id:
            refSeq = seq

    if gbRef not in fasta_dict.keys():
        print("{} was not found in {}".format(gbRef, pa.args.alignment))
        exit()

    # Create 'positions' dict that associates the original coords with the adjusted ones
    pos = 0
    positions = {}

    for i, nuc in enumerate(refSeq):
        if nuc != "-":
            positions[pos] = i
            pos += 1
    # print(positions)

    new_coords = {}
    loc = []

    for key in coordinates:

        if pa.args.gene in key:
            # check if coordinates are joined (eg for ORF1ab)
            if type(coordinates[key]) == CompoundLocation:

                for i, f in enumerate(coordinates[key].parts):

                    loc.append(FeatureLocation(positions[coordinates[key].parts[i].start],
                                               positions[coordinates[key].parts[i].end]))
                new_coords[key] = CompoundLocation(loc)

                print("Extracting gene {}, region {}".format(pa.args.gene,
                                                               new_coords[key]))

            else:

                new_coords[key] = FeatureLocation(positions[coordinates[key].parts[0].start],
                                                  positions[coordinates[key].parts[0].end])

                print("Extracting gene {}, region {}".format(pa.args.gene,
                                                               new_coords[key]))

    for key in new_coords:
        # Custom file name provided by the user
        if pa.args.outfile:
            fileName = "{}.fasta".format(pa.args.outfile)
        # Default file name
        else:
            fileName = '{}.fasta'.format(key)

        print("Writing to file {}".format(fileName))

        # Write to file
        with open(fileName, "w") as f:
            for sequence in fasta_dict.values():
                subsequence = new_coords[key].extract(sequence)
                SeqIO.write(subsequence, f, "fasta")


# Slice out fasta using custom coords (-c)
def extract_customCoords():

    # Split gene coords
    coords = pa.args.coordinates.split(':')
    start = int(coords[0])
    end = int(coords[1])

    # Create list to append extracted sequences
    sliced_genes = []
    # * Parse fasta
    fastaFile = SeqIO.parse(pa.args.alignment, "fasta")
    # Use coordinates to extract subsequences from each fasta record
    for record in fastaFile:
        sliced_genes.append(record[start:end])
    # Create MultipleSeqAlignment object
    alignment = MultipleSeqAlignment(sliced_genes)

    if pa.args.outfile:
        # Custom file name
        fileName = "{}.fasta".format(pa.args.outfile)
    else:
        # Default file name
        fileName = '{}_{}.fasta'.format(start, end)

    # Write to file
    with open(fileName, "w") as f:
        for sequences in alignment:
            SeqIO.write(sequences, f, "fasta")
