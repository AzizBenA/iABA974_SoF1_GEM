from Bio import SeqIO

# you can change the path and name of the input and output file.
gbff_filename = r'output\merged_gbff.gbff'
faa_filename = r'output\FAA_merged.faa'

input_handle = open(gbff_filename, "r")
output_handle = open(faa_filename, "w")
# Biopython is very valuable package in order to handle several biological files for instance GBFF.
for seq_record in SeqIO.parse(input_handle, "genbank"):
    print(f"Dealing with GenBank record {seq_record.id}")
    for seq_feature in seq_record.features:
        # for loop to iterate over each feature of the GBFF file.
        if seq_feature.type == "CDS":
            locus_tag = seq_feature.qualifiers['locus_tag'][0]  # This variable is to get the locus tag
            product = seq_feature.qualifiers['product'][0] if 'product' in seq_feature.qualifiers else "Unknown"  # This variable is to get the gene product
            translation = seq_feature.qualifiers['translation'][0] if 'translation' in seq_feature.qualifiers else ""  # This variable is to get the AA sequence
            output_handle.write(f">{locus_tag} {product}\n{translation}\n") # writing over the fasta file

output_handle.close()
input_handle.close()
