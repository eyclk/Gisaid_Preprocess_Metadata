import sys

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from preprocess_metadata import get_uniform_2000_sample_for_a_variant

SEQ_MIN_LEN = 1235
SEQ_MAX_X_COUNT = 13

CURRENT_ENV = "C:/Users/MONSTER/Desktop/GISAID Grup Projesi/Nextstrain için gerekli dosyalar/"

sampled_metadata_files_path = CURRENT_ENV + "preprocessed_metadata/2000_samples_per_variant/"
fasta_prot_sample_output_path = CURRENT_ENV + "preprocessed_spikeprot/2000_sample_spikeprot_"
all_fasta_prot_file_path = CURRENT_ENV + "hCoV-19_spikeprot0519/spikeprot0519/spikeprot0519.fasta"
filtered_specific_variant_path = CURRENT_ENV + "preprocessed_metadata/"
spike_prot_all_filtered_path = CURRENT_ENV + \
                               "preprocessed_spikeprot/filtered_spikeprot_for_length_and_X/all_filtered_spikeprot_"

selected_cols = ['Virus name', 'Accession ID']


def acquire_all_spike_prot_for_a_variant(variant_name="Omicron"):
    metadata_df = pd.read_table(filtered_specific_variant_path + "filtered_metadata_" + variant_name + ".tsv",
                                usecols=selected_cols)
    metadata_ids = metadata_df['Accession ID'].tolist()

    desired_prot_list = []
    fasta = open(all_fasta_prot_file_path, encoding="iso-8859-1")

    for f_record in SeqIO.parse(fasta, 'fasta'):
        temp_id_header = f_record.description
        seq = f_record.seq[:-1]
        temp_id = temp_id_header.split("|")[3]

        if len(seq) <= SEQ_MIN_LEN or seq.count('X') >= SEQ_MAX_X_COUNT:
            continue

        if temp_id[:7] != "EPI_ISL":
            temp_id = temp_id_header.split("|")[2]  # Sometimes, the date is skipped and 3rd element is the id.
        if temp_id[:7] != "EPI_ISL":
            temp_id = temp_id_header.split("|")[4]  # Sometimes, 5th element is the id.

        temp_id_header = temp_id_header.replace('\x9c', '-').replace('\x87', '-')
        if temp_id in metadata_ids:
            desired_prot_list.append(SeqRecord(seq, id=temp_id_header, description=""))
            metadata_ids.remove(temp_id)

    print("\n=========> Number of appropriate spike prots for the variant " + variant_name + ": "
          + str(len(desired_prot_list)) + ".\n")
    output_fasta = open(spike_prot_all_filtered_path + variant_name + ".fasta", "w", encoding="iso-8859-1")
    SeqIO.write(desired_prot_list, spike_prot_all_filtered_path + variant_name + ".fasta", "fasta")
    fasta.close()
    output_fasta.close()


def acquire_sample_fasta_prot_for_variant(variant_name="Omicron"):
    metadata_df = pd.read_table(sampled_metadata_files_path + "sampled_metadata_" + variant_name + ".tsv",
                                usecols=selected_cols)
    metadata_sample_ids = metadata_df['Accession ID'].tolist()
    print("\n***** For", variant_name, "*****")
    print("metadata_sample_ids[:5] --> ", metadata_sample_ids[:5], "\n")

    desired_prot_list = []
    fasta = open(all_fasta_prot_file_path, encoding="iso-8859-1")

    filtered_out_count = 0

    found_count = 0
    for f_record in SeqIO.parse(fasta, 'fasta'):
        temp_id_header = f_record.description
        seq = f_record.seq[:-1]
        temp_id = temp_id_header.split("|")[3]

        if temp_id[:7] != "EPI_ISL":
            temp_id = temp_id_header.split("|")[2]  # Sometimes, the date is skipped and 3rd element is the id.
        if temp_id[:7] != "EPI_ISL":
            temp_id = temp_id_header.split("|")[4]  # Sometimes, 5th element is the id.

        temp_id_header = temp_id_header.replace('\x9c', '-').replace('\x87', '-')
        if temp_id in metadata_sample_ids:

            # **********************************************
            if len(seq) <= SEQ_MIN_LEN or seq.count('X') >= SEQ_MAX_X_COUNT:
                filtered_out_count += 1
                continue

            found_count += 1
            desired_prot_list.append(SeqRecord(seq, id=temp_id_header, description=""))
            metadata_sample_ids.remove(temp_id)
        if found_count == 2000:
            break

    print("======> Remaining unused sample count: ", len(metadata_sample_ids) - filtered_out_count)

    if len(desired_prot_list) != 2000:
        print("\n\nERROR! There are less than 2000 spikeprots found. Current count -->", len(desired_prot_list))
        print("\nRemaining metadata_sample_ids ==> ", metadata_sample_ids)
        sys.exit(-1)

    output_fasta = open(fasta_prot_sample_output_path+variant_name+".fasta", "w", encoding="iso-8859-1")
    SeqIO.write(desired_prot_list, fasta_prot_sample_output_path+variant_name+".fasta", "fasta")
    fasta.close()
    output_fasta.close()


print("\n~~~~~~ For Alpha ~~~~~~\n")
get_uniform_2000_sample_for_a_variant(variant_name='Alpha', extra_sample_amount=100)
acquire_sample_fasta_prot_for_variant(variant_name="Alpha")

print("\n~~~~~~ For Beta ~~~~~~\n")
get_uniform_2000_sample_for_a_variant(variant_name='Beta', extra_sample_amount=2400)
acquire_sample_fasta_prot_for_variant(variant_name="Beta")

print("\n~~~~~~ For Delta ~~~~~~\n")
get_uniform_2000_sample_for_a_variant(variant_name='Delta', extra_sample_amount=500)
acquire_sample_fasta_prot_for_variant(variant_name="Delta")

print("\n~~~~~~ For Gamma ~~~~~~\n")
get_uniform_2000_sample_for_a_variant(variant_name='Gamma', extra_sample_amount=60)
acquire_sample_fasta_prot_for_variant(variant_name="Gamma")

print("\n~~~~~~ For Omicron ~~~~~~\n")
get_uniform_2000_sample_for_a_variant(variant_name='Omicron', extra_sample_amount=400)
acquire_sample_fasta_prot_for_variant(variant_name="Omicron")
