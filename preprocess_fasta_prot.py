import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

working_environment_path = "C:/Users/MONSTER/Desktop/GISAID Grup Projesi/Nextstrain için gerekli dosyalar/"
sampled_metadata_files_path = working_environment_path + "preprocessed_metadata/2000_samples_per_variant/"
fasta_prot_sample_output_path = working_environment_path + "preprocessed_spikeprot/2000_sample_spikeprot_"
all_fasta_prot_file_path = working_environment_path + "hCoV-19_spikeprot0519/spikeprot0519/spikeprot0519.fasta"

selected_cols = ['Virus name', 'Accession ID']


def acquire_sample_fasta_prot_for_variant(variant_name="Omicron"):
    metadata_df = pd.read_table(sampled_metadata_files_path + "sampled_metadata_" + variant_name + ".tsv",
                                usecols=selected_cols)
    metadata_sample_ids = metadata_df['Accession ID'].tolist()
    print("\n***** For", variant_name, "*****")
    print("metadata_sample_ids[:5] --> ", metadata_sample_ids[:5], "\n")

    desired_prot_list = []
    fasta = open(all_fasta_prot_file_path, encoding="iso-8859-1")

    found_count = 0
    for f_record in SeqIO.parse(fasta, 'fasta'):
        temp_id_header = f_record.description
        # temp_id_header = f_record.id + " " + " ".join(f_record.description.split(" ")[1:])
        seq = f_record.seq
        temp_id = temp_id_header.split("|")[3]

        if temp_id[:7] != "EPI_ISL":
            # print("temp_id_header --->>>", temp_id_header)
            temp_id = temp_id_header.split("|")[2]  # Sometimes, the date is skipped and 3rd element is the id.
        if temp_id[:7] != "EPI_ISL":
            temp_id = temp_id_header.split("|")[4]  # Sometimes, 5th element is the id.
            # print("temp_id --> ", temp_id, "\n")
        if temp_id[:7] != "EPI_ISL":
            print("**** There are still incorrectly selected ids from the header. ****")

        if temp_id in metadata_sample_ids:
            found_count += 1
            desired_prot_list.append(SeqRecord(seq, id=temp_id_header.replace('\x9c', '-'), description=""))
            metadata_sample_ids.remove(temp_id)
        if found_count == 2000:
            break

    if len(desired_prot_list) != 2000:
        print("\n\nERROR! There are less than 2000 spikeprots found. Current count -->", len(desired_prot_list))
        print("\nRemaining metadata_sample_ids ==> ", metadata_sample_ids)

    output_fasta = open(fasta_prot_sample_output_path+variant_name+".fasta", "w", encoding="iso-8859-1")
    SeqIO.write(desired_prot_list, fasta_prot_sample_output_path+variant_name+".fasta", "fasta")
    fasta.close()
    output_fasta.close()


acquire_sample_fasta_prot_for_variant(variant_name="Alpha")
acquire_sample_fasta_prot_for_variant(variant_name="Beta")
acquire_sample_fasta_prot_for_variant(variant_name="Delta")

acquire_sample_fasta_prot_for_variant(variant_name="Gamma")
acquire_sample_fasta_prot_for_variant(variant_name="Omicron")
