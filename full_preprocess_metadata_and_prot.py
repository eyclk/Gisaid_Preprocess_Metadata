import pandas as pd
from Bio import SeqIO

CURRENT_ENV = "./"

SEQ_MIN_LEN = 1235
SEQ_MAX_X_COUNT = 13

metadata_path = CURRENT_ENV + "metadata_tsv_2023_05_22/metadata.tsv"
simplified_metadata_path = CURRENT_ENV + "preprocessed_metadata/simplified_metadata.tsv"
filtered_complete_metadata_path = CURRENT_ENV + "preprocessed_metadata/filtered_all_metadata.tsv"
filtered_specific_variant_path = CURRENT_ENV + "preprocessed_metadata/filtered_metadata_"

all_fasta_prot_file_path = CURRENT_ENV + "hCoV-19_spikeprot0519/spikeprot0519/spikeprot0519.fasta"
complete_filtered_prot_csv = CURRENT_ENV + "preprocessed_spikeprot_csv/complete_filtered_prot_"
sampled_filtered_prot_csv = CURRENT_ENV + "preprocessed_spikeprot_csv/sampled_prot_"

selected_cols = ['Virus name', 'Accession ID', 'Collection date', 'Host', 'Variant', 'Is complete?',
                 'Is high coverage?']
simplified_cols = ['Accession ID', 'Collection date']


def simplify_metadata():
    df = pd.read_table(metadata_path, usecols=selected_cols)
    print("--> Table with selected columns has been read.")
    print("******** Total number of rows =", len(df))

    for i in range(len(df)):
        if isinstance(df.at[i, 'Variant'], str):
            df.at[i, 'Variant'] = df.at[i, 'Variant'].split(' ')[1]
    print("--> Variant explanations simplified to just variant name.")
    df.to_csv(simplified_metadata_path, sep="\t")
    print("--> Simplified table written out to processed_metadata/simplified_metadata.tsv")
    return df


def filter_complete_metadata():
    df = pd.read_table(simplified_metadata_path, usecols=selected_cols)

    indices_to_drop = []
    for i in range(len(df)):
        if str(df.at[i, 'Host']) != 'Human' or (str(df.at[i, 'Is complete?']) != 'True') or (str(df.at[i, 'Is high coverage?']) != 'True') or not((df.at[i, 'Variant'] == 'Omicron') or (df.at[i, 'Variant'] == 'Alpha') or (df.at[i, 'Variant'] == 'Beta') or (df.at[i, 'Variant'] == 'Gamma') or (df.at[i, 'Variant'] == 'Delta')):
            indices_to_drop.append(i)
    df = df.drop(indices_to_drop)

    print("\n\n--> Table filtered out rows for given conditions including for when variant names are something other "
          "than the main 5 variants.")
    df.to_csv(filtered_complete_metadata_path, sep="\t")
    print("--> Filtered complete table written out to processed_metadata/filtered_complete_metadata_path.tsv")
    return df


def get_metadata_for_specific_variant(variant_name='Omicron'):
    df = pd.read_table(filtered_complete_metadata_path, usecols=selected_cols)
    df = df.loc[df['Variant'] == variant_name]

    df.to_csv(filtered_specific_variant_path + variant_name + '.tsv', sep="\t")
    print("\n\n--> Filtered table for", variant_name, "written out to processed_metadata\n")
    return df


# New functions for spikeprots start from here! *****************************************************
def get_spikeprot_for_filtered_metadata(variant_name='Omicron'):
    print("================== For", variant_name, "==================\n")
    df_meta = pd.read_table(filtered_specific_variant_path + variant_name + '.tsv', usecols=simplified_cols)

    date_dict = {}
    for i in range(len(df_meta)):
        date_dict[df_meta.at[i, 'Accession ID']] = df_meta.at[i, 'Collection date']

    id_list = []
    prot_seq_list = []
    prot_fasta = open(all_fasta_prot_file_path, encoding="iso-8859-1")

    fasta_parse_progress = 0
    for f_record in SeqIO.parse(prot_fasta, 'fasta'):
        temp_id_header = f_record.description
        temp_id = temp_id_header.split("|")[3]
        if temp_id[:7] != "EPI_ISL":
            temp_id = temp_id_header.split("|")[2]  # Sometimes, the date is skipped and 3rd element is the id.
        if temp_id[:7] != "EPI_ISL":
            temp_id = temp_id_header.split("|")[4]  # Sometimes, 5th element is the id.

        temp_seq = "".join(f_record.seq[:-1])
        if len(temp_seq) <= SEQ_MIN_LEN or temp_seq.count('X') >= SEQ_MAX_X_COUNT:
            continue

        id_list.append(temp_id)
        prot_seq_list.append(temp_seq)

        if fasta_parse_progress % 1000000 == 0:
            print("--> Fasta file parse progress count =", fasta_parse_progress)
        fasta_parse_progress += 1

    print("\n***** Transforming the contents of fasta into a dataframe:\n")
    df_prot = pd.DataFrame({'accession_id': id_list, 'sequence': prot_seq_list}, columns=['accession_id', 'sequence'])
    print(df_prot)

    df_prot['date'] = df_prot['accession_id'].map(date_dict)

    print("\n***** New id and date columns have been added to the prot dataframe. "
          "Proceeding to drop rows with NaN values.")

    df_prot = df_prot.dropna()
    print("\n***** Smaller prot dataframe with id and date columns (filtered for sequence lengths and X counts):\n")
    print(df_prot, "\n")

    df_prot.to_csv(complete_filtered_prot_csv + variant_name + '.csv', sep=",", index=False)


def try_parsing_date(text):
    for fmt in ('%Y-%m-%d', '%Y-%m', '%Y'):
        try:
            return pd.to_datetime(text, format=fmt)
        except ValueError:
            pass
    print(text, "\n")
    raise ValueError('no valid date format found')


def sort_df_by_date(df):
    for i in range(len(df)):
        df.at[i, 'date'] = try_parsing_date(df.at[i, 'date'])
    df.sort_values(by='date', inplace=True)
    return df


def sample_from_complete_filtered_prot_csv(sample_size=2000, variant_name='Omicron'):
    df_filtered_prot = pd.read_table(complete_filtered_prot_csv + variant_name + '.csv', delimiter=",")
    df_filtered_prot = sort_df_by_date(df_filtered_prot)

    skip_amount = int(len(df_filtered_prot) / sample_size)

    df_sample = df_filtered_prot[::skip_amount]
    df_sample = df_sample[:sample_size]

    df_sample.to_csv(sampled_filtered_prot_csv + variant_name + '_' + str(sample_size) + '.csv', sep=",", index=False)
    print("\n--> Sampled dataset for", variant_name, "written out to processed_spikeprot_csv/sampled_prot_" +
          variant_name + '_' + str(sample_size) + '.csv')


pd.set_option('display.max_columns', None)

"""simplified_df = simplify_metadata()
print(simplified_df)

filtered_complete_df = filter_complete_metadata()
print(filtered_complete_df)

filtered_omicron_df = get_metadata_for_specific_variant(variant_name='Omicron')
print(filtered_omicron_df)"""

"""get_metadata_for_specific_variant(variant_name='Alpha')
get_metadata_for_specific_variant(variant_name='Beta')
get_metadata_for_specific_variant(variant_name='Delta')
get_metadata_for_specific_variant(variant_name='Gamma')"""


get_spikeprot_for_filtered_metadata(variant_name='Beta')
get_spikeprot_for_filtered_metadata(variant_name='Alpha')
get_spikeprot_for_filtered_metadata(variant_name='Delta')
get_spikeprot_for_filtered_metadata(variant_name='Gamma')
get_spikeprot_for_filtered_metadata(variant_name='Omicron')


sample_from_complete_filtered_prot_csv(sample_size=2000, variant_name='Beta')
sample_from_complete_filtered_prot_csv(sample_size=2000, variant_name='Alpha')
sample_from_complete_filtered_prot_csv(sample_size=2000, variant_name='Delta')
sample_from_complete_filtered_prot_csv(sample_size=2000, variant_name='Gamma')
sample_from_complete_filtered_prot_csv(sample_size=2000, variant_name='Omicron')
