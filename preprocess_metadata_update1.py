import pandas as pd

working_environment_path = "C:/Users/MONSTER/Desktop/GISAID Grup Projesi/Nextstrain için gerekli dosyalar/"
metadata_path = working_environment_path + "metadata_tsv_2023_05_22/metadata.tsv"
simplified_metadata_path = working_environment_path + "preprocessed_metadata/simplified_metadata.tsv"
filtered_complete_metadata_path = working_environment_path + "preprocessed_metadata/filtered_all_metadata.tsv"
filtered_specific_variant_path = working_environment_path + "preprocessed_metadata/"
selected_2000_for_variant_path = working_environment_path + "preprocessed_metadata/2000_samples_per_variant/"

selected_cols = ['Virus name', 'Accession ID', 'Collection date', 'Host', 'Variant', 'Is complete?',
                 'Is high coverage?']


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

    df.to_csv(filtered_specific_variant_path + 'filtered_metadata_' + variant_name + '.tsv', sep="\t")
    print("\n\n--> Filtered table for", variant_name, "written out to processed_metadata\n")
    return df


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
        df.at[i, 'Collection date'] = try_parsing_date(df.at[i, 'Collection date'])
    # df['Collection date'] = pd.to_datetime(df['Collection date'], format="mixed")
    df.sort_values(by='Collection date', inplace=True)
    return df


def get_uniform_2000_sample_for_a_variant(variant_name='Omicron', extra_sample_amount=0):
    df = pd.read_table(filtered_specific_variant_path + 'filtered_metadata_' + variant_name + '.tsv',
                       usecols=selected_cols)
    df = sort_df_by_date(df)

    sample_amount = 2000 + extra_sample_amount
    skip_amount = int(len(df) / sample_amount)
    print("\nskip_amount =", skip_amount)

    df_sample = df[::skip_amount]
    df_sample = df_sample[:sample_amount]

    df_sample.to_csv(selected_2000_for_variant_path + 'sampled_metadata_' + variant_name + '.tsv', sep="\t")
    print("\n\n--> Sampled table for", variant_name, "written out to processed_metadata/2000_samples_per_variant\n")
    return df_sample


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

# sampled_beta_df = get_uniform_2000_sample_for_a_variant(variant_name='Beta')
# print("\n", sampled_beta_df, "\n")

get_uniform_2000_sample_for_a_variant(variant_name='Beta')
get_uniform_2000_sample_for_a_variant(variant_name='Omicron', extra_sample_amount=10)
get_uniform_2000_sample_for_a_variant(variant_name='Alpha')
get_uniform_2000_sample_for_a_variant(variant_name='Delta')
get_uniform_2000_sample_for_a_variant(variant_name='Gamma', extra_sample_amount=3)
