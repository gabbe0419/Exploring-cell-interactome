import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

    
def get_dataset_count(interaction_df, lib_df): 
    
    unique_isoforms = _read_csv (interaction_df)
    
    iso_lib_match = iso_lib_matcher(unique_isoforms, lib_df)
    can_lib_match = can_lib_matcher(unique_isoforms, lib_df)
    gene_lib_match = gene_lib_matcher( unique_isoforms, lib_df)
    
    lib_match = addition(iso_lib_match, can_lib_match, gene_lib_match, 'isoform')
    undefined_sorted = subtraction(unique_isoforms, lib_match, 'isoform', 'isoform')
    index = addition(lib_match, undefined_sorted, None, 'isoform')
    
    proteins = isoform_count(index)
    
    return proteins, index
    
def isoform_distribution (proteins):
    df = proteins
    
    count_df = df['isoform_count'].value_counts().reset_index()
    count_df.columns = ['isoform_count', '# of entries']
    count_df = count_df.sort_values(by='isoform_count')
    total_entries = count_df['# of entries'].sum()
    count_df['frequency (%)'] = count_df['# of entries'] / total_entries * 100
    
    fig, ax1 = plt.subplots(figsize=(10, 6))

    ax1.bar(count_df['isoform_count'], count_df['frequency (%)'], color='skyblue')
    ax1.set_title('Isoform Distribution')
    ax1.set_xlabel('Number of Isoforms')
    ax1.set_ylabel('Frequency (%)')
    ax1.grid(True)
    ax1.set_xticks(range(int(count_df['isoform_count'].max()) + 1))
    ax1.tick_params(axis='both', which='major')

    zoomed_data = count_df[count_df['frequency (%)'] < 1]
    ax2 = ax1.inset_axes([0.6, 0.6, 0.3, 0.3])  
    ax2.bar(zoomed_data['isoform_count'], zoomed_data['frequency (%)'], color='skyblue')
    ax2.set_xlabel('Number of Isoforms')
    ax2.set_ylabel('Frequency (%)')
    ax2.set_xticks(zoomed_data['isoform_count'])
    ax2.tick_params(axis='both', which='major')
    ax2.set_ylim(0, 1.2 * zoomed_data['frequency (%)'].max())
    
    fig.tight_layout()
    plt.show()
    
    return None

def _read_csv(data_df, check_duplicates = False):

    if "UniprotA" in data_df.columns and "UniprotB" in data_df.columns:
        uniprotA_df = data_df[["UniprotA", "SymbolA"]].dropna().rename(columns={"UniprotA": "isoform", "SymbolA": "gene_name"})
        uniprotB_df = data_df[["UniprotB", "SymbolB"]].dropna().rename(columns={"UniprotB": "isoform", "SymbolB": "gene_name"})
        merged_df = pd.concat([uniprotA_df, uniprotB_df], ignore_index=True)
        unique_isoforms_df = merged_df.drop_duplicates(subset=["isoform"]).dropna()
    
    else:
        uniprotA_293t = data_df[["UniprotA-293T", "SymbolA-293T"]].dropna().rename(columns={"UniprotA-293T": "isoform", "SymbolA-293T": "gene_name"})
        uniprotB_293t = data_df[["UniprotB-293T", "SymbolB-293T"]].dropna().rename(columns={"UniprotB-293T": "isoform", "SymbolB-293T": "gene_name"})
        uniprotA_hct116 = data_df[["UniprotA-HCT116", "SymbolA-HCT116"]].dropna().rename(columns={"UniprotA-HCT116": "isoform", "SymbolA-HCT116": "gene_name"})
        uniprotB_hct116 = data_df[["UniprotB-HCT116", "SymbolB-HCT116"]].dropna().rename(columns={"UniprotB-HCT116": "isoform", "SymbolB-HCT116": "gene_name"})        
        merged_df = pd.concat([uniprotA_293t, uniprotB_293t, uniprotA_hct116, uniprotB_hct116], ignore_index=True)
        unique_isoforms_df = merged_df.drop_duplicates(subset=["isoform"]).dropna()

    if check_duplicates == True:
        check_duplicates(unique_isoforms_df)

    return unique_isoforms_df

def iso_lib_matcher (dataset_df, lib_df, check_duplicates = False):

    lib_column_2 = "unique_identifier"
    lib_column_4 = "isoform_of"
    lib_column_7 = "gene_name"
    dataset_column = "isoform"
    
    filtered_df = lib_df[lib_df[lib_column_2].isin(dataset_df[dataset_column].dropna())][[lib_column_4, lib_column_2, lib_column_7]]
    filtered_df.columns = ["canonical_protein", "isoform", lib_column_7]

    if check_duplicates == True: 
        check_duplicates(filtered_df)

    return filtered_df
    
def can_lib_matcher(dataset_df, lib_df, check_duplicates = False):


    dataset_column_1 = "isoform"
    dataset_column_2 = "gene_name_x"
    lib_column_4 = "isoform_of"
    merged_df = pd.merge(dataset_df, lib_df, how="inner", left_on=dataset_column_1, right_on=lib_column_4)
    filtered_df = merged_df[[lib_column_4, dataset_column_1, dataset_column_2]].drop_duplicates(subset=[dataset_column_1])
    filtered_df.columns = ["canonical_protein", dataset_column_1, "gene_name"]

    if check_duplicates == True: 
        check_duplicates(filtered_df)

    return filtered_df

def gene_lib_matcher (dataset_df, lib_df, check_duplicates = False):

    dataset_column_1 = "isoform"
    dataset_column_2 = "gene_name"
    lib_column_4 = "isoform_of"
    lib_column_7 = "gene_name"
    merged_df = pd.merge(dataset_df, lib_df, how="inner", left_on=dataset_column_2, right_on=lib_column_7)
    filtered_df = merged_df[[lib_column_4, dataset_column_1, dataset_column_2]].drop_duplicates(subset=[dataset_column_1])
    filtered_df.columns = ["canonical_protein", dataset_column_1, dataset_column_2]

    if check_duplicates == True: 
        check_duplicates(filtered_df)

    return filtered_df

def addition(dataset_1, dataset_2, dataset_3, column_drop_duplicates, check_duplicates = False):
    merged_df = pd.concat([dataset_1, dataset_2, dataset_3], ignore_index=True)
    merged_df = merged_df.drop_duplicates(subset=[column_drop_duplicates])

    if check_duplicates == True: 
        check_duplicates(merged_df)
        
    return merged_df

def subtraction (add_df, subtract_df, add_column, subtract_column):
    sum = add_df[~add_df[add_column].isin(subtract_df[subtract_column])]
    
    data_column_1 = "isoform"
    data_column_2 = "gene_name"
    sum_sorted = sum.copy()
    sum_sorted['new_column'] = sum_sorted[data_column_2]
    sum_sorted = sum_sorted.loc[:, ['new_column', data_column_1, data_column_2]]
    sum_sorted.columns = ["canonical_protein", data_column_1, data_column_2]
    
    
    return sum_sorted.sort()    

def isoform_count(dataset, print_total = False):
    isoform_dict = {}

    for index, row in dataset.iterrows():
        data_column_1 = row['canonical_protein']
        data_column_2 = row['isoform']
        if data_column_1 in isoform_dict:
            isoform_dict[data_column_1].append(data_column_2)
        else:
            isoform_dict[data_column_1] = [data_column_2]

    for data_column_1, data_column_2_values in isoform_dict.items():
        for i, data_column_2_value in enumerate(data_column_2_values, start=1):
            dataset.loc[dataset['canonical_protein'] == data_column_1, f'isoform_{i}'] = data_column_2_value

    dataset.drop_duplicates(subset=['canonical_protein'], inplace=True)
    dataset_count = dataset[['canonical_protein'] + [col for col in dataset.columns if col.startswith('isoform_')]].copy()
    dataset_count.reset_index(drop=True, inplace=True)
    dataset_count['isoform_count'] = dataset_count.filter(like='isoform_').notna().sum(axis=1)

    if print_total == True:
        total_isoform_count = dataset_count['isoform_count'].sum()
        print("Total Isoform Count:", total_isoform_count)

    return dataset_count

def check_duplicates (df):
        duplicates_exist = df.duplicated(subset=["isoform"]).any()
        if duplicates_exist:
            print("DETECTED duplicates in 'isoform'.")
        else:
            print("There're NO duplicates in 'isoform'.")
        return duplicates_exist

