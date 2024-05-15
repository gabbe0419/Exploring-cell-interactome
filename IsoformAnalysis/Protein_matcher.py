import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def get_uniqueppi_isoforms (interaction_df, proteins_main, proteins_compare, index_main):

    protein_match = protein_matcher(proteins_main, proteins_compare, "canonical_protein", "canonical_protein")
    extracted_isoforms = isoform_extractor(protein_match)
    processed_ppi = process_ppi(interaction_df, extracted_isoforms)
    
    unique_ppi_isoforms = fill_list(index_main, processed_ppi, "isoform", "isoform")
    
    return unique_ppi_isoforms

    
def protein_matcher (dataset, protein_list, dataset_column, protein_list_column, check_duplicates = False):
    
    filtered_df = dataset[dataset[dataset_column].isin(protein_list[protein_list_column].dropna())]
    
    if check_duplicates == True:
        check_duplicates(filtered_df)
        
    return filtered_df

def isoform_extractor ( dataframe):
    isoform_columns = [col for col in dataframe.columns if col.startswith("isoform_")]
    isoform_values = []
    for col in isoform_columns:
        isoform_values.extend(dataframe[col].dropna().unique())
    isoform_df = pd.DataFrame({"isoform": isoform_values})
    isoform_df = isoform_df[~isoform_df['isoform'].astype(str).str.isdigit()]
    unique_isoforms_df = isoform_df.drop_duplicates(subset=["isoform"]).reset_index(drop=True)
    
    return unique_isoforms_df

def process_ppi(dataframe, isoform_df):
    dataframe = dataframe.dropna(subset=['UniprotA', 'UniprotB'])
    
    match_list = []
    unique_isoforms = set(isoform_df['isoform'].dropna())
    for index, row in dataframe.iterrows():
        if row['UniprotA'] in unique_isoforms and row['UniprotB'] in unique_isoforms:
            match_list.append(row)

    dataframe = pd.DataFrame(match_list)
    unique_uniprot_values = pd.concat([dataframe['UniprotA'], dataframe['UniprotB']]).dropna().unique()
    return pd.DataFrame({'isoform': unique_uniprot_values})

def fill_list(dataset, protein_list, dataset_column, protein_list_column, check_duplicates = False):

    dataset_column_1 = "canonical_protein"
    dataset_column_2 = "isoform"
    dataset_column_3 = "gene_name"
    filtered_df = dataset[dataset[dataset_column].isin(protein_list[protein_list_column].dropna())][[dataset_column_1, dataset_column_2, dataset_column_3]]

    if check_duplicates == True:
        check_duplicates(filtered_df)
        
    return filtered_df

def check_duplicates( df):
    duplicates_exist = df.duplicated(subset=["canonical_protein"]).any()
    if duplicates_exist:
        print("DETECTED duplicates in 'canonical_protein'.")
    else:
        print("There're NO duplicates in 'canonical_protein'.")
    return duplicates_exist

    
    


