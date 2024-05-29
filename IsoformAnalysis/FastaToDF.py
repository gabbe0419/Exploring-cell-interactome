import pandas as pd
from Bio import SeqIO # type: ignore


def read_fasta( file_path):
    file_path = file_path

    data = []
    for record in SeqIO.parse(file_path, "fasta"):
        header_info = record.id.split("|")
        description = record.description
        data.append({
            "db": header_info[0],
            "unique_identifier": header_info[1],
            "entry_name": header_info[2],
            "isoform_of": extract_isoform_of(description),
            "functional_traits": extract_functional_traits(description),
            "organism_name": extract_organism_name(description),
            "gene_name": extract_gene_name(description),
            "protein_sequence": str(record.seq)
        })

    fasta = pd.DataFrame(data)  
    
    return fasta 

def extract_isoform_of( description: str) -> str:
    """
    Extracts protein of origin from the description.

    :param description: The description part of the fasta entry.
    :type description: str
    :return: Native protein.
    :rtype: str
    """
    isoform_parts = description.split("Isoform of")
    if len(isoform_parts) > 1:
        isoform_info = isoform_parts[1].strip().split(",")[0]
        return isoform_info.strip()
    return None

def extract_functional_traits( description: str) -> str:
    """
    Extracts functional traits from the description.

    :param description: The description part of the fasta entry.
    :type description: str
    :return: Functional traits.
    :rtype: str
    """
    comma_index = description.find(",")
    os_index = description.find("OS=")
    if comma_index != -1 and os_index != -1:
        return description[comma_index + 1:os_index].strip()
    else:
        return None


def extract_organism_name( description: str) -> str:
    """
    Extracts organism name from the description.

    :param description: The description part of the fasta entry.
    :type description: str
    :return: Organism name.
    :rtype: str
    """
    organism_start = description.find("OS=")
    if organism_start != -1:
        organism_end = description.find(" ", organism_start + 3)
        if organism_end != -1:
            organism_end = description.find(" ", organism_end + 1)
            if organism_end == -1:
                organism_end = len(description)
            return description[organism_start + 3:organism_end].strip()
    return None


def extract_gene_name( description: str) -> str:
    """
    Extracts gene name from the description.

    :param description: The description part of the fasta entry.
    :type description: str
    :return: Gene name.
    :rtype: str
    """
    gene_start = description.find("GN=")
    if gene_start != -1:
        gene_end = description.find(" ", gene_start)
        if gene_end == -1:
            gene_end = len(description)
        return description[gene_start + 3:gene_end].strip()
    else:
        return None
