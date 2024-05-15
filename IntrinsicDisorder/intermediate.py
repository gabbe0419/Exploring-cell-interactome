import pandas
import numpy
from Bio import SeqIO # type: ignore
import subprocess
from concurrent.futures import ThreadPoolExecutor
import os


class Intermediate:
    def __init__(self, input, write=False):
        self.data_293T = pandas.read_csv(
            input[0], sep="\t", header=0
        )  # Requires correct file path for bioplex3 293T-newtork
        self.data_HCT116 = pandas.read_csv(
            input[1], sep="\t", header=0
        )  # Requires correct file path for bioplex3 HCT116-newtork

        # This 'xlsb' file type require pip install of something, read the error message
        data_merged_genes = pandas.read_excel(
            input[2],
            sheet_name="E. Merged Interaction List",
        )  # Requires correct file path for bioplex3 Table S1
        # Removes unnescesary columns
        self.data_merged = data_merged_genes.drop(
            columns=[
                "SymbolA",
                "SymbolB",
                "293T (2%)",
                "293T (5%)",
                "HCT116 (2%)",
                "HCT116 (5%)",
            ],
        ).copy()

        self.unique_293T, self.unique_HCT116, self.shared = self.filter_interactions()

        self.proteins_293T, self.proteins_HCT116, self.proteins_shared = (
            self.unique_proteins()
        )

        self.sequence_df = self.read_fasta(input[3])

        if write == True:
            # Writes result to csv
            self.unique_293T.to_csv(
                "Data/Interactions/Unique_interactions_293T.csv", index=False
            )
            self.unique_HCT116.to_csv(
                "Data/Interactions/Unique_interactions_HCT116.csv", index=False
            )
            self.shared.to_csv("Data/Interactions/Shared_interactions.csv", index=False)

    def filter_interactions(self):
        data_293T = self.data_293T
        data_HCT116 = self.data_HCT116
        data_merged_genes = self.data_merged

        # Function to sort bioplex-data by gene that maintans the relation of gene-uniprot-symbol identifiers
        def sort_by_gene(df):
            df_sort = numpy.sort(df[["GeneA", "GeneB"]])
            filt = df[["GeneA", "GeneB"]].ne(df_sort, axis=1)
            data_flipped = numpy.array(
                df.loc[
                    filt["GeneA"],
                    ["GeneB", "GeneA", "UniprotB", "UniprotA", "SymbolB", "SymbolA"],
                ]
            )
            data_unflipped = numpy.array(
                df.loc[
                    filt["GeneA"] == False,
                    ["GeneA", "GeneB", "UniprotA", "UniprotB", "SymbolA", "SymbolB"],
                ]
            )
            return numpy.concatenate([data_unflipped, data_flipped])

        data_293T[["GeneA", "GeneB", "UniprotA", "UniprotB", "SymbolA", "SymbolB"]] = (
            sort_by_gene(data_293T)
        )
        data_HCT116[
            ["GeneA", "GeneB", "UniprotA", "UniprotB", "SymbolA", "SymbolB"]
        ] = sort_by_gene(data_HCT116)

        # Sorts genes in merged Data
        data_merged_genes[["GeneA", "GeneB"]] = numpy.sort(
            data_merged_genes[["GeneA", "GeneB"]]
        )

        # Filters specified database by genes and 'shared-boolean' in merged-interaction data
        def filter_by_gene(df, trueiness):
            data = pandas.merge(
                left=df, right=data_merged_genes, on=["GeneA", "GeneB"], how="inner"
            )
            data = data[data["Shared"] == trueiness]
            return data

        unique_293T = filter_by_gene(data_293T, False).replace("UNKNOWN", numpy.nan)
        unique_HCT116 = filter_by_gene(data_HCT116, False).replace("UNKNOWN", numpy.nan)

        # Combines the shared interactions from both bioplex-data, since not all shared interactionsare present in both original datasets
        # This is further complicated by not all uniprot labels being equal for genes of both cell-lines (They are listed as their isoforms)
        shared_293T = (
            filter_by_gene(data_293T, True)
            .drop(columns=["pW", "pNI", "pInt", "Shared"])
            .replace("UNKNOWN", numpy.nan)
        )
        shared_HCT116 = (
            filter_by_gene(data_HCT116, True)
            .drop(columns=["pW", "pNI", "pInt", "Shared"])
            .replace("UNKNOWN", numpy.nan)
        )
        shared_interactions = pandas.merge(
            left=shared_293T,
            right=shared_HCT116,
            on=["GeneA", "GeneB"],
            how="outer",
            suffixes=["-293T", "-HCT116"],
        ).replace("UNKNOWN", numpy.nan)

        return unique_293T, unique_HCT116, shared_interactions

    def read_fasta(self, input):
        identifiers = []
        sequences = []
        fasta_text = []

        for seq_record in SeqIO.parse(input, "fasta"):
            i = 0
            id = ""
            found_pipe = False
            for a in seq_record.id:
                if a == "|" and not found_pipe:
                    found_pipe = True
                    continue
                elif a == "|" and found_pipe:
                    break
                id += a
            identifiers.append(id.replace("sp", ""))
            sequences.append(seq_record.seq)
            fasta_text.append(f">{seq_record.id}\n{seq_record.seq}")

        sequence_dict = {"ID": identifiers, "Sequence": sequences, "Fasta": fasta_text}
        return pandas.DataFrame(sequence_dict)

    def unique_proteins(self):
        # Lists all unique proteins in bioplex3 data
        proteins_unique_293T = pandas.concat(
            [self.unique_293T["UniprotA"], self.unique_293T["UniprotB"]]
        ).drop_duplicates(keep="first")
        proteins_unique_HCT116 = pandas.concat(
            [self.unique_HCT116["UniprotA"], self.unique_HCT116["UniprotB"]]
        ).drop_duplicates(keep="first")
        proteins_shared = pandas.concat(
            [
                self.shared["UniprotA-293T"],
                self.shared["UniprotB-293T"],
                self.shared["UniprotB-HCT116"],
                self.shared["UniprotB-HCT116"],
            ]
        ).drop_duplicates(keep="first")

        return proteins_unique_293T, proteins_unique_HCT116, proteins_shared

    def protein_grouping(self):
        unique_293T = set(self.proteins_293T)
        unique_HCT116 = set(self.proteins_HCT116)
        shared = set(self.proteins_shared)
        all_proteins = (
            pandas.concat(
                [self.proteins_293T, self.proteins_HCT116, self.proteins_shared]
            )
            .drop_duplicates()
            .dropna()
        )

        # Proteins appearing in unique interactions of both cell lines
        uniqueppi_allcells_inclusive = unique_293T.union(unique_HCT116)

        # Proteins appearing in both unique AND shared interactions of both cell lines
        allppi_allcells_inclusive = uniqueppi_allcells_inclusive.intersection(shared)

        # Proteins only appearing in unique interactions of botch cell lines
        uniqueppi_allcells_exclusive = (
            uniqueppi_allcells_inclusive - allppi_allcells_inclusive
        )

        # Proteins only appearing in shared interactions
        sharedppi_allcells_exclusive = shared - allppi_allcells_inclusive

        # Proteins in unique interactions also unique to their respective cell line
        uniqueppi_293T_exclusive = uniqueppi_allcells_exclusive - unique_HCT116
        uniqueppi_HCT116_exclusive = uniqueppi_allcells_exclusive - unique_293T

        # Proetins shared between cell lines but in unique interactions
        uniqueppi_common = (
            uniqueppi_allcells_exclusive
            - uniqueppi_293T_exclusive
            - uniqueppi_HCT116_exclusive
        )

        protein_lists = [
            uniqueppi_allcells_inclusive,
            shared,
            sharedppi_allcells_exclusive,
            allppi_allcells_inclusive,
            uniqueppi_common,
            uniqueppi_293T_exclusive,
            uniqueppi_HCT116_exclusive,
        ]

        lists = [
            "proteins_uniqueppi",
            "proteins_sharedppi",
            "proteins_common_sharedppi",
            "proteins_common_allppi",
            "proteins_common_uniqueppi",
            "proteins_unique_uniqueppi_293T",
            "proteins_unique_uniqueppi_HCT116",
        ]

        for i, prot_list in enumerate(protein_lists):
            filt = all_proteins.isin(prot_list)
            proteins = all_proteins.loc[filt]
            proteins.to_csv(f"Data/Protein_lists/{lists[i]}.csv")
