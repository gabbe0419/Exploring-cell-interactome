# Exploring-cell-interactome

## For isoform analysis
    Run "isoform_analysis.ipynb"

## For intrinsic disorder prediciton and results
    1. Requires installation of Iupred 3 which should be placed into "iupred3/" folder
        Available from: https://iupred.elte.hu/download_new (Downloaded - 2024-03-17 )

    2. Run "disorder_prediction.ipynb" (Uncomment last cell) to generate Iupred3 prediction
        OBS! Long process (ca. 2h)
        
    3. Run "disorder_plotting_output.ipynb" to generate plots

## Data used (And specific file path that code is built upon)
    Data\BioPlex_293T_Network_10K_Dec_2019.tsv ( https://bioplex.hms.harvard.edu  )
    
    Data\BioPlex_HCT116_Network_5.5K_Dec_2019.tsv ( https://bioplex.hms.harvard.edu  )
    
    Data\Huttlin_BioPlex3_Table_S1.xlsb (Huttlin EL, Bruckner RJ, Navarrete-Perea J, Cannon JR, Baltier K, Gebreab F, et al. Dual proteome-scale networks reveal cell-specific remodeling of the human interactome. Cell. 2021 May;184(11):3022–40. doi:10.1016/j.cell.2021.04.011)
    
    Data\uniprotkb_AND_reviewed_true_AND_model_o_2024_02_28.fasta (https://www.uniprot.org/uniprotkb?query=reviewed%3Atrue&facets=reviewed%3Atrue%2Cmodel_organism%3A9606 )
    
    Data\UP000005640_9606_additional.fasta (https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD024364 )

## Environment used in result available in "BB103X-2024-05-15.yml" (python version: 3.11.5):

