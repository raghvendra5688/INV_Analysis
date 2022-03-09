In the scripts folder, just run each R script after installing all the required R packages to obtain each result figure.

Each result figure will appear in Results/Paper_Figures/Figure_i.pdf , where i represents figure number.

Text results are provided in Results/Paper_Text_Results/ folder. Similarly, figures are provided in Results/Paper_Figures/ folder.

The data used for all the experiments is also available at: Mall, Raghvendra (2020), “Transcriptomic Dataset for Network based identification of key Master Regulators for Immunologic Constant of Rejection”, Mendeley Data, V3, doi: 10.17632/d9ffb7kkzt.3

#Data should be downloaded and all its content should be added to the 'Data' folder within main INV_Analysis repository for running R scripts.

In the Data/Others/ folder, please gunzip the me_net_full.Rdata.gz file.

In the Data/PRECOG/ folder, please put the es_list1.rds and es_list2.rds files for performing the validation on PRECOG repository.

Run the misc_figures.R in scripts/ folder for validation results on the PRECOG repository.

Run the INV_Consensus_Classification.R followed by High_Medium_Low_INV_classification.R in scripts/ folder to get the INV-High, INV-Medium and INV-Low labels for each cancer sample per cancer type.
