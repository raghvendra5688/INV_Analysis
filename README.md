In the scripts folder, just run each R script after installing all the required R packages to obtain each result figure.

Each result figure will appear in Results/Revised_Figures/Revised_Figure_(v2)_i.pdf , where i represents figure number.

Text results are provided in Results/Revised_Text_Results/ folder. Similarly, figures are provided in Results/Revised_Figures/ folder.

The data used for all the experiments is also available at: Mall, Raghvendra (2020), “Transcriptomic Dataset for Network based identification of key Master Regulators for Immunologic Constant of Rejection”, Mendeley Data, V3, doi: 10.17632/d9ffb7kkzt.3

#Data should be downloaded and all its content should be added to the 'Data' folder within main ICR_Analysis repository for running R scripts.

In the Data/Others/ folder, please unzip the me_net_full.Rdata.gz file.

In the Data/PRECOG/ folder, please put the es_list1.rds and es_list2.rds files (where are these files present/how do we get the said files?) for performing the validation on PRECOG repository.

Run the misc_figures.R in scripts/ folder for validation results on the PRECOG repository.

Run the ICR_Consensus_Classification.R followed by High_Medium_Low_ICR_classification.R in scripts/ folder to get the ICR-High, ICR-Medium and ICR-Low labels for each cancer sample per cancer type.
