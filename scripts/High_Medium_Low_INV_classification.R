#################################################################
###
### This script re-classifies patients to High-Medium-Low INV
### clusters as specified in 

### Input- & Outputfiles:
### Cancer,".Normcounts.INV.reps5000/",
### "Data/",Cancer,"/",Cancer, "_INV_cluster_assignment_k2-6.Rdata"
### (Rdata files get updated using this script)
#################################################################

manual_annotation <- function(Cluster_file,Cancer) 
{
  load(Cluster_file)
  
  if (optimal.calinsky == 2){
    table_cluster_assignment$HML_cluster = as.character(table_cluster_assignement$INV_cluster_k2)
    table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "INV1"] = "INV Low"
    table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "INV2"] = "INV High"
  }
  
  if(optimal.calinsky == 3){
    table_cluster_assignment$HML_cluster = as.character(table_cluster_assignment$INV_cluster_k3)
    table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "INV1"] = "INV Low"
    table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "INV2"] = "INV Medium"
    table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "INV3"] = "INV High"
  }
  if(optimal.calinsky == 4){
    table_cluster_assignment$HML_cluster = as.character(table_cluster_assignment$INV_cluster_k4)
    table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "INV1"] = "INV Low"
    table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "INV2"] = "INV Medium"
    table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "INV3"] = "INV Medium"
    table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "INV4"] = "INV High"
  }
  if(optimal.calinsky == 5){
    table_cluster_assignment$HML_cluster = as.character(table_cluster_assignment$INV_cluster_k5)
    table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "INV1"] = "INV Low"
    table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "INV2"] = "INV Medium"
    table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "INV3"] = "INV Medium"
    table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "INV4"] = "INV Medium"
    table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "INV5"] = "INV High"
  }
  return(table_cluster_assignment)
}