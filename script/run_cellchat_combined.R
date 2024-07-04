library(CellChat)
library(patchwork)
library(tidyverse)
library(Seurat)

endo_cluster <- readRDS("data/fixed_annotated_aux.seurat.rds")
#meta_table <- readRDS("data/final_endo_meta_table.RDS")

#Making group_id(s) for cell chat comparison 


run_subset_cell_chat <- function(cell_chat_data,group_id){
  meta_table <- cell_chat_data@meta.data
  cell_chat_data <- GetAssayData(cell_chat_data, assay = "RNA", slot = "data")
  cell_chat_data <- createCellChat(object = cell_chat_data, meta = meta_table, group.by = group_id)
  
  print("Cells start Chatting !")
  
  cell_chat_data@DB <- CellChatDB.human
  
  #pre-processing
  print("pre-prcoessing is now commenced")
  cell_chat_data <- cell_chat_data %>%
    subsetData %>%
    identifyOverExpressedGenes %>%
    identifyOverExpressedInteractions #%>% 
    #projectData(PPI.human) 
  
  print("pre-processsing done -let's start processing!")
  #post-processing
  
  cell_chat_data <- cell_chat_data %>%  
    computeCommunProb(.,raw.use = TRUE,population.size = TRUE) %>%
    filterCommunication(.,min.cells = 10) %>%
    computeCommunProbPathway %>%
    aggregateNet
  

}

make_n_merge_cell_chat <- function(seurat_obj,subset_var,group_id) {
  seurat_obj <- SplitObject(seurat_obj,split.by =subset_var)
  object.list <- lapply(seurat_obj,run_subset_cell_chat,group_id)
  cellchat <- mergeCellChat(object.list, add.names = names(seurat_obj))
  cellchat <- c(cellchat, object.list)                         
}

#select non-published EnS and EnEpi cells 

barcode_filter <- endo_cluster@meta.data %>%
                  filter( (is.na(published_epithelial_subtype) &!is.na(epi_cluster)) | (is.na(published_mesenchymal_subtype) & !is.na(mesen_cluster))) %>%
                  # some EnEpi clusters are reasigned to EnS for reclustering - further filtering is needed 
                  filter(is.na(published_mesenchymal_subtype)) %>% 
                  rownames()
endo_cluster <- endo_cluster %>%
                subset( x=., cells = barcode_filter,invert = TRUE)

endo_cluster@meta.data <-  endo_cluster@meta.data %>% 
                           mutate( harmonized_major_plus_immune_cell_type = str_remove_all(harmonized_major_plus_immune_cell_type, "\\([0-9]\\)$") %>% str_trim() )

cell_chat_data <- run_subset_cell_chat(endo_cluster,"harmonized_major_plus_immune_cell_type")
saveRDS(cell_chat_data,"data/combined/CellChat_all.rds") 

# endo_cluster_subset <- endo_cluster %>%
#                 subset( x=.,subset=LA_Status %in% c("High", "Low")) %>%
#                 subset(x=.,subset= Class== "Peritoneal endometriosis") %>%
#                 subset( x=.,subset= harmonized_major_plus_immune_cell_type != "Exclude")
# 
# cell_chat_data <- make_n_merge_cell_chat(endo_cluster,"LA_Status", "harmonized_major_plus_immune_cell_type")
# saveRDS(cell_chat_data,"data/combined/CellChat_endo_LA_final_exclude.rds")
# 
# endo_cluster_subset <- endo_cluster %>%
# subset( x=.,subset=Class %in% c("Peritoneal endometriosis", "Endometrioma")) %>%
# subset( x=.,subset= harmonized_major_plus_immune_cell_type != "Exclude")
# 
# cell_chat_data <- make_n_merge_cell_chat(endo_cluster_subset,"Class", "harmonized_major_plus_immune_cell_type")
# saveRDS(cell_chat_data,"data/combined/CellChat_endo_class_PE_Endo_exclude.rds")
# 
# endo_cluster_subset <- endo_cluster %>%
#   subset( x=.,subset=Class %in% c( "Endometrioma","Eutopic endometrium")) %>%
#   subset( x=.,subset= harmonized_major_plus_immune_cell_type != "Exclude")
# cell_chat_data <- make_n_merge_cell_chat(endo_cluster_subset,"Class", "harmonized_major_plus_immune_cell_type")
# saveRDS(cell_chat_data,"data/combined/CellChat_endo_class_Endo_Eutopic_exclude.rds")
# 
# endo_cluster_subset <- endo_cluster %>%
#   subset( x=.,subset=Class %in% c( "Peritoneal endometriosis","Eutopic endometrium")) %>%
#   subset( x=.,subset= harmonized_major_plus_immune_cell_type != "Exclude")
# cell_chat_data <- make_n_merge_cell_chat(endo_cluster_subset,"Class", "harmonized_major_plus_immune_cell_type")
# saveRDS(cell_chat_data,"data/combined/CellChat_endo_class_PE_Eutopic_exclude.rds")

# endo_cluster@meta.data <- endo_cluster@meta.data %>%
#                           mutate( Class=fct_recode(Class,"PE+Endo" ="Endometrioma","PE+Endo"="Peritoneal endometriosis"))
# endo_cluster_subset <- endo_cluster %>%
#   subset( x=.,subset=Class %in% c( "PE+Endo","Eutopic endometrium")) %>%
#   subset( x=.,subset= harmonized_major_plus_immune_cell_type != "Exclude")
# cell_chat_data <- make_n_merge_cell_chat(endo_cluster_subset,"Class", "harmonized_major_plus_immune_cell_type")
# saveRDS(cell_chat_data,"data/combined/CellChat_endo_class_Combine_endometrium_exclude.rds")




