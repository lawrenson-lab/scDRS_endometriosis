####### RNA-seq analysis on fresh epithelial and stromal cells #####

##### IL-1 gene expression for M2 paper ########

#### Previous PCA analysis and covariate performed in Fresh only analysis - leading to the exclusion of three samples ###

#### Load the appropriate package and set the Working directory ####
library(dplyr)
library(edgeR)
library(Seurat)
library(ggplot2)
library(lsr)
library(limma)
library(tidyverse)
library(ggrepel)
library(magrittr)

setwd("/Volumes/ENDOOMICBM-Q4641")


##### Load the expression data ####
counts <- read.csv("/Volumes/ENDOOMICBM-Q4641/M2_paper/ins/IMB_Organoid_Cell_gene_count_matrix_100bp_ENSEMBL_Dec2021.csv", header=T, check.names = F)
head(counts)
rownames(counts) <- counts[,1]
counts[,1] <- NULL


###################
# LOWLY EXPRESSED GENES FILTERING 
thres <- as.numeric(cpm(10,mean(colSums(counts)))) # check the average cpm to get 10 counts in each gene (on average)
keep = rowSums(cpm(counts) > thres) >= (0.9*ncol(counts)) # keep only genes with min 10 counts in at least N samples 
summary(keep)
counts = counts[keep,] # CPM-filtered transcripts

counts$ID <- rownames(counts)
counts$ID<-gsub("\\|.*","",counts$ID) #this helped to remove the string after the |

#True = 18,998
#False = 39,304

##### Load the sample info ####
info <- read.csv("/Volumes/ENDOOMICBM-Q4641/M2_paper/ins/Sample_Information.csv", header=T)  
index <- which(info$ID %in% colnames(counts))  
info <- info[index,]
info <- info[order(info$Cell_Type_2),]

##### Load gene info
gene_lists <- read.csv("/Volumes/ENDOOMICBM-Q4641/M2_paper/ins/gene_lists.csv", header = T)
gene_symbol=read.csv('/Volumes/ENDOOMICBM-Q4641/M2_paper/ins/Gene_information_ENSEMBL.csv',
                     stringsAsFactors = F)
gene_symbol=gene_symbol[,6:7]

#### Subset out fresh samples

index <- which(info$Passage == "F")
info2 <- info[index,]

# Get the IDs from the info2 data frame
ids_to_select <- info2$ID

# Select columns from count dataframe based on IDs from info dataframe
fresh <- counts[, c("ID", ids_to_select)]


### Get rid of the first column (with Gene ID and Gene name)
rownames(fresh) <- fresh[,1]
fresh[,1] <- NULL
class(fresh)
fresh <-as.matrix(fresh)
#fresh$ID <- rownames(fresh)


##### Make generic ID's for sample type ######

# Get ids for filtered fresh stromal 
index <- info2$Cell_Type == "ESC"
ids <- info2[index,]
ids <- ids[order(ids$Study_ID),]

# Make new IDs
ids$New <- 1:nrow(ids)
ids$New2 <- paste("Stromal_",ids$New,sep = "")


# Get ids for filtered fresh epithelial
index2 <- info2$Cell_Type == "EEC"
ids_E <- info2[index2,]
ids_E <- ids_E[order(ids_E$Study_ID),]

# Make new IDs
ids_E$New <- 1:nrow(ids_E)
ids_E$New2 <- paste("Epithelial_",ids_E$New,sep = "")

## merge the stromal and epithelial
info3 <- full_join(ids, ids_E)


# Replace old with new sample ids
index <- which(colnames(fresh) %in% info3$ID)
fresh2 <- fresh[,index]
fresh2 <- fresh2[,order(colnames(fresh2))]
info3 <- info3[order(info3$ID),]
index <- which(info3$ID==colnames(fresh2))
colnames(fresh2) <- info3$New2

################## Export the sample Info #################################
write.csv(info3, "/Volumes/ENDOOMICBM-Q4641/M2_paper/outs/Sample_info1.csv") 


####################################### Exclude problem samples and re-perform ##########################################################

#### Exclude samples that do not cluster appropriately based on PCA

values_to_exclude <- c('Stromal_4','Stromal_3', 'Stromal_9')

info4 <- info3 %>% 
  filter(!(New2 %in% values_to_exclude))

################## Export the sample Info #################################
write.csv(info4, "/Volumes/ENDOOMICBM-Q4641/M2_paper/outs/Sample_info2.csv")

# Exclude columns "Stromal_4", "Stromal_3", "Stromal_9"
fresh3 <- fresh2[, !(colnames(fresh2) %in% c("Stromal_4", "Stromal_3", "Stromal_9"))]

####################################################################################################################################################
############################################################ TMM NORMALISATION ###############################################################
####################################################################################################################################################

y = DGEList(counts=fresh3)
y = calcNormFactors(y, method = "TMM") # estimate the normalisation factor using TMM
# Get normalised expression
Norm <- as.data.frame(cpm(y,normalized.lib.sizes=TRUE, log = T))
nc = cpm(y, normalized.lib.sizes = TRUE) # extract the normalised count (in CPM unit)
range(nc)
range(Norm)


##### Assign gene names
info=read.csv('ESC_bulk/Gene_information_ENSEMBL.csv',
              stringsAsFactors = F)
info=info[,6:7]
nc <- as.data.frame(nc)
nc$ID <- rownames(nc)


############# normalised cell counts data frame ##################
rna_Seq_all3=merge(nc,info, by.x='ID',by.y='geneid')

#########################################################################################################################################
################################### Principle component analysis ##################################################################
##########################################################################################################################################

# Get counts into appropriate format
class(rna_Seq_all3)
sapply(rna_Seq_all3, class)
any(is.na(rna_Seq_all3))

rna_Seq_all3 <- rna_Seq_all3[!duplicated(rna_Seq_all3$gene_name),]
rownames(rna_Seq_all3) <- rna_Seq_all3$gene_name
rna_Seq_all3$gene_name <- NULL
rna_Seq_all3$ID <- NULL

# Calculate Principle Components (PCA Analysis)
pca <- prcomp(t(as.matrix(rna_Seq_all3)))

# Calculate percentage variation
percentage1 <- round(pca$sdev / sum(pca$sdev) * 100, 2)
percentage1 <- paste(colnames(pca$x),"(", paste(as.character(percentage1),"%",sep=""),")", sep="")
pca2 <- as.data.frame(pca$x)
pca2



#######################################################################################################################################
###################################### Differential gene expression analysis ##########################################################
#######################################################################################################################################


################## EEC v ESC comparison ##########################
# COVARIATES FITTING In MODEL/DESIGN

stage = as.factor(info4$Stage_Histology)
dis = as.factor(info4$Endo_Current)
cell = as.factor(info4$Cell_Type)

########################################
##########  limma   ####################
########################################

design = model.matrix(~0 + cell + stage + dis)
head(design)

v <- voom(y, design, plot=F)
fit.v <- lmFit(v, design)


################## EEC v ESC comparison ##########################
con = makeContrasts(cellESC - cellEEC, levels = design)
fit.v2 <- contrasts.fit(fit.v, contrasts=con)
efit <- eBayes(fit.v2)
summary(decideTests(efit))
res.v <- topTable(efit,  adjust.method = "BH", sort.by = "P", n = nrow(fresh2))
head(res.v)

res.v.genes <- data.frame(ENSEMBL = rownames(res.v), res.v)
gene_symbol=read.csv('/Volumes/ENDOOMICBM-Q4641/ESC_bulk/Gene_information_ENSEMBL.csv',
                     stringsAsFactors = F)
gene_symbol=gene_symbol[,6:7]
res.v.genes=merge(res.v.genes,gene_symbol, by.x='ENSEMBL',by.y='geneid')
res.v.genes <- res.v.genes[order(res.v.genes$adj.P.Val), ]

##################################################

write.csv(res.v.genes,"/Volumes/ENDOOMICBM-Q4641/M2_paper/outs/DEG_EEC_v_ESC_all.csv", row.names = F)
res.v.sig = res.v.genes[(res.v.genes$adj.P.Val < 0.05),]
write.csv(res.v.sig,"/Volumes/ENDOOMICBM-Q4641/M2_paper/outs/DEG_EEC_v_ESC_sig.csv", row.names = F)

######## Subset Top table for IL genes

### Selecting out the IL genes

IL_genes <- gene_lists$IL_reactome
IL_genes <- gene_symbol[gene_symbol$gene_name %in% IL_genes, ]
efit_IL <- efit[row.names(efit) %in% IL_genes$geneid, ]

res.v <- topTable(efit_IL,  adjust.method = "BH", sort.by = "P", n = nrow(fresh2))
head(res.v)
summary(decideTests(efit_IL))

res.v.genes <- data.frame(ENSEMBL = rownames(res.v), res.v)
res.v.genes=merge(res.v.genes,gene_symbol, by.x='ENSEMBL',by.y='geneid')
res.v.genes <- res.v.genes[order(res.v.genes$adj.P.Val), ]

write.csv(res.v.genes,"/Volumes/ENDOOMICBM-Q4641/M2_paper/outs/DEG_EEC_v_ESC_IL_all.csv", row.names = F)
res.v.sig = res.v.genes[(res.v.genes$adj.P.Val < 0.05),]
write.csv(res.v.sig,"/Volumes/ENDOOMICBM-Q4641/M2_paper/outs/DEG_EEC_v_ESC_IL_sig.csv", row.names = F)

##################### EPITHELIAL Case v control  ###########################################

#### Subset out Epithelial samples from gene counts
epithelial_columns <- grep("^Epithelial", names(nc), value = TRUE)
epithelial_nc <- nc[,epithelial_columns]

#### Subset out Epithelial samples from sample info
index <- grep("^Epithelial", info4$New2)
info_epithelial <- info4[index, ]  

# COVARIATES FITTING In MODEL/DESIGN
stage = as.factor(info_epithelial$Stage_Histology)
dis = as.factor(info_epithelial$Endo_Current)
#cell = as.factor(info4$Cell_Type)

epithelial_nc <- as.matrix(epithelial_nc)
y = DGEList(counts=epithelial_nc)

##################
####  limma   ####
##################
design = model.matrix(~0 + dis + stage)
head(design)

v <- voom(y, design, plot=F)
fit.v <- lmFit(v, design)

########## Endo v Non-Endo ########
con = makeContrasts(disNo - disEndo, levels = design)
fit.v2 <- contrasts.fit(fit.v, contrasts=con)
efit <- eBayes(fit.v2)
summary(decideTests(efit))
res.v <- topTable(efit,  adjust.method = "BH", sort.by = "P", n = nrow(fresh2))
head(res.v)


res.v.genes <- data.frame(ENSEMBL = rownames(res.v), res.v)
res.v.genes=merge(res.v.genes,gene_symbol, by.x='ENSEMBL',by.y='geneid')
res.v.genes <- res.v.genes[order(res.v.genes$adj.P.Val), ]

write.csv(res.v.genes,"/Volumes/ENDOOMICBM-Q4641/M2_paper/outs/DEG_EEC_case_v_control_all.csv", row.names = F)
res.v.sig = res.v.genes[(res.v.genes$adj.P.Val < 0.05),]
write.csv(res.v.sig,"/Volumes/ENDOOMICBM-Q4641/M2_paper/outs/DEG_EEC_case_v_control_sig.csv", row.names = F)

######## Subset out Top table

# Select out the IL genes
IL_genes <- gene_lists$IL_genes
IL_genes <- gene_symbol[gene_symbol$gene_name %in% IL_genes, ]
efit_IL <- efit[row.names(efit) %in% IL_genes$geneid, ]

res.v <- topTable(efit_IL,  adjust.method = "BH", sort.by = "P", n = nrow(fresh2))
head(res.v)
summary(decideTests(efit_IL))

res.v.genes <- data.frame(ENSEMBL = rownames(res.v), res.v)
res.v.genes=merge(res.v.genes,gene_symbol, by.x='ENSEMBL',by.y='geneid')
res.v.genes <- res.v.genes[order(res.v.genes$adj.P.Val), ]

write.csv(res.v.genes,"/Volumes/ENDOOMICBM-Q4641/M2_paper/outs/DEG_EEC_case_v_control_IL_all.csv", row.names = F)
res.v.sig = res.v.genes[(res.v.genes$adj.P.Val < 0.05),]
write.csv(res.v.sig,"/Volumes/ENDOOMICBM-Q4641/M2_paper/outs/DEG_EEC_case_v_control_IL_sig.csv", row.names = F)


############################################################################################################################################################################################################
##############################################################  Bar plots of Gene Expression comparisons #################################################################
#####################################################################################################################################################################

#### CORRECT FOR COVARIATES - LINEAR MODEL ####

## Following nc = cpm(y, normalized.lib.sizes = TRUE, log = T)

# Below is the function, you can add covariates into line 29 but don't change anything else
norm_S4.fun <- function(x, y) {
  # x: expression data	
  # y: sample information
  
  # corrected expression data
  out1 <- array(0, c(nrow(x), ncol(x)))
  
  # chip and position effects		
  out2 <- array(0, c(nrow(x), 3))
  
  for(i in 1:nrow(x)) {
    i1 <- !is.na(as.numeric(x[i,]))
    i2 <- is.na(as.numeric(x[i,]))
    
    if(length(which(as.numeric(i1)==0)) > ncol(x)*0.8) {
      out1[i,] <- NA
      out2[i,1] <- NA
      out2[i,2] <- NA
      out2[i,3] <- NA
    }
    
    else {	
      
      fit <- summary(lm(as.numeric(x[i,]) ~ as.factor(y$Stage_Histology) + as.factor(y$Endo_Ever)))
      out1[i,i1] <- as.numeric(fit$residuals)
      out1[i,i2] <- NA		
      
      fit2 <- summary(lm(out1[i,]~as.numeric(x[i,])))
      out2[i,1] <- fit2$coefficients[2]
      out2[i,2] <- fit2$r.squared	
      
      out2[i,3] <- fit$r.squared
      
    }		
  }
  out2 <- as.data.frame(out2)
  names(out2) <- c("Mod_fit_est", "Mod_fit_R2", "Correction_r^2")
  out <- list(out1,out2)
  return(out)
}

# nc is the normalised expression matrix and sampleinfo is a dataframe containing the covariates
tmp <- norm_S4.fun(rna_Seq_all3, info4)
corrected_expression <- as.data.frame(round(tmp[[1]],4))
rm(tmp)
colnames(corrected_expression) <- colnames(rna_Seq_all3)
rownames(corrected_expression) <- rownames(rna_Seq_all3)



############################################################# IL genes  #################################################################################
index <- which(rownames(corrected_expression) %in% gene_lists$IL_reactome)
IL_genes <- rna_Seq_all3[index,]
IL_genes <- as.data.frame(t(IL_genes))
rownames(IL_genes) <- info4$New2


# Add cycle phase and disease information alongside gene expression
IL_genes$Cell_Type <- info4$Cell_Type
IL_genes$Endometriosis <- info4$Endo_Current
IL_genes$Cycle_stage <- info4$Stage_Histology

# Create data frame to plot in cycle order
dat <- IL_genes
index <- which(dat$Cell_Type=="EEC")
box_1 <- dat[index,]
index <- which(dat$Cell_Type=="ESC")
box_2 <- dat[index,]
mydf <- rbind(box_1,box_2)
mydf$x2 <- factor(mydf$Cell_Type, levels = unique(mydf$Cell_Type))

####### Changing to a long data frame to enable grid plotting #####
mydf_long6 <- gather(mydf, key = "Gene", value = "Expression", -c(Cell_Type, x2, Endometriosis, Cycle_stage))

# Calculate label positions and values
label_data <- mydf_long6 %>%
  group_by(x2, Gene) %>%
  summarize(LabelPos = median(Expression *1.5), .groups = 'keep') %>%
  ungroup()


############################################# Plot and save Cell type comparison figure

png(file=paste("/Volumes/ENDOOMICBM-Q4641/M2_paper/outs/barplot_IL_EEC_v_ESC.png", sep = ""), width=6000, height=3000, res = 300)
ggplot(mydf_long6, aes(x=factor(x2), y= Expression, fill = Gene)) +
  geom_boxplot(color = "black", fill = "white") +
  scale_x_discrete(name = "Cell_Type") +
  scale_y_continuous(name = "Normalised Expression") +
  theme_classic() +
  theme(text = element_text(size = 20),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        legend.position = "none") +
  geom_point(position = position_dodge(width = 0.75), aes(group = Cell_Type), size = 1) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 7)
dev.off() 

##### New chart
# Replace values in the Cell_Type column
mydf_long6$Cell_Type <- recode(mydf_long6$Cell_Type, 
                               "EEC" = "EnEpi", 
                               "ESC" = "EnS")

###Plot 1
pdf("/Volumes/ENDOOMICBM-Q4641/M2_paper/outs/barplot_IL_EEC_ESC_v2.pdf", width = 24, height = 4)
ggplot(mydf_long6, aes(x = factor(Cell_Type), y = Expression)) +
  geom_boxplot(color = "black", fill = "white") +
  scale_x_discrete(name = "Cell_Type") +
  scale_y_continuous(name = "Normalized Expression") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        legend.position = "none",
        panel.border = element_blank(),          
        strip.background = element_blank(),      
        strip.text.x = element_text(size = 14)) +
  geom_point(position = position_dodge(width = 0.01), aes(group = Endometriosis), size = 0.3, shape = 21, fill = "black") +
  facet_wrap(~ Gene, scales = "free_y", ncol = 13)
dev.off()


############################################### Plot and save endometriosis figure 

# Plot 2
pdf("/Volumes/ENDOOMICBM-Q4641/M2_paper/outs/barplot_IL_Case_Control_v2.pdf", width = 24, height = 4)
ggplot(mydf_long6, aes(x = factor(Cell_Type), y = Expression, fill = Endometriosis)) +
  geom_boxplot() +
  scale_fill_manual(values = c("No" = "blue", "Endo" = "orange")) +  # Set custom colors
  scale_x_discrete(name = "Cell_Type") +
  scale_y_continuous(name = "Normalized Expression") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.title = element_blank(),
        panel.border = element_blank(),          
        strip.background = element_blank(),      
        strip.text.x = element_text(size = 14)) +
  geom_point(position = position_dodge(width = 0.75), aes(group = Endometriosis), size = 0.3) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 13)
dev.off()


png(file=paste("/Volumes/ENDOOMICBM-Q4641/M2_paper/outs/barplot_IL_Case_v_Control.png", sep = ""), width=6000, height=2000, res = 300)
ggplot(mydf_long6, aes(x=factor(x2), y= Expression, fill = Endometriosis)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c("No" = "blue", "Endo" = "orange")) +  # Set custom colors
  scale_x_discrete(name = "Cell_Type") +
  scale_y_continuous(name = "Normalised Expression") +
  theme_classic() + 
  theme(text = element_text(size = 20),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        legend.position = "none") +
  geom_point(position = position_dodge(width = 0.75), aes(group = Endometriosis), size = 0.3) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 7)
dev.off() 

############################################################# Menstrual cycle comparison #######################################################
################################################# Removed cycle stage as a covariate ###############################

#### CORRECT FOR COVARIATES - LINEAR MODEL ####

## Following nc = cpm(y, normalized.lib.sizes = TRUE, log = T)

# Below is the function, you can add covariates into line 29 but don't change anything else
norm_S4.fun <- function(x, y) {
  # x: expression data	
  # y: sample information
  
  # corrected expression data
  out1 <- array(0, c(nrow(x), ncol(x)))
  
  # chip and position effects		
  out2 <- array(0, c(nrow(x), 3))
  
  for(i in 1:nrow(x)) {
    i1 <- !is.na(as.numeric(x[i,]))
    i2 <- is.na(as.numeric(x[i,]))
    
    if(length(which(as.numeric(i1)==0)) > ncol(x)*0.8) {
      out1[i,] <- NA
      out2[i,1] <- NA
      out2[i,2] <- NA
      out2[i,3] <- NA
    }
    
    else {	
      
      fit <- summary(lm(as.numeric(x[i,]) ~ as.factor(y$Endo_Ever)))
      out1[i,i1] <- as.numeric(fit$residuals)
      out1[i,i2] <- NA		
      
      fit2 <- summary(lm(out1[i,]~as.numeric(x[i,])))
      out2[i,1] <- fit2$coefficients[2]
      out2[i,2] <- fit2$r.squared	
      
      out2[i,3] <- fit$r.squared
      
    }		
  }
  out2 <- as.data.frame(out2)
  names(out2) <- c("Mod_fit_est", "Mod_fit_R2", "Correction_r^2")
  out <- list(out1,out2)
  return(out)
}

# nc is the normalised expression matrix and sampleinfo is a dataframe containing the covariates
tmp <- norm_S4.fun(rna_Seq_all3, info4)
corrected_expression <- as.data.frame(round(tmp[[1]],4))
rm(tmp)
colnames(corrected_expression) <- colnames(rna_Seq_all3)
rownames(corrected_expression) <- rownames(rna_Seq_all3)



# Plot and save figure with facets for each gene
ggplot(mydf_long6, aes(x = factor(x2), y = Expression, fill = Cycle_stage)) +
  geom_boxplot() +
  scale_x_discrete(name = "Cell_Type") +
  scale_y_continuous(name = "Normalized Expression") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 20),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  geom_point(position = position_dodge(width = 0.75), aes(group = Cycle_stage), size = 0.3) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 3)

png(file=paste("/Volumes/ENDOOMICBM-Q4641/M2_paper/outs/barplot_IL_Pro_v_Sec.png", sep = ""), width=6000, height=2000, res = 300)
ggplot(mydf_long6, aes(x=factor(x2), y= Expression, fill = Cycle_stage)) + geom_boxplot() + scale_x_discrete(name = "Cell_Type") +
  scale_y_continuous(name = "Normalised Expression") +
  theme(text = element_text(size = 20),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  geom_point(position = position_dodge(width = 0.75), aes(group = Cycle_stage), size = 0.3) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 7)
dev.off() 

