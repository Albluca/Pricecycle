# Load libraries ----------------------------------------------------------------

options(java.parameters = "-Xmx4g")

library(dplyr)
library(ggplot2)
library(rpgraph)
library(readr)
library(readxl)
library(GO.db)
library(biomaRt)
library(Pricecycle)

# Set geneset ----------------------------------------------------------------


FreemanData <- read_delim("~/Google Drive/Datasets/Gene List/Freeman.gmt",
                          "\t", escape_double = FALSE, col_names = FALSE,
                          trim_ws = TRUE)

Freeman_CC1 <- unlist(FreemanData[1,-c(1,2)], use.names = FALSE)
Freeman_CC2 <- unlist(FreemanData[2,-c(1,2)], use.names = FALSE)
Freeman_G1S_CC4 <- unlist(FreemanData[3,-c(1,2)], use.names = FALSE)
Freeman_G2M_CC6 <- unlist(FreemanData[4,-c(1,2)], use.names = FALSE)
Freeman_CC9 <- unlist(FreemanData[5,-c(1,2)], use.names = FALSE)
Freeman_G1S_CC4A <- unlist(FreemanData[6,-c(1,2)], use.names = FALSE)
Freeman_S_CC4B <- unlist(FreemanData[7,-c(1,2)], use.names = FALSE)
Freeman_G2_CC6A <- unlist(FreemanData[8,-c(1,2)], use.names = FALSE)
Freeman_M_CC6B <- unlist(FreemanData[9,-c(1,2)], use.names = FALSE)

Freeman_CC1_Mouse <- ConvertNames("human", "mouse", Freeman_CC1)
Freeman_CC2_Mouse <- ConvertNames("human", "mouse", Freeman_CC2)
Freeman_G1S_CC4_Mouse <- ConvertNames("human", "mouse", Freeman_G1S_CC4)
Freeman_G2M_CC6_Mouse <- ConvertNames("human", "mouse", Freeman_G2M_CC6)
Freeman_CC9_Mouse <- ConvertNames("human", "mouse", Freeman_CC9)
Freeman_G1S_CC4A_Mouse <- ConvertNames("human", "mouse", Freeman_G1S_CC4A)
Freeman_S_CC4B_Mouse <- ConvertNames("human", "mouse", Freeman_S_CC4B)
Freeman_G2_CC6A_Mouse <- ConvertNames("human", "mouse", Freeman_G2_CC6A)
Freeman_M_CC6B_Mouse <- ConvertNames("human", "mouse", Freeman_M_CC6B)

AllFreman_Mouse <- c(Freeman_CC1_Mouse, Freeman_CC2_Mouse, Freeman_CC9_Mouse,
                     Freeman_G1S_CC4_Mouse, Freeman_G1S_CC4A_Mouse, Freeman_G2_CC6A_Mouse,
                     Freeman_G2M_CC6_Mouse, Freeman_M_CC6B_Mouse, Freeman_S_CC4B_Mouse)

AllFreman_Mouse_CCP <- c(Freeman_G1S_CC4_Mouse, Freeman_G1S_CC4A_Mouse, Freeman_G2_CC6A_Mouse,
                         Freeman_G2M_CC6_Mouse, Freeman_M_CC6B_Mouse, Freeman_S_CC4B_Mouse)


AllTerms <- GOBPOFFSPRING[["GO:0007049"]]
AllTerms <- sort(c("GO:0007049", AllTerms))

MouseGenes_GOCellCycle <- getBM(attributes = c("external_gene_name"),
                                filters = "go",
                                values = AllTerms,
                                mart = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl"))

MouseGenes_GOCellCycle <- unlist(MouseGenes_GOCellCycle, use.names = FALSE)

# AllWit_Mouse <- ConvertNames("human", "mouse", unlist(StageAssociation_Whit_Ext[-c(1,2)], use.names = FALSE))



# Buettner et al ----------------------------------------------------------------

# Load data

BaseDir <- "~/Google Drive/Datasets/Buettner et al. - Murie embrionic stem cells/"

Murine_ESC_G1_Data <- read_delim(paste(BaseDir, "G1_singlecells_counts.txt", sep= ''), "\t", escape_double = FALSE, trim_ws = TRUE)
Murine_ESC_S_Data <- read_delim(paste(BaseDir, "S_singlecells_counts.txt", sep= ''), "\t", escape_double = FALSE, trim_ws = TRUE)
Murine_ESC_G2M_Data <- read_delim(paste(BaseDir, "G2M_singlecells_counts.txt", sep= ''), "\t", escape_double = FALSE, trim_ws = TRUE)

GeneToConsider <- which(!is.na(Murine_ESC_G1_Data$AssociatedGeneName))

Murine_ESC_All <- cbind(Murine_ESC_G1_Data[GeneToConsider, ],
                        Murine_ESC_S_Data[GeneToConsider, -c(1:4)],
                        Murine_ESC_G2M_Data[GeneToConsider, -c(1:4)])

GNames <- Murine_ESC_All$AssociatedGeneName

GNames[duplicated(GNames)] <- paste(GNames[duplicated(GNames)], 1:sum(duplicated(GNames)), sep ="_R")

Murine_ESC_All <- Murine_ESC_All[,-c(1:4)]
rownames(Murine_ESC_All) <- GNames

Murine_ESC_Stages <- c(rep(x = "G1", ncol(Murine_ESC_G1_Data)-4), rep(x = "S", ncol(Murine_ESC_S_Data)-4), rep(x = "G2M", ncol(Murine_ESC_G2M_Data)-4))
names(Murine_ESC_Stages) <- colnames(Murine_ESC_All)


Murine_ESC_Stages <- factor(Murine_ESC_Stages, levels = c("G1", "S", "G2M"))

Factors <- scran::computeSumFactors(x = data.matrix(Murine_ESC_All), positive = TRUE)
Normalized_Exp <- t(t(data.matrix(Murine_ESC_All))/Factors)


GeneExprMat <- Murine_ESC_All
StartSet <- MouseGenes_GOCellCycle
Categories <- Murine_ESC_Stages



Output.Buet <- SelectGenesOnGraph(DataSet = GeneExprMat,
                                  StartSet = MouseGenes_GOCellCycle,
                                  Categories = Categories,
                                  PCACenter = TRUE, PCAProjCenter = TRUE)


# BuettData <- read_rds("~/Google Drive/Datasets/Buettner et al. - Murie embrionic stem cells/rPG_Last.rds")

write_rds(x = list(Analysis = Output.Buet,
                   ExpMat = log10(GeneExprMat+1),
                   Cats = Categories),
          path = "~/Google Drive/Datasets/Buettner et al. - Murie embrionic stem cells/rPG_Last.rds")




GeneExprMat <- Normalized_Exp[, Factors>0]
StartSet <- MouseGenes_GOCellCycle
Categories <- Murine_ESC_Stages[Factors>0]


Output.Buet <- SelectGenesOnGraph(DataSet = GeneExprMat,
                                  StartSet = MouseGenes_GOCellCycle,
                                  Categories = Categories,
                                  PCACenter = TRUE, PCAProjCenter = TRUE)

write_rds(x = list(Analysis = Output.Buet,
                   ExpMat = log10(GeneExprMat+1),
                   Cats = Categories),
          path = "~/Google Drive/Datasets/Buettner et al. - Murie embrionic stem cells/rPG_Last-Filtered.rds")








# Kowalczyk et al ----------------------------------------------------------------

# Load data

BaseDir <- "~/Google Drive/Datasets/Kowalczyk et al - Murine hematopoietic stem cells/"

Murine_HSC_Data <- read_excel(paste(BaseDir, "GSE59114_C57BL6_GEO_all.xlsx", sep = ''))
Murine_HSC_Samples <- read_excel(paste(BaseDir, "Table_S2.xlsx", sep = ''))

CCVect <- factor(Murine_HSC_Samples$`Estimated phase`, levels = c("G0", "G1(early)", "G1(late)", "S", "G2/M"))
names(CCVect) <- Murine_HSC_Samples$`cell ID`

ProgVect <- Murine_HSC_Samples$`Progression Rank`
names(ProgVect) <- Murine_HSC_Samples$`cell ID`

SelectedSamples <- Murine_HSC_Samples$`cell ID`[Murine_HSC_Samples$population %in% c("'old LT-HSC'", "'old LT-HSC(rep)'", "'old ST-HSC'",
                                                                                     "'old ST-HSC(rep)'", "'young LT-HSC'", "'young ST-HSC'",
                                                                                     "'young MPP'", "'old MPP'")]

SelectedSamples <- intersect(SelectedSamples, colnames(Murine_HSC_Data))

Org_Categories <- CCVect[paste(SelectedSamples)]

AllData <- data.matrix(Murine_HSC_Data[-1,paste(SelectedSamples)])
rownames(AllData) <- gsub("'", '', unlist(Murine_HSC_Data[-1, 1]))

Factors <- scran::computeSumFactors(x = data.matrix(2^AllData - 1), positive = TRUE)
Normalized_Exp <- t(t(data.matrix(2^AllData - 1))/Factors)


GeneExprMat <- data.matrix(2^AllData - 1)
StartSet <- MouseGenes_GOCellCycle
Categories <- Org_Categories


Output.Kowa <- SelectGenesOnGraph(DataSet = GeneExprMat, StartSet = MouseGenes_GOCellCycle, Categories = Categories,
                                  PCACenter = TRUE, PCAProjCenter = TRUE)


write_rds(x = list(Analysis = Output.Kowa,
                   ExpMat = log10(GeneExprMat + 1),
                   Cats = Categories),
          path = "~/Google Drive/Datasets/Kowalczyk et al - Murine hematopoietic stem cells/rPG_Last.rds")







GeneExprMat <- Normalized_Exp[, Factors>0]
StartSet <- MouseGenes_GOCellCycle
Categories <- Org_Categories[Factors>0]


# GeneExprMat <- data.matrix(AllData)
# StartSet <- MouseGenes_GOCellCycle
# Categories <- Categories


Output.Kowa.Norm <- SelectGenesOnGraph(DataSet = GeneExprMat, StartSet = MouseGenes_GOCellCycle, Categories = Categories,
                                  PCACenter = TRUE, PCAProjCenter = TRUE)


write_rds(x = list(Analysis = Output.Kowa.Norm,
                   ExpMat = log10(GeneExprMat + 1),
                   Cats = Categories),
          path = "~/Google Drive/Datasets/Kowalczyk et al - Murine hematopoietic stem cells/rPG_Last-Filtered.rds")




# Kowalczyk et al - G0 (PriceCycle) ----------------------------------------------------------------

Data.Kowa <- read_rds("~/Google Drive/Datasets/Kowalczyk et al - Murine hematopoietic stem cells/rPG_Last.rds")
Data.Kowa.Norm <- read_rds("~/Google Drive/Datasets/Kowalczyk et al - Murine hematopoietic stem cells/rPG_Last-Filtered.rds")

# Plot data and construct auxiliary structures

TargetStruct <- Data.Kowa$Analysis$PGStructs[[length(Data.Kowa$Analysis$PGStructs)]]
Proc.Exp.Kowa <- PlotOnStages(Structure = "Circle",
                              Categories = TargetStruct$Categories,
                              nGenes = 2,
                              TaxonList = TargetStruct$TaxonList[[length(TargetStruct$TaxonList)]],
                              PrinGraph = TargetStruct$PrinGraph,
                              Net = TargetStruct$Net[[length(TargetStruct$Net)]],
                              SelThr = .29,
                              ComputeOverlaps = TRUE,
                              ExpData = TargetStruct$FiltExp,
                              RotatioMatrix = TargetStruct$PCAData$rotation[,1:TargetStruct$nDims],
                              PCACenter = TargetStruct$PCAData$center,
                              PointProjections = TargetStruct$ProjPoints[[length(TargetStruct$ProjPoints)]],
                              OrderOnCat = TRUE,
                              SmoothPoints = 1,
                              MinCellPerNode = 2,
                              Title = 'Kowalczyk et al')



TargetStruct <- Data.Kowa.Norm$Analysis$PGStructs[[length(Data.Kowa.Norm$Analysis$PGStructs)]]
Proc.Exp.Kowa.Norm <- PlotOnStages(Structure = "Circle",
                                   Categories = TargetStruct$Categories,
                                   nGenes = 2,
                                   TaxonList = TargetStruct$TaxonList[[length(TargetStruct$TaxonList)]],
                                   PrinGraph = TargetStruct$PrinGraph,
                                   Net = TargetStruct$Net[[length(TargetStruct$Net)]],
                                   SelThr = .33,
                                   ComputeOverlaps = TRUE,
                                   ExpData = TargetStruct$FiltExp,
                                   RotatioMatrix = TargetStruct$PCAData$rotation[,1:TargetStruct$nDims],
                                   PCACenter = TargetStruct$PCAData$center,
                                   PointProjections = TargetStruct$ProjPoints[[length(TargetStruct$ProjPoints)]],
                                   OrderOnCat = TRUE,
                                   SmoothPoints = 1,
                                   MinCellPerNode = 2,
                                   Title = 'Kowalczyk et al')


SelStageInfo.Strong.Norm <- read_rds("~/Desktop/SelStageInfo.Strong.Norm.rds")
SelStageInfo.Strong <- read_rds("~/Desktop/SelStageInfo.Strong.rds")


Staged.Proc.Exp.Kowa <- StageWithPeaks(DataStruct = Data.Kowa,
                                         ProcStruct = Proc.Exp.Kowa,
                                         ComputeG0 = TRUE,
                                         FiltMax = 0,
                                         Thr = .8,
                                         QuantSel = .8,
                                         StageInfo = SelStageInfo.Strong.Norm,
                                         MinNodes = 2,
                                         Mode = 1,
                                         G0Level = .9,
                                         Title = 'Kowalczyk et al')


Staged.Proc.Exp.Kowa.Norm <- StageWithPeaks(DataStruct = Data.Kowa.Norm,
                                       ProcStruct = Proc.Exp.Kowa.Norm,
                                       ComputeG0 = TRUE,
                                       FiltMax = 0,
                                       Thr = .8,
                                       QuantSel = .8,
                                       StageInfo = SelStageInfo.Strong.Norm,
                                       MinNodes = 2,
                                       Mode = 1,
                                       G0Level = .9,
                                       Title = 'Kowalczyk et al')


Inf_NonG0 <- names(Staged.Proc.Exp.Kowa$CellStages_Ext[Staged.Proc.Exp.Kowa$CellStages != "G0"])
Inf_NonG0.Norm <- names(Staged.Proc.Exp.Kowa.Norm$CellStages_Ext[Staged.Proc.Exp.Kowa.Norm$CellStages != "G0"])

ToKeep <- union(Inf_NonG0, Inf_NonG0.Norm)
ToKeep.2 <- intersect(Inf_NonG0, Inf_NonG0.Norm)


# Load data

BaseDir <- "~/Google Drive/Datasets/Kowalczyk et al - Murine hematopoietic stem cells/"

Murine_HSC_Data <- read_excel(paste(BaseDir, "GSE59114_C57BL6_GEO_all.xlsx", sep = ''))
Murine_HSC_Samples <- read_excel(paste(BaseDir, "Table_S2.xlsx", sep = ''))

CCVect <- factor(Murine_HSC_Samples$`Estimated phase`, levels = c("G0", "G1(early)", "G1(late)", "S", "G2/M"))
names(CCVect) <- Murine_HSC_Samples$`cell ID`

ProgVect <- Murine_HSC_Samples$`Progression Rank`
names(ProgVect) <- Murine_HSC_Samples$`cell ID`

SelectedSamples <- Murine_HSC_Samples$`cell ID`[Murine_HSC_Samples$population %in% c("'old LT-HSC'", "'old LT-HSC(rep)'", "'old ST-HSC'",
                                                                                     "'old ST-HSC(rep)'", "'young LT-HSC'", "'young ST-HSC'",
                                                                                     "'young MPP'", "'old MPP'")]

SelectedSamples <- intersect(SelectedSamples, colnames(Murine_HSC_Data))

Org_Categories <- CCVect[paste(SelectedSamples)]

AllData <- data.matrix(Murine_HSC_Data[-1,paste(SelectedSamples)])
rownames(AllData) <- gsub("'", '', unlist(Murine_HSC_Data[-1, 1]))



AllData.1 <- AllData[,ToKeep]
Org_Categories.1 <- Org_Categories[ToKeep]

AllData.2 <- AllData[,ToKeep.2]
Org_Categories.2 <- Org_Categories[ToKeep.2]

# Org_Cat_Names <- names(Org_Categories)
# Org_Categories <- factor(as.character(Org_Categories), levels = levels(Org_Categories)[-1])
# names(Org_Categories) <- Org_Cat_Names

Factors.1 <- scran::computeSumFactors(x = data.matrix(2^AllData.1 - 1), positive = TRUE)
Normalized_Exp.1 <- t(t(data.matrix(2^AllData.1 - 1))/Factors.1)

Factors.2 <- scran::computeSumFactors(x = data.matrix(2^AllData.2 - 1), positive = TRUE)
Normalized_Exp.2 <- t(t(data.matrix(2^AllData.2 - 1))/Factors.2)


GeneExprMat <- data.matrix(2^AllData.1 - 1)
StartSet <- MouseGenes_GOCellCycle
Categories <- Org_Categories.1


Output.Kowa.NoG0 <- SelectGenesOnGraph(DataSet = GeneExprMat, StartSet = MouseGenes_GOCellCycle, Categories = Categories,
                                       PCACenter = TRUE, PCAProjCenter = TRUE)


write_rds(x = list(Analysis = Output.Kowa.NoG0,
                   ExpMat = log10(GeneExprMat + 1),
                   Cats = Categories),
          path = "~/Google Drive/Datasets/Kowalczyk et al - Murine hematopoietic stem cells/rPG_Last.NoG0.PC.rds")


GeneExprMat <- Normalized_Exp.1[, Factors.1>0]
StartSet <- MouseGenes_GOCellCycle
Categories <- Org_Categories.1[Factors.1>0]


# GeneExprMat <- data.matrix(AllData)
# StartSet <- MouseGenes_GOCellCycle
# Categories <- Categories


Output.Kowa.NoG0 <- SelectGenesOnGraph(DataSet = GeneExprMat, StartSet = MouseGenes_GOCellCycle, Categories = Categories,
                                       PCACenter = TRUE, PCAProjCenter = TRUE)


write_rds(x = list(Analysis = Output.Kowa.NoG0,
                   ExpMat = log10(GeneExprMat + 1),
                   Cats = Categories),
          path = "~/Google Drive/Datasets/Kowalczyk et al - Murine hematopoietic stem cells/rPG_Last.NoG0.PC-Filtered.rds")














GeneExprMat <- data.matrix(2^AllData.2 - 1)
StartSet <- MouseGenes_GOCellCycle
Categories <- Org_Categories.2


Output.Kowa.NoG0 <- SelectGenesOnGraph(DataSet = GeneExprMat, StartSet = MouseGenes_GOCellCycle, Categories = Categories,
                                       PCACenter = TRUE, PCAProjCenter = TRUE)


write_rds(x = list(Analysis = Output.Kowa.NoG0,
                   ExpMat = log10(GeneExprMat + 1),
                   Cats = Categories),
          path = "~/Google Drive/Datasets/Kowalczyk et al - Murine hematopoietic stem cells/rPG_Last.NoG0.PC.2.rds")


GeneExprMat <- Normalized_Exp.2[, Factors.2>0]
StartSet <- MouseGenes_GOCellCycle
Categories <- Org_Categories.2[Factors.2>0]


# GeneExprMat <- data.matrix(AllData)
# StartSet <- MouseGenes_GOCellCycle
# Categories <- Categories


Output.Kowa.NoG0 <- SelectGenesOnGraph(DataSet = GeneExprMat, StartSet = MouseGenes_GOCellCycle, Categories = Categories,
                                       PCACenter = TRUE, PCAProjCenter = TRUE)


write_rds(x = list(Analysis = Output.Kowa.NoG0,
                   ExpMat = log10(GeneExprMat + 1),
                   Cats = Categories),
          path = "~/Google Drive/Datasets/Kowalczyk et al - Murine hematopoietic stem cells/rPG_Last.NoG0.PC.2-Filtered.rds")























# Sasagawa et al ----------------------------------------------------------------

BaseDir <- "~/Google Drive/Datasets/Sasagawa et al - Murine Stem Cells (EB5 cell line)/"

FilesToRead <- list.files(path = paste(BaseDir, "ES/", sep = ''), full.names = TRUE, pattern = "txt.gz")
LoadedData <- list()

FilesToReadShort <- list.files(path = paste(BaseDir, "ES/", sep = ''), full.names = FALSE, pattern = "txt.gz")

ID <- unlist(lapply(strsplit(FilesToRead, "_"), "[[", 2))
StageVect <- rep("S_G1", length(ID))
StageVect[ID == "ESS"] <- "S_S"
StageVect[ID == "ESM"] <- "S_G2/M"

StageVect <- factor(StageVect, levels = c("S_G1", "S_S", "S_G2/M"))

for(i in 1:length(FilesToRead)){
  Data <- read_delim(FilesToRead[i], "\t", escape_double = FALSE, trim_ws = TRUE)
  LoadedData[[i]] <- list(Name = FilesToRead[i], Data = Data, Stage = StageVect[i])
}

AllGenes <- NULL

for(i in 1:length(LoadedData)){
  AllGenes <- unique(c(AllGenes, LoadedData[[i]]$Data$gene.symbol))
}

AllGenesEns <- rep(NA, length(AllGenes))
names(AllGenesEns) <- AllGenes

for(i in 1:length(LoadedData)){
  AllGenesEns[LoadedData[[i]]$Data$gene.symbol] <- LoadedData[[i]]$Data$id
}

ValMat <- matrix(rep(0, length(AllGenes)*length(FilesToRead)), ncol = length(AllGenes))

colnames(ValMat) <- AllGenes
rownames(ValMat) <- FilesToRead

for(i in 1:length(LoadedData)){
  ValMat[i,LoadedData[[i]]$Data$gene.symbol] <- LoadedData[[i]]$Data$fpkm
}

Factors <- scran::computeSumFactors(x = t(data.matrix(ValMat)), sizes=c(5, 10, 15), positive = TRUE)
Normalized_Exp <- t(data.matrix(ValMat)/Factors)


# GeneExprMat <- (2^AllData - 1)[, Factors>0]
# StartSet <- MouseGenes_GOCellCycle
# Categories <- Categories[Factors>0]
#

GeneExprMat <- t(data.matrix(ValMat))
StartSet <- MouseGenes_GOCellCycle
Categories <- factor(StageVect)


Output.Sasa <- SelectGenesOnGraph(DataSet = GeneExprMat, StartSet = MouseGenes_GOCellCycle, Categories = Categories,
                                  PCACenter = TRUE, PCAProjCenter = TRUE)


write_rds(x = list(Analysis = Output.Sasa,
                   ExpMat = log10(GeneExprMat + 1),
                   Cats = Categories),
          path = "~/Google Drive/Datasets/Sasagawa et al - Murine Stem Cells (EB5 cell line)/rPG_Last.rds")






GeneExprMat <- Normalized_Exp[, Factors>0]
StartSet <- MouseGenes_GOCellCycle
Categories <- factor(StageVect)[Factors>0]


Output.Sasa <- SelectGenesOnGraph(DataSet = GeneExprMat, StartSet = MouseGenes_GOCellCycle, Categories = Categories,
                                  PCACenter = TRUE, PCAProjCenter = TRUE)


write_rds(x = list(Analysis = Output.Sasa,
                   ExpMat = log10(GeneExprMat + 1),
                   Cats = Categories),
          path = "~/Google Drive/Datasets/Sasagawa et al - Murine Stem Cells (EB5 cell line)/rPG_Last-Filtered.rds")






















# Kowalczyk et al - G0 (comp) removed ----------------------------------------------------------------

# Load data

BaseDir <- "~/Google Drive/Datasets/Kowalczyk et al - Murine hematopoietic stem cells/"

Murine_HSC_Data <- read_excel(paste(BaseDir, "GSE59114_C57BL6_GEO_all.xlsx", sep = ''))
Murine_HSC_Samples <- read_excel(paste(BaseDir, "Table_S2.xlsx", sep = ''))

CCVect <- factor(Murine_HSC_Samples$`Estimated phase`, levels = c("G0", "G1(early)", "G1(late)", "S", "G2/M"))
names(CCVect) <- Murine_HSC_Samples$`cell ID`

ProgVect <- Murine_HSC_Samples$`Progression Rank`
names(ProgVect) <- Murine_HSC_Samples$`cell ID`

SelectedSamples <- Murine_HSC_Samples$`cell ID`[Murine_HSC_Samples$population %in% c("'old LT-HSC'", "'old LT-HSC(rep)'", "'old ST-HSC'",
                                                                                     "'old ST-HSC(rep)'", "'young LT-HSC'", "'young ST-HSC'",
                                                                                     "'young MPP'", "'old MPP'")]

SelectedSamples <- intersect(SelectedSamples, colnames(Murine_HSC_Data))

Org_Categories <- CCVect[paste(SelectedSamples)]

AllData <- data.matrix(Murine_HSC_Data[-1,paste(SelectedSamples)])
rownames(AllData) <- gsub("'", '', unlist(Murine_HSC_Data[-1, 1]))


Prolif <- names(Org_Categories[Org_Categories != "G0"])


AllData <- AllData[,Prolif]
Org_Categories <- Org_Categories[Prolif]
Org_Cat_Names <- names(Org_Categories)
Org_Categories <- factor(as.character(Org_Categories), levels = levels(Org_Categories)[-1])
names(Org_Categories) <- Org_Cat_Names


Factors <- scran::computeSumFactors(x = data.matrix(2^AllData - 1), positive = TRUE)
Normalized_Exp <- t(t(data.matrix(2^AllData - 1))/Factors)


GeneExprMat <- data.matrix(2^AllData - 1)
StartSet <- MouseGenes_GOCellCycle
Categories <- Org_Categories


Output.Kowa.NoG0 <- SelectGenesOnGraph(DataSet = GeneExprMat, StartSet = MouseGenes_GOCellCycle, Categories = Categories,
                                  PCACenter = TRUE, PCAProjCenter = TRUE)


write_rds(x = list(Analysis = Output.Kowa.NoG0,
                   ExpMat = log10(GeneExprMat + 1),
                   Cats = Categories),
          path = "~/Google Drive/Datasets/Kowalczyk et al - Murine hematopoietic stem cells/rPG_Last.NoG0.rds")


GeneExprMat <- Normalized_Exp[, Factors>0]
StartSet <- MouseGenes_GOCellCycle
Categories <- Org_Categories[Factors>0]


# GeneExprMat <- data.matrix(AllData)
# StartSet <- MouseGenes_GOCellCycle
# Categories <- Categories


Output.Kowa.NoG0 <- SelectGenesOnGraph(DataSet = GeneExprMat, StartSet = MouseGenes_GOCellCycle, Categories = Categories,
                                  PCACenter = TRUE, PCAProjCenter = TRUE)


write_rds(x = list(Analysis = Output.Kowa.NoG0,
                   ExpMat = log10(GeneExprMat + 1),
                   Cats = Categories),
          path = "~/Google Drive/Datasets/Kowalczyk et al - Murine hematopoietic stem cells/rPG_Last.NoG0-Filtered.rds")








