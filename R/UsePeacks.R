#' Title
#'
#' @param DataStruct
#' @param ProcStruct
#' @param MinNodes
#' @param FiltMax
#' @param Thr
#' @param QuantSel
#' @param SinglePeack
#'
#' @return
#' @export
#'
#' @examples
GetGenesWithPeaks <- function(DataStruct,
                              ProcStruct,
                              AllGenes = FALSE,
                              MinNodes = 2,
                              FiltMax = .2,
                              Thr = .9,
                              QuantSel = .75,
                              MedianCVThr = 1,
                              StagesNames = c("G1", "S", "G2")) {

  WorkPG <- DataStruct$Analysis$PGStructs[[length(DataStruct$Analysis$PGStructs)]]
  TaxList <- WorkPG$TaxonList[[length(WorkPG$TaxonList)]]
  TaxVect <- rep(NA, nrow(WorkPG$Data))
  names(TaxVect) <- rownames(WorkPG$Data)

  lapply(1:length(TaxList), function(j){
    TaxVect[ TaxList[[j]] ] <<- j
  })


  if(AllGenes){

    GeneExprMat.DF <- data.frame(t(DataStruct$ExpMat[rowSums(DataStruct$ExpMat)>0, names(ProcStruct$CellsPT)]))
    colnames(GeneExprMat.DF) <- rownames(DataStruct$ExpMat)[rowSums(DataStruct$ExpMat)>0]

    SelCell <- names(TaxVect)

    GeneExprMat.DF.Split <- split(GeneExprMat.DF[SelCell, ], TaxVect[SelCell])
    GeneExprMat.DF.Split.Mean <- lapply(GeneExprMat.DF.Split, colMeans)
    GeneExprMat.DF.Split.Sd <- lapply(GeneExprMat.DF.Split, function(x){apply(x, 2, sd)})

    GeneExprMat.DF.Split.Mean.Bind <- sapply(GeneExprMat.DF.Split.Mean, cbind)
    GeneExprMat.DF.Split.Sd.Bind <- sapply(GeneExprMat.DF.Split.Sd, cbind)

    rownames(GeneExprMat.DF.Split.Mean.Bind) <- colnames(GeneExprMat.DF)
    colnames(GeneExprMat.DF.Split.Mean.Bind) <- names(GeneExprMat.DF.Split.Mean)

    rownames(GeneExprMat.DF.Split.Sd.Bind) <- colnames(GeneExprMat.DF)
    colnames(GeneExprMat.DF.Split.Sd.Bind) <- names(GeneExprMat.DF.Split.Sd)

    Reord <- ProcStruct$ExtPath[-length(ProcStruct$ExtPath)]
    Reord <- Reord[Reord %in% colnames(GeneExprMat.DF.Split.Mean.Bind)]

    GeneExprMat.DF.Split.Mean.Bind <- GeneExprMat.DF.Split.Mean.Bind[,Reord]
    GeneExprMat.DF.Split.Sd.Bind <- GeneExprMat.DF.Split.Sd.Bind[,Reord]

    MedianCV <- apply(GeneExprMat.DF.Split.Sd.Bind/GeneExprMat.DF.Split.Mean.Bind, 1, median, na.rm=TRUE)

    # NormExp <- ProcStruct$NodesExp

    NormExp <- GeneExprMat.DF.Split.Mean.Bind
    # dim(NormExp)

  } else {

    NormExp <- ProcStruct$NodesExp
    colnames(NormExp) <- ProcStruct$ExtPath

    MedianCV <- rep(0, nrow(ProcStruct$NodesExp))

  }


  MedianCV <- MedianCV[rowSums(NormExp>0) > MinNodes]
  NormExp <- NormExp[rowSums(NormExp>0) > MinNodes, ]

  # dim(NormExp)

  SDVect <- apply(NormExp, 1, sd)

  SDVect[is.na(SDVect)] <- 0
  MedianCV[is.na(MedianCV)] <- 0

  NormExp <- NormExp[SDVect > quantile(SDVect, QuantSel) & MedianCV < MedianCVThr,]

  NormExp[NormExp <= FiltMax] <- 0
  NormExp <- NormExp[rowSums(NormExp>0) > MinNodes, ]

  NormExp <- (NormExp - apply(NormExp, 1, min))/(apply(NormExp, 1, max) - apply(NormExp, 1, min))

  pheatmap::pheatmap(NormExp, cluster_cols = FALSE)

  FiltStagesOnNodes <- list()

  for(stg in StagesNames){
    Selected <- ProcStruct$StageOnNodes[grep(stg, names(ProcStruct$StageOnNodes))]
    FiltStagesOnNodes[[stg]] <- unlist(Selected, use.names = FALSE)
    names(FiltStagesOnNodes[[stg]]) <- unlist(lapply(Selected, names), use.names = FALSE)
  }


  StageGenes <- lapply(FiltStagesOnNodes, function(x){

    SelX <- intersect(names(x), colnames(NormExp))
    UnSelX <- setdiff(colnames(NormExp), SelX)

    if(length(SelX)>0){
      if(length(SelX)==1){
        NormExp[, UnSelX] < Thr & NormExp[, SelX] > Thr
      } else {
        apply(NormExp[, UnSelX] < Thr, 1, all) & apply(NormExp[, SelX] > Thr, 1, any)
      }
    } else {
      NA
    }

  })

  StageGenes.Names.UP <- lapply(StageGenes, function(x){
    names(which(x))
  })


  tt <- lapply(as.list(1:length(StageGenes.Names.UP)), function(i) {

    if(length(StageGenes.Names.UP[[i]])>1){
      pheatmap::pheatmap(NormExp[StageGenes.Names.UP[[i]],], cluster_cols = FALSE, main = names(StageGenes.Names.UP)[i])
    }

    if(length(StageGenes.Names.UP[[i]])==1){
      pheatmap::pheatmap(NormExp[StageGenes.Names.UP[[i]],], cluster_cols = FALSE, cluster_rows = FALSE, main = names(StageGenes.Names.UP)[i])
    }

  })

  barplot(unlist(lapply(StageGenes.Names.UP, length)), las = 2, horiz = FALSE, main = "Up genes")













  StageGenes <- lapply(FiltStagesOnNodes, function(x){

    SelX <- intersect(names(x), colnames(NormExp))
    UnSelX <- setdiff(colnames(NormExp), SelX)

    if(length(SelX)>0){
      if(length(SelX)==1){
        NormExp[, UnSelX] > 1 - Thr & NormExp[, SelX] < 1 - Thr
      } else {
        apply(NormExp[, UnSelX] > 1- Thr, 1, all) & apply(NormExp[, SelX] < 1 - Thr, 1, any)
      }
    } else {
      NA
    }

  })

  StageGenes.Names.DOWN <- lapply(StageGenes, function(x){
    names(which(x))
  })


  tt <- lapply(as.list(1:length(StageGenes.Names.DOWN)), function(i) {

    if(length(StageGenes.Names.DOWN[[i]])>1){
      pheatmap::pheatmap(NormExp[StageGenes.Names.DOWN[[i]],], cluster_cols = FALSE, main = names(StageGenes.Names.DOWN)[i])
    }

    if(length(StageGenes.Names.DOWN[[i]])==1){
      pheatmap::pheatmap(NormExp[StageGenes.Names.DOWN[[i]],], cluster_cols = FALSE, cluster_rows = FALSE, main = names(StageGenes.Names.DOWN)[i])
    }

  })

  barplot(unlist(lapply(StageGenes.Names.DOWN, length)), las = 2, horiz = FALSE, main = "DOWN genes")









  return(list(UP = StageGenes.Names.UP, DOWN = StageGenes.Names.DOWN))

}
















#' Title
#'
#' @param StageMatrix
#' @param NodePenalty
#'
#' @return
#' @export
#'
#' @examples
FitStagesCirc <- function(StageMatrix, NodePenalty, Mode = 1) {

  NormStageMatrix <- StageMatrix

  # Find which columns are associated with at least one stage

  ToAnalyze <- which(colSums(StageMatrix)>0)

  if(length(ToAnalyze)<nrow(StageMatrix)){

    # Not enough columns. Padding around them

    PaddedIdxs <- c(
      (min(ToAnalyze) - nrow(StageMatrix)):(min(ToAnalyze)-1),
      ToAnalyze,
      (max(ToAnalyze)+1):(max(ToAnalyze) + nrow(StageMatrix))
    )

    PaddedIdxs[PaddedIdxs > ncol(StageMatrix)] <-  PaddedIdxs[PaddedIdxs > ncol(StageMatrix)] - ncol(StageMatrix)
    PaddedIdxs[PaddedIdxs <= 0] <-  PaddedIdxs[PaddedIdxs <= 0] + ncol(StageMatrix)

    ToAnalyze <- sort(unique(PaddedIdxs))

  }

  # Selecting the matrix that will be analysed and saving the indices

  NormStageMatrix <- NormStageMatrix[,ToAnalyze]
  NormStageMatrixIdx <- ToAnalyze

  NormNodePenalty <- NodePenalty[ToAnalyze]

  Possibilities <- combn(c(1:ncol(NormStageMatrix), rep(NA, nrow(NormStageMatrix))), nrow(NormStageMatrix))

  NoChange <- apply(is.na(Possibilities), 2, sum) == nrow(NormStageMatrix)

  ToKeep <- (!apply(Possibilities, 2, is.unsorted, na.rm = TRUE))
  #   (apply(!is.na(Possibilities), 2, sum) != 1) &
  #   (!is.na(Possibilities[nrow(NormStageMatrix),]))

  Possibilities <- Possibilities[, ToKeep | NoChange]
  dim(Possibilities) <- c(length(Possibilities)/(sum(ToKeep | NoChange)),
                          sum(ToKeep | NoChange))

  PathPenality <- function(ChangeNodes, InitialStage, Mode) {

    Sphases <- rep(InitialStage, ncol(NormStageMatrix))

    tStart <- NA
    tEnd <- NA

    for (i in 1:(length(ChangeNodes)-1)) {

      if(!is.na(ChangeNodes[i])){
        tStart <- ChangeNodes[i]
      }

      tEnd <- ChangeNodes[i+1]

      if(is.na(tStart) | is.na(tEnd)){
        next()
      }

      Sphases[(tStart+1):(tEnd)] <- InitialStage + i
    }

    Sphases[Sphases>length(ChangeNodes)] <- Sphases[Sphases>length(ChangeNodes)] - length(ChangeNodes)

    if(Mode == 1){
      # Squared distance from the maximum
      return(
        sum(
          NormNodePenalty*(apply(NormStageMatrix, 2, sum) - diag(NormStageMatrix[Sphases,]))
        )
      )
    }

    if(Mode == 2){
      # Sum of "off-stage" contributions
      return(
        sum(
          NormNodePenalty*(apply(NormStageMatrix, 2, max) - mapply("[[", apply(NormStageMatrix, 2, as.list), Sphases))
        )
      )
    }

    if(Mode == 3){
      # Sum binary difference from maximum
      return(
        sum(
          NormNodePenalty*(apply(NormStageMatrix, 2, which.max) != Sphases)
        )
      )
    }

    if(Mode == 4){
      # Sum binary difference from maximum
      tDiff <- abs(apply(NormStageMatrix, 2, which.max) - Sphases)
      tDiff[tDiff > nrow(NormStageMatrix)/2] <- nrow(NormStageMatrix) - tDiff[tDiff > nrow(NormStageMatrix)/2]
      return(
        sum(
          NormNodePenalty*tDiff
        )
      )
    }

  }

  CombinedInfo <- NULL

  for (i in 1:nrow(NormStageMatrix)) {
    CombinedInfo <- cbind(CombinedInfo,
                          rbind(rep(i, sum(ToKeep | NoChange)),
                                apply(Possibilities, 2, PathPenality, InitialStage = i, Mode = Mode),
                                1:(sum(ToKeep | NoChange))
                          )
    )
  }

  return(list(Penality = CombinedInfo, Possibilities = Possibilities))

}




















# DataStruct = Data.Buet
# ProcStruct = Proc.Exp.Buet
# FiltMax = 0
# Thr = .7
# QuantSel = .5
# SinglePeack = TRUE
# Mode = "UP"
# MinNodes = 2
# StageInfo <- list(G1 = G1.All, S = S.All, G2M = G2M.All)
# StageInfo <- list(G1 = G1.Sel, S = S.Sel, G2M = G2M.Sel)
# StageInfo <- list(G0 = G0.All, G1 = G1.All, S = S.All, G2M = G2.All)

#' Title
#'
#' @param DataStruct
#' @param ProcStruct
#' @param StageInfo
#' @param MinNodes
#' @param FiltMax
#' @param Thr
#' @param QuantSel
#' @param Mode
#'
#' @return
#' @export
#'
#' @examples
StageWithPeaks <- function(DataStruct,
                           ProcStruct,
                           StageInfo.UP,
                           StageInfo.DOWN,
                           ComputeG0 = TRUE,
                           CCGenes = NULL,
                           MinNodes = 2,
                           FiltMax = .2,
                           Thr = .9,
                           QuantSel = .75,
                           MedianCVThr = 1,
                           G0Level = .6,
                           # SinglePeack = FALSE,
                           Mode = 1,
                           Title = ''
                              ) {

  WorkPG <- DataStruct$Analysis$PGStructs[[length(DataStruct$Analysis$PGStructs)]]
  TaxList <- WorkPG$TaxonList[[length(WorkPG$TaxonList)]]
  TaxVect <- rep(NA, nrow(WorkPG$Data))
  names(TaxVect) <- rownames(WorkPG$Data)

  lapply(1:length(TaxList), function(j){
    TaxVect[ TaxList[[j]] ] <<- j
  })

  GeneExprMat.DF <- data.frame(t(DataStruct$ExpMat[, names(ProcStruct$CellsPT)]))
  colnames(GeneExprMat.DF) <- rownames(DataStruct$ExpMat)

  # SelCell <- names(TaxVect)

  SelCell <- intersect(names(TaxVect), rownames(GeneExprMat.DF))

  GeneExprMat.DF.Split <- split(GeneExprMat.DF[SelCell, ], TaxVect[SelCell])
  GeneExprMat.DF.Split.Mean <- lapply(GeneExprMat.DF.Split, colMeans)
  GeneExprMat.DF.Split.Sd <- lapply(GeneExprMat.DF.Split, function(x){apply(x, 2, sd)})

  GeneExprMat.DF.Split.Mean.Bind <- sapply(GeneExprMat.DF.Split.Mean, cbind)
  GeneExprMat.DF.Split.Sd.Bind <- sapply(GeneExprMat.DF.Split.Sd, cbind)

  rownames(GeneExprMat.DF.Split.Mean.Bind) <- colnames(GeneExprMat.DF)
  colnames(GeneExprMat.DF.Split.Mean.Bind) <- names(GeneExprMat.DF.Split.Mean)

  rownames(GeneExprMat.DF.Split.Sd.Bind) <- colnames(GeneExprMat.DF)
  colnames(GeneExprMat.DF.Split.Sd.Bind) <- names(GeneExprMat.DF.Split.Sd)

  Reord <- ProcStruct$ExtPath[-length(ProcStruct$ExtPath)]
  Reord <- Reord[Reord %in% colnames(GeneExprMat.DF.Split.Mean.Bind)]

  GeneExprMat.DF.Split.Mean.Bind <- GeneExprMat.DF.Split.Mean.Bind[,Reord]
  GeneExprMat.DF.Split.Sd.Bind <- GeneExprMat.DF.Split.Sd.Bind[,Reord]

  MedianCV <- apply(GeneExprMat.DF.Split.Sd.Bind/GeneExprMat.DF.Split.Mean.Bind, 1, median, na.rm=TRUE)

  # NormExp <- ProcStruct$NodesExp

  NormExp <- GeneExprMat.DF.Split.Mean.Bind
  dim(NormExp)

  MedianCV <- MedianCV[rowSums(NormExp>0) > MinNodes]
  NormExp <- NormExp[rowSums(NormExp>0) > MinNodes, ]
  dim(NormExp)

  SDVect <- apply(NormExp, 1, sd)
  SDVect[is.na(SDVect)] <- 0
  MedianCV[is.na(MedianCV)] <- 0

  NormExp <- NormExp[SDVect > quantile(SDVect, QuantSel) & MedianCV < MedianCVThr,]
  dim(NormExp)

  NormExp[NormExp <= FiltMax] <- 0
  NormExp <- NormExp[rowSums(NormExp>0) > MinNodes, ]

  AVGenAct_NoNorm <- apply(NormExp[rownames(NormExp) %in% CCGenes,], 2, mean)
  AVGenAct_NoNorm.Comp <- apply(NormExp[!(rownames(NormExp) %in% CCGenes),], 2, mean)

  NormExp <- (NormExp - apply(NormExp, 1, min))/(apply(NormExp, 1, max) - apply(NormExp, 1, min))









  temp <- lapply(1:length(StageInfo.UP), function(i){

    if(is.null(StageInfo.UP[[i]])){
      return(NULL)
    }

    if(sum(rownames(NormExp) %in% StageInfo.UP[[i]]) == 0){
      return(NULL)
    }

    pheatmap::pheatmap(1*NormExp[rownames(NormExp) %in% StageInfo.UP[[i]],],
                       show_colnames = TRUE,
                       show_rownames = FALSE,
                       cluster_cols = FALSE,
                       cluster_rows = TRUE,
                       main = paste(names(StageInfo.UP)[i], "UP"))
  })








  temp <- lapply(1:length(StageInfo.DOWN), function(i){

    if(is.null(StageInfo.DOWN[[i]])){
      return(NULL)
    }

    if(sum(rownames(NormExp) %in% StageInfo.DOWN[[i]]) == 0){
      return(NULL)
    }

    pheatmap::pheatmap(1*NormExp[rownames(NormExp) %in% StageInfo.DOWN[[i]],],
                       show_colnames = TRUE,
                       show_rownames = FALSE,
                       cluster_cols = FALSE,
                       cluster_rows = TRUE,
                       main = paste(names(StageInfo.DOWN)[i], "DOWN"))
  })





  # PercOnNodes.UP <- sapply(StageInfo.UP, function(x){
  #   # apply(NormExp.Mod[rownames(NormExp.Mod) %in% intersect(x, Selected),], 2, sum)
  #   # apply(NormExp[rownames(NormExp) %in% x,], 2, sum)
  #   tMat <- NormExp[rownames(NormExp) %in% x,]
  #   tMat[tMat < Thr] <- 0
  #   apply(tMat, 2, mean)
  # })
  #
  #
  # PercOnNodes.DOWN <- sapply(StageInfo.DOWN, function(x){
  #   # apply(NormExp.Mod[rownames(NormExp.Mod) %in% intersect(x, Selected),], 2, sum)
  #   # apply(NormExp[rownames(NormExp) %in% x,], 2, sum)
  #   tMat <- 1 - NormExp[rownames(NormExp) %in% x,]
  #   tMat[tMat < Thr] <- 0
  #   apply(tMat, 2, mean)
  # })


  PercOnNodes <- sapply(1:length(StageInfo.UP), function(i){

    if(sum(rownames(NormExp) %in% union(StageInfo.UP[[i]], StageInfo.DOWN[[i]])) <= 1){
      return(rep(0, ncol(NormExp)))
    }

    tMat.UP <- NormExp[rownames(NormExp) %in% StageInfo.UP[[i]],]
    tMat.UP[tMat.UP < Thr] <- 0

    tMat.DOWN <- 1 - NormExp[rownames(NormExp) %in% StageInfo.DOWN[[i]],]
    tMat.DOWN[tMat.DOWN < Thr] <- 0

    CombMat <- rbind(tMat.DOWN, tMat.UP)
    colnames(CombMat) <- colnames(tMat.UP)

    apply(CombMat, 2, mean)
  })

  colnames(PercOnNodes) <- names(StageInfo.UP)

  barplot(t(PercOnNodes), beside = TRUE)

  PercOnNodes <- PercOnNodes[, !apply(PercOnNodes == 0, 2, all)]

  if(ComputeG0){

    AVGenAct <- apply(NormExp[rownames(NormExp) %in% CCGenes,], 2, mean)
    # AVGenAct.Comp <- apply(NormExp[!(rownames(NormExp) %in% CCGenes),], 2, mean)

    barplot(rbind(AVGenAct_NoNorm, AVGenAct_NoNorm.Comp), beside = TRUE)
    abline(h = G0Level)

    ToKeep <- AVGenAct_NoNorm > G0Level

    PercOnNodes.Comb <- t(PercOnNodes)
    PercOnNodes.Comb[] <- 0
    PercOnNodes.Comb[,ToKeep] <- t(PercOnNodes[ToKeep, ])/apply(PercOnNodes[ToKeep, ], 2, quantile, .9)

    # G0Val <- min(apply(PercOnNodes.Comb[,ToKeep], 2, max))
    # PercOnNodes.Comb <- rbind(G0Val*(!ToKeep), PercOnNodes.Comb)

    # PercOnNodes.Comb <- rbind((1 - AVGenAct)*(!ToKeep), PercOnNodes.Comb)

    PercOnNodes.Comb <- rbind(
      (1 - AVGenAct)/quantile(1 - AVGenAct, .9),
      PercOnNodes.Comb
    )

    PercOnNodes.Comb[1, ToKeep] <- 0
    rownames(PercOnNodes.Comb)[1] <- "G0"

    barplot(PercOnNodes.Comb, beside = TRUE)

    OutPos <- apply(PercOnNodes.Comb[, ToKeep], 1, scater::isOutlier) %>%
      apply(., 1, any)

    OutVect <- ToKeep
    OutVect[OutVect] <- OutPos

    NonOutMax <- apply(PercOnNodes.Comb[-1,!OutVect & ToKeep], 1, max)
    NonOutMax[NonOutMax == 0] <- 1

    NormVect <- sapply(1:length(OutVect), function(i){
      if(!OutVect[i]){
        return(1)
      } else {
        return(max(PercOnNodes.Comb[-1, i]/NonOutMax))
      }
    })

    PercOnNodes.Comb <- t(t(PercOnNodes.Comb)/NormVect)

    barplot(PercOnNodes.Comb, beside = TRUE)

  } else {

    PercOnNodes.Comb <- t(PercOnNodes)/apply(PercOnNodes, 2, quantile, .9)

    barplot(PercOnNodes.Comb, beside = TRUE)

  }



  # G0Perc <- colSums(!PlotMat[rownames(PlotMat) %in% unlist(StageInfo, use.names = FALSE),])/sum(rownames(PlotMat) %in% unlist(StageInfo, use.names = FALSE))

  # PercOnNodes.UP <- t(PercOnNodes.UP)/apply(PercOnNodes.UP, 2, quantile, .9)
  # PercOnNodes.DOWN <- t(PercOnNodes.DOWN)/apply(PercOnNodes.DOWN, 2, quantile, .9)
  #
  # PercOnNodes.UP[is.na(PercOnNodes.UP)] <- 0
  # PercOnNodes.DOWN[is.na(PercOnNodes.DOWN)] <- 0

  # barplot(PercOnNodes, beside = TRUE)




  # PercOnNodes.Comb <- t(PercOnNodes)


  # PercOnNodes <- PercOnNodes + min(PercOnNodes)

  # PercOnNodes[PercOnNodes > 1] <- 1



  # LowHigh <- apply(PercOnNodes, 1, quantile, .7)
  #
  # # barplot(PercOnNodes, beside = TRUE)
  #
  # pheatmap::pheatmap((PercOnNodes > LowHigh)*1, cluster_cols = FALSE, cluster_rows = FALSE)
  #
  # lapply(1:nrow(PercOnNodes), function(i){
  #
  #   PercOnNodes[i,] > LowHigh[1,i]
  #
  # })

  pheatmap::pheatmap(PercOnNodes.Comb,
                     cluster_cols = FALSE,
                     cluster_rows = FALSE)

  apply(PercOnNodes.Comb, 2, which.max)

  # PercOnNodes[PercOnNodes < .02] <- 0
  #
  # pheatmap::pheatmap(PercOnNodes,
  #                    cluster_cols = FALSE,
  #                    cluster_rows = FALSE)
  #
  # apply(PercOnNodes, 2, which.max)

  BestStaging <- AssociteNodes(PerMat = PercOnNodes.Comb, Mode = Mode)

  TB <- apply(BestStaging, 2, function(x){table(factor(x, levels = 1:nrow(PercOnNodes.Comb)))})

  rownames(TB) <- rownames(PercOnNodes.Comb)

  pheatmap::pheatmap(TB, cluster_cols = FALSE, cluster_rows = FALSE,
                     main = Title)

  BestTB <- TB

  InferredStages <- rownames(PercOnNodes.Comb)[apply(TB, 2, which.max)]
  names(InferredStages) <- colnames(PercOnNodes.Comb)

  CellStages <- InferredStages[paste(TaxVect)]
  names(CellStages) <- names(TaxVect)

  if(ComputeG0){
    CellStages <- factor(CellStages, levels = c('G0', names(StageInfo.UP)))
  } else {
    CellStages <- factor(CellStages, levels = names(StageInfo.UP))
  }


  if(!is.null(names(DataStruct$Cats))){
    TB1 <- table(DataStruct$Cats[names(CellStages)], CellStages)
  } else {
    TB1 <- table(DataStruct$Cats, CellStages)
  }

  print(TB1/rowSums(TB1))

  pheatmap::pheatmap(TB1/rowSums(TB1), cluster_rows = FALSE, cluster_cols = FALSE, main = Title,
                     color = rev(heat.colors(25)))

  if(nrow(TB1) > 1){

    print(t(t(TB1)/colSums(TB1)))
    pheatmap::pheatmap(t(t(TB1)/colSums(TB1)), cluster_rows = FALSE, cluster_cols = FALSE, main = Title,
                       color = rev(heat.colors(25)))

  }


  ExtStages <- c(InferredStages[length(InferredStages)], InferredStages, InferredStages[1])
  CombStages <- rep(NA, length(InferredStages))
  names(CombStages) <- names(InferredStages)

  for(i in 2:(length(ExtStages)-1)){

    if(ExtStages[i-1] == ExtStages[i] &
       ExtStages[i+1] == ExtStages[i]){
      CombStages[i-1] <- ExtStages[i]
      next()
    }

    if(ExtStages[i-1] == ExtStages[i] &
       ExtStages[i+1] != ExtStages[i]){
      CombStages[i-1] <- paste(ExtStages[i], ExtStages[i+1], sep = " / ")
      next()
    }

    if(ExtStages[i-1] != ExtStages[i] &
       ExtStages[i+1] == ExtStages[i]){
      CombStages[i-1] <- paste(ExtStages[i-1], ExtStages[i], sep = " / ")
      next()
    }

    if(ExtStages[i-1] != ExtStages[i] &
       ExtStages[i+1] != ExtStages[i]){
      CombStages[i-1] <- paste(ExtStages[i-1], ExtStages[i], ExtStages[i+1], sep = " / ")
      next()
    }

  }


  CellStages_Ext <- CombStages[as.character(TaxVect)]
  names(CellStages_Ext) <- names(TaxVect)
  # CellStages_Ext <- factor(CellStages_Ext)


  if(!is.null(names(DataStruct$Cats))){
    TB2 <- table(DataStruct$Cats[names(CellStages_Ext)], CellStages_Ext)
  } else {
    TB2 <- table(DataStruct$Cats, CellStages_Ext)
  }

  print(TB2/rowSums(TB2))

  pheatmap::pheatmap(TB2/rowSums(TB2), cluster_rows = FALSE, cluster_cols = FALSE, main = Title,
                     color = rev(heat.colors(25)))

  if(nrow(TB2) > 1){

    print(t(t(TB2)/colSums(TB2)))
    pheatmap::pheatmap(t(t(TB2)/colSums(TB2)), cluster_rows = FALSE, cluster_cols = FALSE, main = Title,
                       color = rev(heat.colors(25)))


  }


  return(list(Inferred = InferredStages,
              Extinferred = CombStages,
              CellStages = CellStages,
              CellStages_Ext = CellStages_Ext,
              StageMat = PercOnNodes.Comb,
              TB1 = TB1, TB2 = TB2,
              BestTB = BestTB))

}








#' Title
#'
#' @param PerMat
#'
#' @return
#' @export
#'
#' @examples
AssociteNodes <- function(PerMat, Mode = 4) {

  tictoc::tic()
  print("Direct staging")
  Staging <- FitStagesCirc(StageMatrix = PerMat,
                           NodePenalty = rep(1, ncol(PerMat)),
                           Mode = Mode)
  tictoc::toc()

  tictoc::tic()
  print("Reverse staging")
  StagingRev <- FitStagesCirc(StageMatrix = PerMat[, rev(1:ncol(PerMat))],
                              NodePenalty = rep(1, ncol(PerMat)),
                              Mode = Mode)
  tictoc::toc()

  AllPenality <- rbind(cbind(Staging$Penality, StagingRev$Penality), rep(1:2, each=ncol(Staging$Penality)))

  Idxs <- which(AllPenality[2, ] == min(AllPenality[2, ]))

  SelPenality <- AllPenality[,Idxs]
  dim(SelPenality) <- c(4, length(SelPenality)/4)

  DirectPenality <- NULL
  DirectChanges <- NULL

  if(sum(SelPenality[4,] == 1) > 0){

    ExpandStages <- function(idx) {

      ChangeNodes <- Staging$Possibilities[ , SelPenality[3, idx]]

      StageVect <- rep(SelPenality[1, idx], ncol(PerMat))

      tStart <- NA
      tEnd <- NA

      for (i in 1:(length(ChangeNodes)-1)) {

        if(!is.na(ChangeNodes[i])){
          tStart <- ChangeNodes[i]
        }

        tEnd <- ChangeNodes[i+1]

        if(is.na(tStart) | is.na(tEnd)){
          next()
        }

        StageVect[(tStart+1):(tEnd)] <- SelPenality[1, idx] + i
      }

      StageVect[StageVect>length(ChangeNodes)] <- StageVect[StageVect>length(ChangeNodes)] - length(ChangeNodes)

      return(StageVect)

    }

    SelPenIdx <- which(SelPenality[4,]==1)

    DirectPenality <- SelPenality[2,SelPenIdx]
    DirectChanges <- t(sapply(SelPenIdx, ExpandStages))

  }

  ReversePenality <- NULL
  ReverseChanges <- NULL

  if(sum(SelPenality[4,] == 2) > 0){

    ExpandStages <- function(idx) {

      ChangeNodes <- StagingRev$Possibilities[ , SelPenality[3, idx]]

      StageVect <- rep(SelPenality[1, idx], ncol(PerMat))

      tStart <- NA
      tEnd <- NA

      for (i in 1:(length(ChangeNodes)-1)) {

        if(!is.na(ChangeNodes[i])){
          tStart <- ChangeNodes[i]
        }

        tEnd <- ChangeNodes[i+1]

        if(is.na(tStart) | is.na(tEnd)){
          next()
        }

        StageVect[(tStart+1):(tEnd)] <- SelPenality[1, idx] + i
      }

      StageVect[StageVect>length(ChangeNodes)] <- StageVect[StageVect>length(ChangeNodes)] - length(ChangeNodes)

      return(rev(StageVect))

    }

    SelPenIdx <- which(SelPenality[4,]==2)

    ReversePenality <- SelPenality[2,SelPenIdx]
    ReverseChanges <- t(sapply(SelPenIdx, ExpandStages))

  }

  AllStg <- rbind(DirectChanges, ReverseChanges)
  AllPen <- c(DirectPenality, ReversePenality)
  AllDir <- c(rep("Dir", length(DirectPenality)),
              rep("Rev", length(ReversePenality)))
  colnames(AllStg) <- colnames(PerMat)

  return(AllStg)

}







