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
                              SinglePeack = FALSE,
                              Mode = "UP") {

  WorkPG <- DataStruct$Analysis$PGStructs[[length(DataStruct$Analysis$PGStructs)]]
  TaxList <- WorkPG$TaxonList[[length(WorkPG$TaxonList)]]
  TaxVect <- rep(NA, nrow(WorkPG$Data))
  names(TaxVect) <- rownames(WorkPG$Data)

  lapply(1:length(TaxList), function(j){
    TaxVect[ TaxList[[j]] ] <<- j
  })


  if(AllGenes){

    GeneExprMat.DF <- data.frame(t(DataStruct$ExpMat[, names(ProcStruct$CellsPT)]))
    colnames(GeneExprMat.DF) <- rownames(DataStruct$ExpMat)

    SelCell <- names(TaxVect)

    GeneExprMat.DF.Split <- split(GeneExprMat.DF[SelCell, ], TaxVect[SelCell])
    GeneExprMat.DF.Split.Mean <- lapply(GeneExprMat.DF.Split, colMeans)

    GeneExprMat.DF.Split.Mean.Bind <- sapply(GeneExprMat.DF.Split.Mean, cbind)
    rownames(GeneExprMat.DF.Split.Mean.Bind) <- colnames(GeneExprMat.DF)
    colnames(GeneExprMat.DF.Split.Mean.Bind) <- names(GeneExprMat.DF.Split.Mean)

    Reord <- ProcStruct$ExtPath[-length(ProcStruct$ExtPath)]
    Reord <- Reord[Reord %in% colnames(GeneExprMat.DF.Split.Mean.Bind)]

    GeneExprMat.DF.Split.Mean.Bind <- GeneExprMat.DF.Split.Mean.Bind[,Reord]

    # NormExp <- ProcStruct$NodesExp

    NormExp <- GeneExprMat.DF.Split.Mean.Bind
    dim(NormExp)

  } else {

    NormExp <- ProcStruct$NodesExp
    colnames(NormExp) <- ProcStruct$ExtPath

  }



  NormExp <- NormExp[rowSums(NormExp>0) > MinNodes, ]
  dim(NormExp)

  SDVect <- apply(NormExp, 1, sd)
  NormExp <- NormExp[SDVect > quantile(SDVect, QuantSel),]
  dim(NormExp)


  # CorGenes <- cor(t(NormExp))
  #
  # CorGenes[CorGenes < .95] <- 0
  # CorGenes[CorGenes >= .95] <- 1
  #
  # CorGenes <- CorGenes[rowSums(CorGenes) > 1, colSums(CorGenes) > 1]
  #
  # dim(CorGenes)
  #
  # pheatmap::pheatmap(CorGenes)
  #
  # Modules <- apply(CorGenes==1, 1, which)
  #
  # Modules <- Modules[unlist(lapply(Modules, length)) > 3]
  #
  # ModMat <- lapply(Modules, function(x){rownames(CorGenes) %in% x})
  #
  # ModMat <- t(sapply(ModMat, cbind))
  #
  # dim(ModMat)

  NormExp[NormExp <= FiltMax] <- 0
  NormExp <- NormExp[rowSums(NormExp>0) > MinNodes, ]

  if(Mode == "UP"){
    Min <- apply(NormExp, 1, min)
    Range <- apply(NormExp, 1, max) - Min
    NormExp.Bin <- (NormExp > Min + Range*Thr)
  } else {
    Min <- apply(NormExp, 1, min)
    Range <- apply(NormExp, 1, max) - Min
    NormExp.Bin <- (NormExp < Min + Range*Thr)
  }


  NormExp.Bin <- NormExp.Bin[apply(NormExp.Bin, 1, any), ]

  pheatmap::pheatmap(1*NormExp.Bin[order(rowMeans(t(t(NormExp.Bin)*1:ncol(NormExp.Bin)))), ],
                     cluster_cols = FALSE, cluster_rows = FALSE)

  if(SinglePeack){

    SinglePk_UP <- lapply(apply(NormExp.Bin, 1, which), function(x){
      if(all(min(x):max(x) %in% x)){
        return(TRUE)
      } else {
        return(FALSE)
      }
    })

    SinglePk_DOWN <- lapply(apply(!NormExp.Bin, 1, which), function(x){
      if(all(min(x):max(x) %in% x)){
        return(TRUE)
      } else {
        return(FALSE)
      }
    })

    Selected <- c(names(which(unlist(SinglePk_DOWN))),
                  names(which(unlist(SinglePk_UP))))

  } else {

    Selected <- rownames(NormExp.Bin)

  }



  Selected <- unique(Selected)

  PlotMat <- 1*NormExp.Bin[Selected,]

  PlotMat.Numb <- t(t(PlotMat)*1:ncol(PlotMat))
  PlotMat.Numb[PlotMat.Numb == 0] <- NA

  pheatmap::pheatmap(PlotMat[order(apply(PlotMat.Numb, 1, median, na.rm=TRUE)),], show_colnames = FALSE,
                     show_rownames = FALSE, cluster_cols = FALSE, cluster_rows = FALSE)

  barplot(colSums(NormExp.Bin[Selected,]))

  StageGenes <- lapply(ProcStruct$StageOnNodes, function(x){
    SelX <- intersect(names(x), colnames(NormExp.Bin))
    if(length(SelX)>0){
      if(length(SelX)==1){
        NormExp.Bin[Selected, SelX]
      } else {
        apply(NormExp.Bin[Selected, SelX], 1, function(y){
          any(which(y))
        })
      }
    } else {
      NA
    }

  })

  StageGenes.Names <- lapply(StageGenes, function(x){
    names(which(x))
  })

  barplot(unlist(lapply(StageGenes.Names, length)), las = 2, horiz = FALSE)

  return(StageGenes.Names)

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
                           StageInfo,
                           ComputeG0 = TRUE,
                           MinNodes = 2,
                           FiltMax = .2,
                           Thr = .9,
                           QuantSel = .75,
                           G0Level = 1,
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

  SelCell <- names(TaxVect)

  GeneExprMat.DF.Split <- split(GeneExprMat.DF[SelCell, ], TaxVect[SelCell])
  GeneExprMat.DF.Split.Mean <- lapply(GeneExprMat.DF.Split, colMeans)

  GeneExprMat.DF.Split.Mean.Bind <- sapply(GeneExprMat.DF.Split.Mean, cbind)
  rownames(GeneExprMat.DF.Split.Mean.Bind) <- colnames(GeneExprMat.DF)
  colnames(GeneExprMat.DF.Split.Mean.Bind) <- names(GeneExprMat.DF.Split.Mean)

  Reord <- ProcStruct$ExtPath[-length(ProcStruct$ExtPath)]
  Reord <- Reord[Reord %in% colnames(GeneExprMat.DF.Split.Mean.Bind)]

  GeneExprMat.DF.Split.Mean.Bind <- GeneExprMat.DF.Split.Mean.Bind[,Reord]

  # NormExp <- ProcStruct$NodesExp

  NormExp <- GeneExprMat.DF.Split.Mean.Bind
  dim(NormExp)

  NormExp <- NormExp[rowSums(NormExp>0) > MinNodes, ]
  dim(NormExp)

  SDVect <- apply(NormExp, 1, sd)
  NormExp <- NormExp[SDVect > quantile(SDVect, QuantSel),]
  dim(NormExp)

  NormExp[NormExp <= FiltMax] <- 0
  NormExp <- NormExp[rowSums(NormExp>0) > MinNodes, ]


#
#   pheatmap::pheatmap(NormExp[rownames(NormExp) %in% G1.buet,],
#                      cluster_cols = FALSE,
#                      cluster_rows = FALSE)


  # if(Mode == "UP"){
  #   Min <- apply(NormExp, 1, min)
  #   Range <- apply(NormExp, 1, max) - Min
  #   NormExp.Bin <- (NormExp > Min + Range*Thr)
  # } else {
  #   Min <- apply(NormExp, 1, min)
  #   Range <- apply(NormExp, 1, max) - Min
  #   NormExp.Bin <- (NormExp < Min + Range*Thr)
  # }

  Min <- apply(NormExp, 1, min)
  Max <- apply(NormExp, 1, max)
  # NormExp.Norm <- (NormExp - Min)/(Max-Min)

  # NormExp.Bin <- NormExp.Bin[apply(NormExp.Bin, 1, any), ]
  #
  # pheatmap::pheatmap(1*NormExp.Bin[order(rowMeans(t(t(NormExp.Bin)*1:ncol(NormExp.Bin)))), ],
  #                    cluster_cols = FALSE, cluster_rows = FALSE)

  pheatmap::pheatmap(NormExp, cluster_cols = FALSE, cluster_rows = FALSE)


  # if(SinglePeack){
  #
  #   SinglePk_UP <- lapply(apply(NormExp.Bin, 1, which), function(x){
  #     if(all(min(x):max(x) %in% x)){
  #       return(TRUE)
  #     } else {
  #       return(FALSE)
  #     }
  #   })
  #
  #   SinglePk_DOWN <- lapply(apply(!NormExp.Bin, 1, which), function(x){
  #     if(all(min(x):max(x) %in% x)){
  #       return(TRUE)
  #     } else {
  #       return(FALSE)
  #     }
  #   })
  #
  #   Selected <- c(names(which(unlist(SinglePk_DOWN))),
  #                 names(which(unlist(SinglePk_UP))))
  #
  # } else {
  #
  #   Selected <- rownames(NormExp.Bin)
  #
  # }
  #
  #
  #
  # Selected <- unique(Selected)
  #
  # PlotMat <- 1*NormExp.Bin[Selected,]
  #
  # PlotMat.Numb <- t(t(PlotMat)*1:ncol(PlotMat))
  # PlotMat.Numb[PlotMat.Numb == 0] <- NA
  #
  # pheatmap::pheatmap(PlotMat[order(apply(PlotMat.Numb, 1, median, na.rm=TRUE)),], show_colnames = FALSE,
  #                    show_rownames = FALSE, cluster_cols = FALSE, cluster_rows = FALSE)
  #
  # barplot(colSums(NormExp.Bin[Selected,]))

  # temp <- lapply(1:length(StageInfo), function(i){
  #   pheatmap::pheatmap(PlotMat[rownames(PlotMat) %in% StageInfo[[i]],],
  #                      show_colnames = TRUE,
  #                      show_rownames = FALSE,
  #                      cluster_cols = FALSE,
  #                      cluster_rows = FALSE,
  #                      main = names(StageInfo)[i])
  # })

  temp <- lapply(1:length(StageInfo), function(i){
    pheatmap::pheatmap(NormExp[rownames(NormExp) %in% StageInfo[[i]],],
                       show_colnames = TRUE,
                       show_rownames = FALSE,
                       cluster_cols = FALSE,
                       cluster_rows = FALSE,
                       main = names(StageInfo)[i])
  })

  # PercOnNodes <- sapply(StageInfo, function(x){
  #   colSums(PlotMat[rownames(PlotMat) %in% x,])/sum(rownames(PlotMat) %in% x)
  # })

  # PercOnNodes <- sapply(StageInfo, function(x){
  #   apply(NormExp.Norm[rownames(NormExp.Norm) %in% x,], 2, mean)
  # })

  PercOnNodes <- sapply(StageInfo, function(x){
    apply(NormExp[rownames(NormExp) %in% x,], 2, sum)
  })

  # G0Perc <- colSums(!PlotMat[rownames(PlotMat) %in% unlist(StageInfo, use.names = FALSE),])/sum(rownames(PlotMat) %in% unlist(StageInfo, use.names = FALSE))

  PercOnNodes <- t(PercOnNodes)

  PercOnNodes <- PercOnNodes/apply(PercOnNodes, 1, mean)

  if(ComputeG0){
    G0Perc <- rep(G0Level, ncol(PercOnNodes))
    G0Perc[G0Perc < 0] <- 0

    PercOnNodes <- rbind(G0Perc, PercOnNodes)
    rownames(PercOnNodes)[1] <- "G0"
  }

  # PercOnNodes[PercOnNodes > 1] <- 1

  barplot(PercOnNodes, beside = TRUE)

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

  pheatmap::pheatmap(PercOnNodes,
                     cluster_cols = FALSE,
                     cluster_rows = FALSE)

  apply(PercOnNodes, 2, which.max)

  # PercOnNodes[PercOnNodes < .02] <- 0
  #
  # pheatmap::pheatmap(PercOnNodes,
  #                    cluster_cols = FALSE,
  #                    cluster_rows = FALSE)
  #
  # apply(PercOnNodes, 2, which.max)

  BestStaging <- AssociteNodes(PerMat = PercOnNodes, Mode = Mode)

  TB <- apply(BestStaging, 2, function(x){table(factor(x, levels = 1:nrow(PercOnNodes)))})

  rownames(TB) <- rownames(PercOnNodes)

  pheatmap::pheatmap(TB, cluster_cols = FALSE, cluster_rows = FALSE,
                     main = Title)

  BestTB <- TB

  InferredStages <- rownames(PercOnNodes)[apply(TB, 2, which.max)]
  names(InferredStages) <- colnames(PercOnNodes)

  CellStages <- InferredStages[paste(TaxVect)]
  names(CellStages) <- names(TaxVect)

  if(ComputeG0){
    CellStages <- factor(CellStages, levels = c('G0', names(StageInfo)))
  } else {
    CellStages <- factor(CellStages, levels = names(StageInfo))
  }


  if(!is.null(names(DataStruct$Cats))){
    TB1 <- table(DataStruct$Cats[names(CellStages)], CellStages)
  } else {
    TB1 <- table(DataStruct$Cats, CellStages)
  }

  print(TB1/rowSums(TB1))

  pheatmap::pheatmap(TB1/rowSums(TB1), cluster_rows = FALSE, cluster_cols = FALSE, main = Title)

  if(nrow(TB1) > 1){

    print(t(t(TB1)/colSums(TB1)))
    pheatmap::pheatmap(t(t(TB1)/colSums(TB1)), cluster_rows = FALSE, cluster_cols = FALSE, main = Title)

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

  pheatmap::pheatmap(TB2/rowSums(TB2), cluster_rows = FALSE, cluster_cols = FALSE, main = Title)

  if(nrow(TB2) > 1){

    print(t(t(TB2)/colSums(TB2)))
    pheatmap::pheatmap(t(t(TB2)/colSums(TB2)), cluster_rows = FALSE, cluster_cols = FALSE, main = Title)


  }


  return(list(Inferred = InferredStages,
              Extinferred = CombStages,
              CellStages = CellStages,
              CellStages_Ext = CellStages_Ext,
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







