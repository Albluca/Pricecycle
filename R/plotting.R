#################################################################################
# PlotOnStages ------------------------------------------------------
#
################################################################################


#' Title
#'
#' @param Structure
#' @param TaxonList
#' @param Categories
#' @param PrinGraph
#' @param Net
#' @param SelThr
#' @param ComputeOverlaps
#' @param RotatioMatrix
#' @param ExpData
#'
#' @return
#' @export
#'
#' @examples
PlotOnStages <- function(Structure, TaxonList, Categories, PrinGraph, Net, SelThr = NULL,
                         ComputeOverlaps = TRUE, RotatioMatrix, PointProjections,
                         ExpData, nGenes = 10, OrderOnCat = FALSE, SmoothPoints = 1,
                         PCACenter = FALSE, MinCellPerNode = 3, Title = '') {

  if(!is.factor(Categories)){
    stop("Categories must be a factor")
  }

  if(is.null(SelThr)){
    SelThr <- 1.1*(1/length(levels(Categories)))
  }

  # Structure <- "Circle"
  # TaxonList <- BuettInfo$FinalStruct$TaxonList[[length(BuettInfo$FinalStruct$TaxonList)]]
  # Categories <- BuettInfo$FinalStruct$Categories
  # Net <- BuettInfo$FinalStruct$Net[[length(BuettInfo$FinalStruct$Net)]]
  # RotatioMatrix <- BuettInfo$FinalStruct$PCAData$rotation[,1:BuettInfo$FinalStruct$nDims]
  # PointProjections = BuettInfo$FinalStruct$ProjPoints[[length(BuettInfo$FinalStruct$ProjPoints)]]
  # ExpData <- BuettInfo$FinalStruct$FiltExp

  Nodes <- PrinGraph$Nodes
  nNodes <- nrow(PrinGraph$Nodes)

  # Find the best Staging if circular

  TaxVect <- rep(NA, nNodes)

  for(i in 1:length(TaxonList)){
    TaxVect[TaxonList[[i]]] <- i
  }

  TaxVect <- factor(TaxVect, levels = paste(1:nNodes))

  TB <- table(Categories, TaxVect)
  # colnames(TB)

  print("Step I - Finding the best path")

  StepIDone <- FALSE

  if(Structure == 'Circle'){

    StepIDone <- TRUE

    print("Getting all circular subisomorphisms")

    AllPaths <- GetLongestPath(Net = Net, Structure = Structure, Circular = TRUE)

    if(OrderOnCat){

      SummInfo <- NULL

      for(i in 1:nrow(AllPaths$VertNumb)){
        SelPath <- AllPaths$VertNumb[i,]
        # SelPath <- SelPath[SelPath %in% unique(TaxVect)]
        SelPath <- SelPath[-length(SelPath)]

        Reordered <- TaxVect
        for(j in 1:length(SelPath)){
          Reordered[as.character(TaxVect) == SelPath[j]] <- j
        }

        NumReord <- as.numeric(as.character(Reordered))

        AGG <- aggregate(NumReord, by = list(Categories), median)
        AGG2 <- aggregate(NumReord, by = list(Categories), min)
        AGG3 <- aggregate(NumReord, by = list(Categories), max)

        Sorted <- FALSE

        if(S4Vectors::isSorted(AGG[,2])){
          SummInfo <- rbind(SummInfo,
                            c(i, 1, summary(aov(NumReord ~ Categories))[[1]][1,"Pr(>F)"],
                              AGG2[1,2], AGG[1,2])
          )
          Sorted <- TRUE
          # boxplot(NumReord ~ Categories, main = i)
        }

        if(S4Vectors::isSorted(rev(AGG[,2]))){
          SummInfo <- rbind(SummInfo,
                            c(i, 2, summary(aov(NumReord ~ Categories))[[1]][1,"Pr(>F)"],
                              nNodes - AGG3[1,2] + 1,
                              nNodes - AGG[1,2] + 1)
          )
          Sorted <- TRUE
          # boxplot(NumReord ~ Categories, main = paste(i, "rev"))
        }

        if(!Sorted){

          if(cor(AGG[,2], 1:nrow(AGG))>0){
            SummInfo <- rbind(SummInfo,
                              c(i, 3, summary(aov(NumReord ~ Categories))[[1]][1,"Pr(>F)"],
                                AGG2[1,2], AGG[1,2])
            )
          } else {
            SummInfo <- rbind(SummInfo,
                              c(i, 4, summary(aov(NumReord ~ Categories))[[1]][1,"Pr(>F)"],
                                nNodes - AGG3[1,2] + 1,
                                nNodes - AGG[1,2] + 1)
            )
          }

        }
      }

      # print(SummInfo)

      print("Finding the best path")

      if(any(SummInfo[,2] %in% 1:2)){
        SummInfo <- SummInfo[SummInfo[,2] %in% 1:2, ]
      } else {
        print("Unable to find strongly consecutive stages")
      }

      if(length(SummInfo) == 5){
        SummInfo <- matrix(SummInfo, nrow = 1)
      }

      Selected <- SummInfo[SummInfo[, 3] == min(SummInfo[, 3]), ]

      print(Selected)

      if(length(Selected)>5){
        Selected <- Selected[which.min(Selected[,4]),]
      }

      print(Selected)

      SelPath <- AllPaths$VertNumb[Selected[1],]
      # SelPath <- SelPath[SelPath %in% unique(TaxVect)]
      SelPath <- SelPath[-length(SelPath)]

      if(Selected[2] %in% c(2, 4)){
        SelPath <- rev(SelPath)
      }

      Reordered <- TaxVect
      for(j in 1:length(SelPath)){
        Reordered[as.character(TaxVect) == SelPath[j]] <- j
      }

      NumReord <- as.numeric(as.character(Reordered))

      boxplot(NumReord ~ Categories, main = Title,
              ylab = "Position on path")

    } else {
      SelPath <- AllPaths$VertNumb[sample(x = 1:nrow(AllPaths$VertNumb), size = 1),]
      SelPath <- SelPath[-length(SelPath)]
    }


  }

  if(Structure == 'Line'){

    StepIDone <- TRUE

    print("Getting all line subisomorphisms")

    AllPaths <- GetLongestPath(Net = Net, Structure = Structure, Circular = FALSE)

    if(OrderOnCat){

      SummInfo <- NULL

      for(i in 1:nrow(AllPaths$VertNumb)){
        SelPath <- AllPaths$VertNumb[i,]
        # SelPath <- SelPath[SelPath %in% unique(TaxVect)]
        # SelPath <- SelPath[-length(SelPath)]

        Reordered <- TaxVect
        for(j in 1:length(SelPath)){
          Reordered[TaxVect == SelPath[j]] <- j
        }

        NumReord <- as.numeric(as.character(Reordered))

        AGG <- aggregate(NumReord, by = list(Categories), median)
        AGG2 <- aggregate(NumReord, by = list(Categories), min)
        AGG3 <- aggregate(NumReord, by = list(Categories), max)

        if(S4Vectors::isSorted(AGG[,2])){
          SummInfo <- rbind(SummInfo,
                            c(i, 1, summary(aov(NumReord ~ Categories))[[1]][1,"Pr(>F)"],
                              AGG2[1,2], AGG[1,2])
          )

          # boxplot(NumReord ~ Categories, main = i)

        }


        if(S4Vectors::isSorted(rev(AGG[,2]))){
          SummInfo <- rbind(SummInfo,
                            c(i, 2, summary(aov(NumReord ~ Categories))[[1]][1,"Pr(>F)"],
                              nNodes - AGG3[1,2] + 1,
                              nNodes - AGG[1,2] + 1)
          )
          # boxplot(Reordered ~ Categories, main = paste(i, "rev"))
        }
      }

      print("Finding the best path")

      Selected <- SummInfo[SummInfo[, 3] == min(SummInfo[, 3]), ]

      if(length(Selected)>5){
        Selected <- Selected[which.min(Selected[,4]),]
      }

      SelPath <- AllPaths$VertNumb[Selected[1],]
      # SelPath <- SelPath[SelPath %in% unique(TaxVect)]
      SelPath <- SelPath[-length(SelPath)]

      if(Selected[2] == 2){
        SelPath <- rev(SelPath)
      }

      Reordered <- TaxVect
      for(j in 1:length(SelPath)){
        Reordered[as.character(TaxVect) == SelPath[j]] <- j
      }

      NumReord <- as.numeric(as.character(Reordered))

      boxplot(NumReord ~ Categories)

    } else {
      SelPath <- AllPaths$VertNumb[sample(x = 1:nrow(AllPaths$VertNumb), size = 1),]
    }


  }

  if(!StepIDone){
    stop("Unsupported Structure")
  }


  print("Step II")

  Empty <- which(unlist(lapply(lapply(TaxonList, is.na), any)))

  ExtPath <- SelPath

  if(Structure == 'Circle'){
    ExtPath <- c(ExtPath, ExtPath[1])
  }

  ExtendedTB <- TB[,ExtPath]

  if(OrderOnCat){
    barplot(t(t(ExtendedTB)/colSums(ExtendedTB)), col = rainbow(nrow(TB)), beside = TRUE, las = 2,
            legend.text = rownames(ExtendedTB), args.legend = list(x = "top", fill=rainbow(nrow(TB))),
            ylim = c(0, 1.25), yaxt = "n")
    axis(2, seq(from=0, to=1, by=.25), las=2)
  }



  # SelPath <- SelPath[-length(SelPath)]
  SelPathSTG <- rep(NA, length(SelPath))


  PercMat <- t(t(TB)/colSums(TB))
  PercMat[is.na(PercMat)] <- 0
  PercMat[,colSums(TB) < MinCellPerNode] <- 0
  PercMat <- PercMat[,ExtPath]
  BinPercMat <- (PercMat > SelThr)

  # print(1*BinPercMat)

  # if(!any(BinPercMat[,1])){
    # BinPercMat[,1] <- BinPercMat[,2]
  # }

  if(OrderOnCat){

    # Fill gaps with TRUE

    if(Structure == 'Circle'){

      SelGrop <- apply(BinPercMat, 1, function(x){
        range(which(x))})

      for(j in 1:ncol(SelGrop)){
        if(any(is.infinite(SelGrop[,j]))){
          next
        } else {
          Diff <- SelGrop[1,j] + ncol(BinPercMat) + 1 - SelGrop[2,j] - 2

          if(Diff > 0 & Diff <= SmoothPoints){
            BinPercMat[colnames(SelGrop)[j],1:(SelGrop[1,j]-1)] <- TRUE
            BinPercMat[colnames(SelGrop)[j],(SelGrop[2,j]+1):ncol(BinPercMat)] <- TRUE
          }
        }
      }

    }

    WhichStage <- apply(BinPercMat, 1, which)

    WhichStage.Diff <- lapply(WhichStage, function(x){x[-1] - x[-length(x)]})

    lapply(as.list(1:length(WhichStage.Diff)), function(j){
      GapsToFill <- which(WhichStage.Diff[[j]] > 1 & WhichStage.Diff[[j]] <= SmoothPoints + 1)

      for(k in GapsToFill){
        BinPercMat[names(WhichStage)[j], WhichStage[[j]][k]:WhichStage[[j]][k+1]] <<- TRUE
      }

    })



    # Fill gaps with FALSE

    if(Structure == 'Circle'){

      SelGrop <- apply(BinPercMat, 1, function(x){
        range(which(!x))})

      for(j in 1:ncol(SelGrop)){
        if(any(is.infinite(SelGrop[,j]))){
          next
        } else {
          Diff <- SelGrop[1,j] + ncol(BinPercMat) + 1 - SelGrop[2,j] - 2

          if(Diff > 0 & Diff <= SmoothPoints){
            BinPercMat[colnames(SelGrop)[j],1:(SelGrop[1,j]-1)] <- FALSE
            BinPercMat[colnames(SelGrop)[j],(SelGrop[2,j]+1):ncol(BinPercMat)] <- FALSE
          }
        }
      }

    }

    WhichStage <- apply(BinPercMat, 1, function(x){which(!x)})

    WhichStage.Diff <- lapply(WhichStage, function(x){x[-1] - x[-length(x)]})

    lapply(as.list(1:length(WhichStage.Diff)), function(j){
      GapsToFill <- which(WhichStage.Diff[[j]] > 1 & WhichStage.Diff[[j]] <= SmoothPoints + 1)

      for(k in GapsToFill){
        BinPercMat[names(WhichStage)[j], WhichStage[[j]][k]:WhichStage[[j]][k+1]] <<- FALSE
      }

    })






    print(1*BinPercMat)


    if(ComputeOverlaps){

      # constructing stage overlaps

      OvefLapCat <- list()

      for(i in 1:nrow(BinPercMat)){

        if(i == nrow(BinPercMat)){
          if(Structure == 'Circle'){
            OvefLapCat[[i]] <- BinPercMat[i,] & BinPercMat[1,]
          }
        } else {
          OvefLapCat[[i]] <- BinPercMat[i,] & BinPercMat[i+1,]
        }

      }

      BinPercMatExt <- NULL

      for(i in 1:nrow(BinPercMat)){

        if(i == 1){
          if(Structure == 'Circle'){
            BinPercMatExt <- rbind(BinPercMatExt, BinPercMat[i,] & !OvefLapCat[[1]] & !OvefLapCat[[nrow(BinPercMat)]],
                                   OvefLapCat[[i]])
          } else {
            BinPercMatExt <- rbind(BinPercMatExt, BinPercMat[i,] & !OvefLapCat[[1]],
                                   OvefLapCat[[i]])
          }
          next()
        }

        if(i == nrow(BinPercMat)){
          if(Structure == 'Circle'){
            BinPercMatExt <- rbind(BinPercMatExt, BinPercMat[i,] & !OvefLapCat[[i-1]] & !OvefLapCat[[i]],
                                   OvefLapCat[[i]])
          } else {
            BinPercMatExt <- rbind(BinPercMatExt, BinPercMat[i,] & !OvefLapCat[[i-1]])
          }
          next()
        }

        BinPercMatExt <- rbind(BinPercMatExt, BinPercMat[i,] & !OvefLapCat[[i]] & !OvefLapCat[[i-1]],
                               OvefLapCat[[i]])

      }

      OverlapCat <- paste(levels(Categories)[-length(levels(Categories))],
                          levels(Categories)[-1], sep = "+")

      if(Structure == 'Circle'){
        OverlapCat <- c(OverlapCat,
                        paste(levels(Categories)[length(levels(Categories))],
                              levels(Categories)[1], sep = "+")
        )
      }

      if(Structure == 'Circle'){
        rownames(BinPercMatExt) <- as.vector(rbind(rownames(BinPercMat), OverlapCat))
      } else {
        rownames(BinPercMatExt) <- as.vector(rbind(rownames(BinPercMat), c(OverlapCat, NA)))[-nrow(BinPercMat)*2]
      }

    } else {
      BinPercMatExt <- BinPercMat
    }

    print(1*BinPercMat)

    AllCat <- rownames(BinPercMatExt)

    if(Structure == "Circle"){
      ExtPath <- ExtPath[-length(ExtPath)]
      BinPercMatExt <- BinPercMatExt[, -ncol(BinPercMatExt)]
    }

    LowStages <- 1
    if(ComputeOverlaps){
      LowStages <- 1:2
    }

    while(any(BinPercMatExt[LowStages,ncol(BinPercMatExt)]) &
          all(!BinPercMatExt[setdiff(1:nrow(BinPercMatExt), LowStages),ncol(BinPercMatExt)])){
      BinPercMatExt <- cbind(BinPercMatExt[, ncol(BinPercMatExt)], BinPercMatExt[, -ncol(BinPercMatExt)])
      ExtPath <- c(ExtPath[length(ExtPath)], ExtPath[-length(ExtPath)])
    }

    HighStages <- nrow(BinPercMatExt)
    if(ComputeOverlaps & Structure == 'Circle'){
      HighStages <- c(nrow(BinPercMatExt)-1, nrow(BinPercMatExt))
    }

    while(any(BinPercMatExt[HighStages,1]) &
          all(!BinPercMatExt[setdiff(1:nrow(BinPercMatExt), HighStages),1])){
      BinPercMatExt <- cbind(BinPercMatExt[, -1], BinPercMatExt[, 1])
      ExtPath <- c(ExtPath[-1], ExtPath[1])
    }

    if(Structure == 'Circle'){
      ExtPath <- c(ExtPath, ExtPath[1])
      BinPercMatExt <- cbind(BinPercMatExt,
                             rep(FALSE, nrow(BinPercMatExt)))
    }

  } else {
    BinPercMatExt <- BinPercMat
    AllCat <- rownames(BinPercMatExt)
  }

  if(!is.null(dim(BinPercMatExt))){
    Idxs <- apply(BinPercMatExt, 1, which)
  } else {
    Idxs <- list(Cat = which(BinPercMatExt))
  }

  print("Step III")
  print(Idxs)

  Bond <- lapply(Idxs[lapply(Idxs, length) > 1], range)

  if(!S4Vectors::isSorted(unlist(Bond))){
    warning("Stages are not sequential!")
  }

  print(1*BinPercMat)

  print(unlist(Bond))

  if(any(PCACenter != FALSE)){
    NodeOnGenes <- t(Nodes %*% t(RotatioMatrix)) + PCACenter
  } else {
    NodeOnGenes <- t(Nodes %*% t(RotatioMatrix))
  }


  OrderedPoints <- OrderOnPath(PrinGraph = PrinGraph, Path = as.numeric(ExtPath),
                               PointProjections = PointProjections)

  if(OrderOnCat){
    SelIdxs <- order(apply(NodeOnGenes, 1, var), decreasing = TRUE)[1:nGenes]
  } else {
    SelIdxs <- 1:2
  }


  for(Idx in SelIdxs){

    p <- ggplot2::ggplot(data = data.frame(x=cumsum(OrderedPoints$PathLen),
                                           y=NodeOnGenes[Idx,as.numeric(ExtPath)]),
                         mapping = ggplot2::aes(x = x, y = y, color="PC")) +
      ggplot2::labs(x = "Pseudotime", y="Gene expression",
                    title = paste(Title, rownames(NodeOnGenes)[Idx]), sep=' / ') +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

    RecCoord <- NULL

    if(OrderOnCat){

      for(i in 1:length(Idxs)){
        LowIds <- Idxs[[i]]-1
        LowIds[LowIds == 0] <- NA

        HiIds <- Idxs[[i]]+1
        HiIds[HiIds == length(ExtPath) + 1] <- NA

        LowCoord <- cumsum(OrderedPoints$PathLen)[LowIds]
        MidCoord <- cumsum(OrderedPoints$PathLen)[Idxs[[i]]]
        HighCoord <- cumsum(OrderedPoints$PathLen)[HiIds]

        RecCoord <- rbind(RecCoord, cbind(
          colMeans(rbind(LowCoord, MidCoord)),
          MidCoord,
          colMeans(rbind(MidCoord, HighCoord)),
          rep(names(Idxs)[i], length(MidCoord))
        )
        )

      }

      colnames(RecCoord) <- c("Min", "Med", "Max", "Stage")
      RecCoord <- data.frame(RecCoord)
      RecCoord$Min <- as.numeric(as.character(RecCoord$Min))
      RecCoord$Med <- as.numeric(as.character(RecCoord$Med))
      RecCoord$Max <- as.numeric(as.character(RecCoord$Max))


      RecCoord$Stage <- factor(as.character(RecCoord$Stage), levels = AllCat)

      RecCoord$Min[is.na(RecCoord$Min)] <- RecCoord$Med[is.na(RecCoord$Min)]
      RecCoord$Max[is.na(RecCoord$Max)] <- RecCoord$Med[is.na(RecCoord$Max)]

      p <- p + ggplot2::geom_rect(data = RecCoord, mapping = ggplot2::aes(fill=Stage, xmin=Min, xmax=Max),
                                  ymin = -Inf, ymax = Inf, inherit.aes = FALSE, alpha=.4) +
        ggplot2::geom_point(data = data.frame(x=OrderedPoints$PositionOnPath, y=ExpData[,Idx]),
                            mapping = ggplot2::aes(x=x, y=y, color="Data"), inherit.aes = FALSE, alpha=.5) +
        ggplot2::geom_point() + ggplot2::geom_line() + ggplot2::scale_color_manual(values = c("blue", "black"))

      print(p)

    } else {

      p <- p + ggplot2::geom_point(data = data.frame(x=OrderedPoints$PositionOnPath, y=ExpData[,Idx]),
                            mapping = ggplot2::aes(x=x, y=y, color="Data"), inherit.aes = FALSE, alpha=.5) +
        ggplot2::geom_point() + ggplot2::geom_line() + ggplot2::scale_color_manual(values = c("blue", "black"))

    }



  }

  CellsPT <- OrderedPoints$PositionOnPath
  names(CellsPT) <-  rownames(ExpData)

  return(list(Structure = Structure, ExtPath = ExtPath,
              NodesPT = cumsum(OrderedPoints$PathLen),
              NodesExp = NodeOnGenes[order(apply(NodeOnGenes, 1, var), decreasing = TRUE),as.numeric(ExtPath)],
              CellsPT = CellsPT, BinPercMatExt = BinPercMatExt, BinPercMat = BinPercMat,
              CellExp = t(ExpData[order(apply(ExpData, 1, var), decreasing = TRUE),]),
              RecCoord = RecCoord,
              StageOnNodes = Idxs
              )
         )

}












#' Title
#'
#' @param WorkStruct
#' @param Expression
#' @param Name
#' @param gName
#' @param SpanVal
#' @param CatOrder
#'
#' @return
#' @export
#'
#' @examples
PlotOnPseudotime <- function(WorkStruct,
                             Expression,
                             Name = '',
                             gName,
                             SpanVal=.3,
                             CatOrder = NULL,
                             legend.position = "bottom") {

  ReOrd.Sel <- match(colnames(WorkStruct$CellExp), names(WorkStruct$CellsPT))
  ReOrd.Sel.Mat <- match(colnames(Expression), names(WorkStruct$CellsPT))

  if(!is.null(CatOrder)){
    WorkStruct$RecCoord$Stage <- factor(as.character(WorkStruct$RecCoord$Stage), levels = CatOrder)
  }

  if(gName %in% rownames(WorkStruct$NodesExp)){

    SmoothData <- data.frame(x=c(WorkStruct$CellsPT[ReOrd.Sel],
                                 WorkStruct$CellsPT[ReOrd.Sel] + max(WorkStruct$NodesPT),
                                 WorkStruct$CellsPT[ReOrd.Sel] - max(WorkStruct$NodesPT)),
                             y=c(WorkStruct$CellExp[gName,],
                                 WorkStruct$CellExp[gName,],
                                 WorkStruct$CellExp[gName,]))

    p1 <- ggplot2::ggplot(data = data.frame(x=WorkStruct$NodesPT,
                                            y=WorkStruct$NodesExp[gName,]),
                          mapping = ggplot2::aes(x = x, y = y, color="PC")) +
      ggplot2::labs(x = "Pseudotime", y="Gene expression", title = paste(Name, "/", gName)) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = legend.position) +
      ggplot2::geom_rect(data = WorkStruct$RecCoord,
                         mapping = ggplot2::aes(fill=Stage, xmin=Min, xmax=Max),
                         ymin = -Inf, ymax = Inf, inherit.aes = FALSE, alpha=.4) +
      ggplot2::geom_smooth(data = SmoothData,
                           mapping = ggplot2::aes(x=x, y=y, color="Data"),
                           inherit.aes = FALSE, span = SpanVal, method = "loess") +
      ggplot2::geom_point(data = data.frame(x=WorkStruct$CellsPT[ReOrd.Sel],
                                            y=WorkStruct$CellExp[gName,]),
                          mapping = ggplot2::aes(x=x, y=y, color="Data"), inherit.aes = FALSE, alpha=.5) +
      ggplot2::geom_point() + ggplot2::geom_line() + ggplot2::scale_color_manual(values = c(Data = "blue", PC = "black")) +
      ggplot2::coord_cartesian(xlim = c(0, max(WorkStruct$NodesPT)))

  } else {

    if(gName %in% rownames(Expression)){

      SmoothData <- data.frame(x=c(WorkStruct$CellsPT[ReOrd.Sel.Mat],
                                   WorkStruct$CellsPT[ReOrd.Sel.Mat] + max(WorkStruct$NodesPT),
                                   WorkStruct$CellsPT[ReOrd.Sel.Mat] - max(WorkStruct$NodesPT)),
                               y=c(unlist(Expression[gName,]),
                                   unlist(Expression[gName,]),
                                   unlist(Expression[gName,])))

      p1 <- ggplot2::ggplot(data = data.frame(x=WorkStruct$CellsPT[ReOrd.Sel.Mat],
                                              y=unlist(Expression[gName,])),
                            mapping = ggplot2::aes(x=x, y=y, color="Data")) +
        ggplot2::labs(x = "Pseudotime", y="Gene expression", title = paste(Name, "/", gName)) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = legend.position) +
        ggplot2::geom_rect(data = WorkStruct$RecCoord,
                           mapping = ggplot2::aes(fill=Stage, xmin=Min, xmax=Max),
                           ymin = -Inf, ymax = Inf, inherit.aes = FALSE, alpha=.4) +
        ggplot2::geom_smooth(data = SmoothData,
                             mapping = ggplot2::aes(x=x, y=y, color="Data"),
                             inherit.aes = FALSE, span = SpanVal, method = "loess") +
        ggplot2::geom_point(alpha=.5) + ggplot2::scale_color_manual(values = c(Data = "blue", PC = "black")) +
        ggplot2::coord_cartesian(xlim = c(0, max(WorkStruct$NodesPT)))

    } else {

      p1 <- ggplot2::ggplot(data = data.frame(x=WorkStruct$CellsPT[ReOrd.Sel.Mat],
                                              y=unlist(Expression[1,])),
                            mapping = ggplot2::aes(x=x, y=y, color="Data")) +
        ggplot2::labs(x = "Pseudotime", y="Gene expression", title = paste(Name, "/", gName)) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = legend.position) +
        ggplot2::geom_rect(data = WorkStruct$RecCoord,
                           mapping = ggplot2::aes(fill=Stage, xmin=Min, xmax=Max),
                           ymin = -Inf, ymax = Inf, inherit.aes = FALSE, alpha=.4) +
        ggplot2::scale_color_manual(values = c(Data = "blue", PC = "black")) +
        ggplot2::coord_cartesian(xlim = c(0, max(WorkStruct$NodesPT)))

    }

  }

  return(p1)

}













