
################################################################################
#
# Function used to project cells on a circle or lasso ------------------------------------------------------
#
################################################################################


#' Title
#'
#' @param DataSet
#' @param GeneSet
#' @param OutThr
#' @param VarThr
#' @param nNodes
#' @param Log
#' @param Categories
#' @param Filter
#' @param GraphType
#' @param PlanVarLimit
#' @param PlanVarLimitIC
#' @param MinBranDiff
#' @param InitStructNodes
#' @param ForceLasso
#' @param DipPVThr
#' @param MinProlCells
#' @param PCAFilter
#' @param OutThrPCA
#' @param EstProlif
#' @param QuaThr
#' @param NonG0Cell
#' @param PCACenter
#' @param PCAProjCenter
#' @param PlotDebug
#' @param PlotIntermediate
#'
#' @return
#' @export
#'
#' @examples
#'
ProjectAndCompute <- function(DataSet,
                              GeneSet = NULL,
                              VarThr,
                              nNodes = 30,
                              Log = TRUE,
                              Categories = NULL,
                              Filter = TRUE,
                              OutThr = 5,
                              PCAFilter = FALSE,
                              OutThrPCA = 5,
                              GraphType = 'Circle',
                              PlanVarLimit = .9,
                              PlanVarLimitIC = NULL,
                              MinBranDiff = 2,
                              InitStructNodes = 15,
                              ForceLasso = FALSE,
                              EstProlif = "MeanPerc",
                              QuaThr = .5,
                              NonG0Cell = NULL,
                              DipPVThr = 1e-3,
                              MinProlCells = 20,
                              PCACenter = TRUE,
                              PCAProjCenter = FALSE,
                              PlotDebug = FALSE,
                              PlotIntermediate = FALSE){

  if(is.null(PlanVarLimitIC)){
    PlanVarLimitIC <- PlanVarLimit + .5*(1-PlanVarLimit)
  }

  DataMat <- t(DataSet)

  if(is.null(Categories)){
    Categories <- factor(rep("NoG", nrow(DataMat)))
  }

  if(!is.null(GeneSet)){
    DataMat <- DataMat[, colnames(DataMat) %in% GeneSet]
  }

  DataMat <- DataMat[, apply(DataMat > 0, 2, sum) > 3]


  if(PlotDebug){
    hist(apply(DataMat > 0, 1, sum), main = "Genes per cell", xlab = "Genes count",
         freq = TRUE, ylab = "Number of cells")

    hist(apply(DataMat, 1, sum), main = "Reads per cell", xlab = "Reads count",
         freq = TRUE, ylab = "Number of cells")

    hist(apply(DataMat>0, 2, sum), main = "Transcripts per cell", xlab = "Reads count",
         freq = TRUE, ylab = "Number of genes identified")
  }


  if(Filter){
    OutExpr <- scater::isOutlier(rowSums(DataMat>0), nmads = OutThr)
    OutCount <- scater::isOutlier(rowSums(DataMat), nmads = OutThr)
  } else {
    OutExpr <- rep(FALSE, nrow(DataMat))
    OutCount <- rep(FALSE, nrow(DataMat))
  }

  if(Log){
    DataMat <- log10(DataMat[!(OutExpr & OutCount),] + 1)
  } else {
    DataMat <- DataMat[!(OutExpr & OutCount),]
  }

  Categories <- Categories[!(OutExpr & OutCount)]

  if(PCAFilter){
    PCAData <- prcomp(DataMat, retx = TRUE, center = PCACenter, scale. = FALSE)
    Centroid <- colMeans(PCAData$x)

    Dists <- as.matrix(dist(rbind(Centroid, PCAData$x)))
    DistFromCent <- Dists[1,-1]

    PCAFil <- scater::isOutlier(DistFromCent, nmads = OutThrPCA)
  } else {
    PCAFil <- rep(FALSE, nrow(DataMat))
  }

  DataMat <- DataMat[!PCAFil, ]
  Categories <- Categories[!PCAFil]

  Categories <- Categories[rowSums(DataMat)>0]
  DataMat <- DataMat[rowSums(DataMat)>0, colSums(DataMat)>0]

  print("Estimating most likely proliferative cells")
  RankedData <- apply(DataMat, 1, rank)


  if(is.null(NonG0Cell)){

    if(EstProlif == "Quantile"){
      print("Using quantile separation")
      NonG0Cell <-  rownames(DataMat)[apply(DataMat, 1, median) > quantile(DataMat, QuaThr)]
    }

    if(EstProlif == "PercQuant"){
      print("Using quantile ordering")
      NonG0Cell <-  rownames(DataMat)[order(apply(DataMat, 1, quantile, QuaThr), decreasing = TRUE)[1:round(nrow(DataMat)/10)]]
    }

    if(EstProlif == "KmeansPerc"){
      print("Using kmeans on quantile data")
      Cellvect <- apply(DataMat/rowSums(DataMat), 1, quantile, QuaThr)

      KM <- kmeans(x = Cellvect, centers = range(Cellvect))
      if(PlotDebug){
        boxplot(apply(DataMat/rowSums(DataMat), 1, median) ~ KM$cluster)
      }

      if(KM$centers[1] < KM$centers[2]){
        NonG0Cell <-  rownames(DataMat)[KM$cluster == 2]
      } else {
        NonG0Cell <-  rownames(DataMat)[KM$cluster == 1]
      }

    }

    if(EstProlif == "MeanPerc"){
      print("Using mean of quantiles")
      Cellvect <- apply(DataMat/rowSums(DataMat), 1, quantile, QuaThr)
      NonG0Cell <- rownames(DataMat)[Cellvect > mean(Cellvect)]
      boxplot(Cellvect)
      abline(h=mean(Cellvect))
    }

    if(EstProlif == "AdaptiveMeanPerc"){
      print("Using adaptive mean of quantiles")
      Cellvect <- apply(DataMat/rowSums(DataMat), 1, quantile, QuaThr)
      AdjFact = 0
      while((sum(Cellvect != 0) < MinProlCells) & (QuaThr + AdjFact < 1)){
        AdjFact <- AdjFact + .01
        Cellvect <- apply(DataMat/rowSums(DataMat), 1, quantile, QuaThr + AdjFact)
      }
      boxplot(Cellvect, main=paste("Q =", QuaThr + AdjFact))
      abline(h=mean(Cellvect))
      NonG0Cell <- rownames(DataMat)[Cellvect > mean(Cellvect)]
    }

  } else {
    NonG0Cell <- intersect(NonG0Cell, rownames(DataMat))
  }

  print(paste(length(NonG0Cell), "strongly proliferative cells inferred"))

  G0Cell <-  setdiff(rownames(DataMat), NonG0Cell)

  if(!is.null(NonG0Cell)){
    TB <- rbind(table(Categories[rownames(DataMat) %in% NonG0Cell]), table(Categories[rownames(DataMat) %in% G0Cell]))
    barplot(TB/rowSums(TB), ylab = "Percentage of cells", xlab = "Category", beside = TRUE,
            legend.text = c("Strongly Proliferative", "Not strongly proliferative"))
    barplot(TB, ylab = "Number of cells", xlab = "Category", beside = TRUE,
            legend.text = c("Strongly Proliferative", "Not strongly proliferative"))

  }



  if(length(NonG0Cell)>MinProlCells){

    print("Using strongly proliferative cells for initial fitting")
    UseTree <- TRUE

  } else {

    print("Unable to find a sufficent number of strongly proliferative cells")

    NonG0Cell <- rownames(DataMat)
    UseTree <- FALSE

  }


  nDims <- min(dim(DataMat))

  print("Transforming data using PCA")

  PCAData <- prcomp(DataMat, retx = TRUE, center = PCACenter, scale.=FALSE)
  ExpVar <- PCAData$sdev^2/sum(PCAData$sdev^2)

  if(VarThr<1){
    nDims <- max(min(which(cumsum(ExpVar) > VarThr)), 2)
  }

  Data <- PCAData$x[, 1:nDims]


  if(GraphType == 'Lasso') {
    FitData <- FitLasso(Data, NonG0Cell, InitStructNodes, PlanVarLimitIC, PlanVarLimit, UseTree, ForceLasso, nNodes, MinBranDiff)
  }


  if(GraphType == 'Circle') {
    FitData <- FitCircle(Data, NonG0Cell, InitStructNodes, PlanVarLimitIC, PlanVarLimit, nNodes)
  }


  if(GraphType == 'Line') {
    FitData <- FitLine(Data, NonG0Cell, InitStructNodes, PlanVarLimitIC, PlanVarLimit, nNodes)
  }


  CombData <- FitData[[length(FitData)]]
  CombData$Report <- FitData[[1]]$Report

  for(i in 2:length(FitData)){
    CombData$Report <- rbind(CombData$Report, FitData[[i]]$Report)
  }

  Net <- list()
  TaxonList <- list()
  InfoData <- list()
  ProjPoints <- list()
  PCAPrGraph <- list()

  for(i in 1:length(FitData)){

    print(paste("Constructing accessory structures - round", i))

    Net[[i]] <- ConstructGraph(Results = FitData[[i]], DirectionMat = NULL, Thr = 0.05)
    TaxonList[[i]] <- getTaxonMap(Results = FitData[[i]], Data = Data)

    ProjPoints[[i]] <- projectPoints(Results = FitData[[i]], Data = Data, TaxonList = TaxonList[[i]],
                                     UseR = TRUE,
                                     method = 'PCALin', Dims = nDims, Debug = FALSE)

    if(i != length(FitData)){
      InfoData[[i]] <- plotPieNet(Results = FitData[[i]], Data = Data, NodeSizeMult = 4,
                                  Categories = Categories, PlotNet = FALSE,
                                  Graph = Net[[i]], TaxonList = TaxonList[[i]], LayOut =, Main = "Pincipal graph")
    } else {
      InfoData[[i]] <- plotPieNet(Results = FitData[[i]], Data = Data, NodeSizeMult = 4,
                                  Categories = Categories, PlotNet = TRUE,
                                  Graph = Net[[i]], TaxonList = TaxonList[[i]], LayOut =, Main = "Pincipal graph")
    }


    PCAPrGraph[[i]] <-  prcomp(FitData[[i]]$Nodes, retx = TRUE, center = FALSE, scale. = FALSE)

    # if(FitData[[i]]$Method == "CircleConfiguration"){
    #   RotatedData <- cbind(Data %*% PCAPrGraph[[i]]$rotation[,1:2], 1:nrow(Data) %in% NonG0Cell,
    #                        as.character(Categories))
    # } else {
    #   RotatedData <- cbind(Data %*% PCAPrGraph[[i]]$rotation[,1:2], rep(TRUE, nrow(Data)),
    #                        as.character(Categories))
    # }
    #
    # colnames(RotatedData) <- c("PC1", "PC2", "NG0", "Cat")
    #
    # RotatedData.DF <- data.frame(RotatedData)
    # RotatedData.DF$PC1 <- as.numeric(as.character(RotatedData.DF$PC1))
    # RotatedData.DF$PC2 <- as.numeric(as.character(RotatedData.DF$PC2))
    # RotatedData.DF$NG0 <- factor(RotatedData.DF$NG0, levels = c("TRUE", "FALSE"))
    #
    #
    # p <- ggplot(data.frame(RotatedData.DF), aes(x=PC1, y=PC2, alpha=NG0, colour=Cat)) + geom_point() +
    #   geom_point(data = data.frame(PCAPrGraph[[i]]$x[,1:2]), mapping = aes(x=PC1, y=PC2),
    #              inherit.aes = FALSE) +
    #   labs(title = paste("Round", i)) + scale_alpha_discrete("Fitted", range = c(1, .1))
    #
    # for(j in 1:nrow(FitData[[i]]$Edges)){
    #   p <- p + geom_path(data = data.frame(PCAPrGraph[[i]]$x[FitData[[i]]$Edges[j,],1:2]),
    #                      mapping = aes(x = PC1, y = PC2), inherit.aes = FALSE)
    # }
    #
    # print(p)

    if(PlotIntermediate & i != length(FitData)){
      if(FitData[[i]]$Method == "CircleConfiguration" & GraphType == 'Lasso'){
        ProjectOnPrincipalGraph(Nodes = FitData[[i]]$Nodes, Edges = FitData[[i]]$Edges, Points = Data,
                                UsedPoints = which(rownames(Data) %in% NonG0Cell), Categories = Categories, Title= paste("Round", i),
                                PCACenter = PCAProjCenter)
      } else {
        ProjectOnPrincipalGraph(Nodes = FitData[[i]]$Nodes, Edges = FitData[[i]]$Edges, Points = Data,
                                UsedPoints = NULL, Categories = Categories, Title= paste("Round", i),
                                PCACenter = PCAProjCenter)
      }
    }

    if(i == length(FitData)){
      ProjectOnPrincipalGraph(Nodes = FitData[[i]]$Nodes, Edges = FitData[[i]]$Edges, Points = Data,
                              UsedPoints = NULL, Categories = Categories, Title= paste("Round", i),
                              PCACenter = PCAProjCenter)
    }

  }


  # legend("center", legend = levels(Categories), fill = InfoData$ColInfo[1:3], cex = 4)

  return(list(Data = Data, FiltExp = DataMat, Categories = Categories, PrinGraph = CombData,
              IntGrahs = FitData, Net = Net, TaxonList = TaxonList, PCAData = PCAData,
              nDims = nDims, ProjPoints = ProjPoints, PCAPrGraph = PCAPrGraph, NonG0Cell = NonG0Cell))
}












################################################################################
#
# Function used to construct and optimize the cycle ------------------------------------------------------
#
################################################################################



#' Construct and optimize cycle
#'
#' @param DataSet
#' @param StartSet
#' @param VarThr
#' @param nNodes
#' @param Log
#' @param Categories
#' @param Filter
#' @param OutThr
#' @param PCAFilter
#' @param OutThrPCA
#' @param GraphType
#' @param PlanVarLimit
#' @param PlanVarLimitIC
#' @param MinBranDiff
#' @param InitStructNodes
#' @param ForceLasso
#' @param EstProlif
#' @param QuaThr
#' @param NonG0Cell
#' @param DipPVThr
#' @param MinProlCells
#' @param PCACenter
#' @param PCAProjCenter
#' @param PlotDebug
#' @param PlotIntermediate
#' @param AddGenePerc
#' @param SelThr1
#' @param SelThr2
#' @param MadsThr
#'
#' @return
#' @export
#'
#' @examples
SelectGenesOnGraph <- function(DataSet,
                               StartSet,
                               VarThr = .99,
                               nNodes = 40,
                               Log = TRUE,
                               Categories = NULL,
                               Filter = TRUE,
                               OutThr = 5,
                               PCAFilter = TRUE,
                               OutThrPCA = 3,
                               GraphType = "Circle",
                               PlanVarLimit = .85,
                               PlanVarLimitIC = .9,
                               MinBranDiff = 3,
                               InitStructNodes = 20,
                               ForceLasso = FALSE,
                               EstProlif = "MeanPerc",
                               QuaThr = .5,
                               NonG0Cell = NULL,
                               DipPVThr = 1e-4,
                               MinProlCells = 50,
                               PCACenter = TRUE,
                               PCAProjCenter = TRUE,
                               PlotDebug = FALSE,
                               PlotIntermediate = FALSE,
                               AddGenePerc = 5,
                               SelThr1 = .95,
                               SelThr2 = .99,
                               MadsThr =  1) {

  # Construct the initial Principal graph

  Steps <- list()
  UsedGenes <- list()

  UsedGenes[[1]] <- intersect(StartSet, rownames(DataSet))

  Steps[[1]] <- ProjectAndCompute(DataSet = DataSet,
                                  GeneSet = UsedGenes[[1]],
                                  VarThr = VarThr,
                                  nNodes = nNodes,
                                  Log = Log,
                                  Categories = Categories,
                                  Filter = Filter,
                                  OutThr = OutThr,
                                  PCAFilter = PCAFilter,
                                  OutThrPCA = OutThrPCA,
                                  GraphType = GraphType,
                                  PlanVarLimit = PlanVarLimit,
                                  PlanVarLimitIC = PlanVarLimitIC,
                                  MinBranDiff = MinBranDiff,
                                  InitStructNodes = InitStructNodes,
                                  ForceLasso = ForceLasso,
                                  EstProlif = EstProlif,
                                  QuaThr = QuaThr,
                                  NonG0Cell = NonG0Cell,
                                  DipPVThr = DipPVThr,
                                  MinProlCells = MinProlCells,
                                  PCACenter = PCACenter,
                                  PCAProjCenter = PCAProjCenter,
                                  PlotDebug = PlotDebug,
                                  PlotIntermediate = PlotIntermediate)


  # Select the genes that are closer to the original curve

  CONVERGED <- FALSE
  i = 2

  DoStuff <- function(TaxVect){

    # SelCell <- intersect(Steps[[1]]$NonG0Cell, names(TaxVect))
    SelCell <- names(TaxVect)

    GeneExprMat.DF.Split <- split(GeneExprMat.DF[SelCell, ], TaxVect[SelCell])

    GeneExprMat.DF.Split.Mean <- lapply(GeneExprMat.DF.Split, colMeans)

    GeneExprMat.DF.MeanRemoved <-
      lapply(as.list(1:length(GeneExprMat.DF.Split)), function(i){
        apply(GeneExprMat.DF.Split[[i]], 2, sd)/GeneExprMat.DF.Split.Mean[[i]]
      })

    GeneExprMat.DF.MeanRemoved.All <- NULL
    for(i in 1:length(GeneExprMat.DF.MeanRemoved)){
      GeneExprMat.DF.MeanRemoved.All <- rbind(GeneExprMat.DF.MeanRemoved.All, GeneExprMat.DF.MeanRemoved[[i]])
    }

    return(apply(abs(GeneExprMat.DF.MeanRemoved.All), 2, mean, na.rm = TRUE))

  }



  while(!CONVERGED){

    TaxList <- Steps[[i-1]]$TaxonList[[length(Steps[[i-1]]$TaxonList)]]
    TaxVect <- rep(NA, nrow(Steps[[i-1]]$Data))
    names(TaxVect) <- rownames(Steps[[i-1]]$Data)

    lapply(1:length(TaxList), function(j){
      TaxVect[ TaxList[[j]] ] <<- j
    })

    GeneExprMat.DF <- data.frame(t(DataSet[UsedGenes[[i-1]], ]), row.names = colnames(DataSet))
    colnames(GeneExprMat.DF) <- rownames(DataSet[UsedGenes[[i-1]], ])

    Base <- DoStuff(TaxVect)
    Base <- Base[is.finite(Base)]

    if(i == 2){
      Thr <- median(Base) + MadsThr*mad(Base)
    }

    UsedGenes[[i]] <- names(which(Base < Thr))

    length(setdiff(UsedGenes[[length(UsedGenes) - 1]], UsedGenes[[length(UsedGenes)]]))
    length(UsedGenes[[length(UsedGenes)]])

    Steps[[i]] <- ProjectAndCompute(DataSet = DataSet,
                                    GeneSet = UsedGenes[[i]],
                                    VarThr = VarThr,
                                    nNodes = nNodes,
                                    Log = Log,
                                    Categories = Categories,
                                    Filter = Filter,
                                    OutThr = OutThr,
                                    PCAFilter = PCAFilter,
                                    OutThrPCA = OutThrPCA,
                                    GraphType = GraphType,
                                    PlanVarLimit = PlanVarLimit,
                                    PlanVarLimitIC = PlanVarLimitIC,
                                    MinBranDiff = MinBranDiff,
                                    InitStructNodes = InitStructNodes,
                                    ForceLasso = ForceLasso,
                                    EstProlif = EstProlif,
                                    QuaThr = QuaThr,
                                    NonG0Cell = Steps[[1]]$NonG0Cell,
                                    DipPVThr = DipPVThr,
                                    MinProlCells = MinProlCells,
                                    PCACenter = PCACenter,
                                    PCAProjCenter = PCAProjCenter,
                                    PlotDebug = PlotDebug,
                                    PlotIntermediate = PlotIntermediate)

    i = i + 1

    if(length(UsedGenes[[i-1]])/length(UsedGenes[[i - 2]]) > SelThr1){
      CONVERGED <- TRUE
    }

  }



  StageIEnd <- i

  # Extend the geneset by including other genes that are closer to the original circle

  CONVERGED2 <- FALSE

  while(!CONVERGED2){

    TaxList <- Steps[[i-1]]$TaxonList[[length(Steps[[i-1]]$TaxonList)]]
    TaxVect <- rep(NA, nrow(Steps[[i-1]]$Data))
    names(TaxVect) <- rownames(Steps[[i-1]]$Data)

    lapply(1:length(TaxList), function(j){
      TaxVect[ TaxList[[j]] ] <<- j
    })

    SelGenes <- rownames(DataSet)

    GeneExprMat.DF <- data.frame(t(DataSet[SelGenes, ]))
    colnames(GeneExprMat.DF) <- rownames(DataSet[SelGenes, ])

    DoStuff <- function(TaxVect){

      # SelCell <- intersect(Steps[[1]]$NonG0Cell, names(TaxVect))
      SelCell <- names(TaxVect)

      GeneExprMat.DF.Split <- split(GeneExprMat.DF[SelCell, ], TaxVect[SelCell])

      GeneExprMat.DF.Split.Mean <- lapply(GeneExprMat.DF.Split, colMeans)

      GeneExprMat.DF.MeanRemoved <-
        lapply(as.list(1:length(GeneExprMat.DF.Split)), function(i){
          apply(GeneExprMat.DF.Split[[i]], 2, sd)/GeneExprMat.DF.Split.Mean[[i]]
        })

      GeneExprMat.DF.MeanRemoved.All <- NULL
      for(i in 1:length(GeneExprMat.DF.MeanRemoved)){
        GeneExprMat.DF.MeanRemoved.All <- rbind(GeneExprMat.DF.MeanRemoved.All, GeneExprMat.DF.MeanRemoved[[i]])
      }

      return(apply(abs(GeneExprMat.DF.MeanRemoved.All), 2, mean, na.rm = TRUE))

    }

    Base <- DoStuff(TaxVect)
    Base <- Base[is.finite(Base)]

    AddGenes <- Base[order(Base)[1:round(AddGenePerc*length(UsedGenes[[i-1]])/100)]]
    AddGenes <- AddGenes[AddGenes < Thr]

    UsedGenes[[i]] <- union(names(AddGenes), UsedGenes[[i-1]])

    length(setdiff(UsedGenes[[i - 1]], UsedGenes[[i]]))

    Steps[[i]] <- ProjectAndCompute(DataSet = DataSet,
                                    GeneSet = UsedGenes[[i]],
                                    VarThr = VarThr,
                                    nNodes = nNodes,
                                    Log = Log,
                                    Categories = Categories,
                                    Filter = Filter,
                                    OutThr = OutThr,
                                    PCAFilter = PCAFilter,
                                    OutThrPCA = OutThrPCA,
                                    GraphType = GraphType,
                                    PlanVarLimit = PlanVarLimit,
                                    PlanVarLimitIC = PlanVarLimitIC,
                                    MinBranDiff = MinBranDiff,
                                    InitStructNodes = InitStructNodes,
                                    ForceLasso = ForceLasso,
                                    EstProlif = EstProlif,
                                    QuaThr = QuaThr,
                                    NonG0Cell = Steps[[1]]$NonG0Cell,
                                    DipPVThr = DipPVThr,
                                    MinProlCells = MinProlCells,
                                    PCACenter = PCACenter,
                                    PCAProjCenter = PCAProjCenter,
                                    PlotDebug = PlotDebug,
                                    PlotIntermediate = PlotIntermediate)

    i = i + 1

    if(length(UsedGenes[[i-2]])/length(UsedGenes[[i-1]]) > SelThr2){
      CONVERGED2 <- TRUE
    }

  }

  plot(unlist(lapply(UsedGenes, length)))
  abline(v=StageIEnd-.5)
  abline(v=1.5)

  return(list(Genes = UsedGenes, PGStructs = Steps))

}








