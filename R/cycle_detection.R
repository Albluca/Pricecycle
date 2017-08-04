
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
#'
#' @return
#' @export
#'
#' @examples
#'
ProjectAndCompute <- function(DataSet, GeneSet = NULL, VarThr, nNodes, Log = TRUE, Categories = NULL,
                              Filter = TRUE, OutThr = 5, PCAFilter = FALSE, OutThrPCA = 5,
                              GraphType = 'Lasso', PlanVarLimit = .9,
                              PlanVarLimitIC = NULL, MinBranDiff = 2, InitStructNodes = 15,
                              ForceLasso = FALSE, EstProlif = "MeanPerc", QuaThr = .5,
                              NonG0Cell = NULL, DipPVThr = 1e-3, MinProlCells = 20, PCACenter = TRUE, PCAProjCenter = FALSE,
                              PlotDebug = FALSE, PlotIntermediate = FALSE){

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

  print("Transforming Data")

  PCAData <- prcomp(DataMat, retx = TRUE, center = PCACenter, scale.=FALSE)
  ExpVar <- PCAData$sdev^2/sum(PCAData$sdev^2)

  if(VarThr<1){
    nDims <- max(min(which(cumsum(ExpVar) > VarThr)), 2)
  } else {
    nDims <- max(dim(PCAData$x))
  }

  Data <- PCAData$x[, 1:nDims]


  if(GraphType == 'Lasso') {

    print("Lasso fitting")
    print("Fitting initial circle")

    # Step I - Construct the base circle

    BasicCircData <- computeElasticPrincipalGraph(Data = Data[rownames(Data) %in% NonG0Cell,], NumNodes = 4,
                                                  Method = 'CircleConfiguration', NodeStep = 1)

    PCAStruct <- prcomp(BasicCircData[[length(BasicCircData)]]$Nodes, center = TRUE, scale. = FALSE, retx = TRUE)
    PlanPerc <- sum(PCAStruct$sdev[1:2]^2)/sum(PCAStruct$sdev^2)

    UsedNodes <- nrow(BasicCircData[[length(BasicCircData)]]$Nodes) - 1

    while(UsedNodes < InitStructNodes & PlanPerc > PlanVarLimitIC){

      print("Expanding initial circle")

      # Contiune to add node untill the circle remains planar

      BasicCircData <- append(BasicCircData, computeElasticPrincipalGraph(Data = Data[rownames(Data) %in% NonG0Cell,],
                                                                          NumNodes = UsedNodes + 1,
                                                                          Method = 'CircleConfiguration',
                                                                          NodesPositions = BasicCircData[[length(BasicCircData)]]$Nodes,
                                                                          Edges = BasicCircData[[length(BasicCircData)]]$Edges))

      UsedNodes <- nrow(BasicCircData[[length(BasicCircData)]]$Nodes)

      PCAStruct <- prcomp(BasicCircData[[length(BasicCircData)]]$Nodes, center = TRUE, scale. = FALSE, retx = TRUE)
      PlanPerc <- sum(PCAStruct$sdev[1:2]^2)/sum(PCAStruct$sdev^2)

    }

    if(UseTree | ForceLasso){

      print("Branching initial circle")

      # Step II - Construct the initial tree

      BasicCircData <- append(BasicCircData, computeElasticPrincipalGraph(Data = Data,
                                                                          NumNodes = UsedNodes + 1,
                                                                          Method = 'DefaultPrincipalTreeConfiguration',
                                                                          NodesPositions = BasicCircData[[length(BasicCircData)]]$Nodes,
                                                                          Edges = BasicCircData[[length(BasicCircData)]]$Edges))


      PCAStruct <- prcomp(BasicCircData[[length(BasicCircData)]]$Nodes, center = TRUE, scale. = FALSE, retx = TRUE)
      PlanPerc <- sum(PCAStruct$sdev[1:2]^2)/sum(PCAStruct$sdev^2)

      UsedNodes <- nrow(BasicCircData[[length(BasicCircData)]]$Nodes) - 1

      while(UsedNodes < nNodes & PlanPerc > PlanVarLimit){

        print("Keep Branching")

        # Step IIa - keep constructing trees

        BasicCircData <- append(BasicCircData, computeElasticPrincipalGraph(Data = Data,
                                                                            NumNodes = UsedNodes + 1,
                                                                            Method = 'DefaultPrincipalTreeConfiguration',
                                                                            NodesPositions = BasicCircData[[length(BasicCircData)]]$Nodes,
                                                                            Edges = BasicCircData[[length(BasicCircData)]]$Edges))


        UsedNodes <- nrow(BasicCircData[[length(BasicCircData)]]$Nodes)

        PCAStruct <- prcomp(BasicCircData[[length(BasicCircData)]]$Nodes, center = TRUE, scale. = FALSE, retx = TRUE)
        PlanPerc <- sum(PCAStruct$sdev[1:2]^2)/sum(PCAStruct$sdev^2)

        Net <- ConstructGraph(Results = BasicCircData[[length(BasicCircData)]], DirectionMat = NULL, Thr = NULL)
        if(max(igraph::degree(Net))>2){
          break()
        }
      }

    }

    # Step III - Using curves

    while(UsedNodes < nNodes & PlanPerc > PlanVarLimit){

      print("Extending circle and branches")

      BasicCircData <- append(BasicCircData, computeElasticPrincipalGraph(Data = Data,
                                                                          NumNodes = nrow(BasicCircData[[length(BasicCircData)]]$Nodes)+1,
                                                                          Method = 'CurveConfiguration',
                                                                          NodesPositions = BasicCircData[[length(BasicCircData)]]$Nodes,
                                                                          Edges = BasicCircData[[length(BasicCircData)]]$Edges))

      Net <- ConstructGraph(Results = BasicCircData[[length(BasicCircData)]], DirectionMat = NULL, Thr = NULL)

      PCAStruct <- prcomp(BasicCircData[[length(BasicCircData)]]$Nodes, center = TRUE, scale. = FALSE, retx = TRUE)
      PlanPerc <- sum(PCAStruct$sdev[1:2]^2)/sum(PCAStruct$sdev^2)

      UsedNodes <- nrow(BasicCircData[[length(BasicCircData)]]$Nodes)

      if(sum(igraph::degree(Net)>2)>1){

        print("Multiple branches detected ... trying to select one")

        # There are two branhces, this is no good ...
        # Look for the biggest circle

        CircSize <- 3

        while(CircSize < igraph::vcount(Net)){

          Template <- igraph::graph.ring(n = CircSize, directed = FALSE, mutual = FALSE, circular = TRUE)

          CircSize <- CircSize + 1

          CircleInData <- igraph::graph.get.subisomorphisms.vf2(Net, Template)

          if(length(CircleInData)>0){
            # We found the cicle

            EndPoint <- igraph::V(Net)[igraph::degree(Net) == 1]
            Branches <- igraph::V(Net)[igraph::degree(Net)>2]
            DistMat <- igraph::distances(Net, EndPoint, Branches)

            BranchesLen <- apply(DistMat, 1, min)

            VertToRemove <- NULL

            if(max(BranchesLen[-which.max(BranchesLen)] - max(BranchesLen)) <= -MinBranDiff){

              print("Dominating brach found ... pruning the shortest one")

              # There is a "dominating"" branch
              ToKeep <- names(which.max(BranchesLen))
              for(i in 1:nrow(DistMat)){

                Source <- rownames(DistMat)[i]
                if(ToKeep == Source){
                  next()
                }

                Target <- names(which.min(DistMat[i,]))
                VertToRemove <- c(VertToRemove, names(igraph::get.shortest.paths(Net, from = Target, to = Source)$vpath[[1]][-1]))
              }

              PrunedStruct <- BasicCircData[[length(BasicCircData)]]

              NodesToRem <- as.integer(unlist(lapply(strsplit(VertToRemove, "V_"), "[[", 2)))
              PrunedStruct$Nodes <- PrunedStruct$Nodes[-NodesToRem, ]
              PrunedStruct$Edges <-
                PrunedStruct$Edges[!(PrunedStruct$Edges[,1] %in% NodesToRem | PrunedStruct$Edges[,2] %in% NodesToRem), ]

              for(VerVal in 1:max(PrunedStruct$Edges)){
                if(any(PrunedStruct$Edges == VerVal)){
                  next()
                } else {
                  PrunedStruct$Edges[PrunedStruct$Edges>VerVal] <-
                    PrunedStruct$Edges[PrunedStruct$Edges>VerVal] - 1
                }

              }

              PrunedStruct$Method <- "ManualPruning"

              BasicCircData[[length(BasicCircData) + 1]] <- PrunedStruct

              PCAStruct <- prcomp(BasicCircData[[length(BasicCircData)]]$Nodes, center = TRUE, scale. = FALSE, retx = TRUE)
              PlanPerc <- sum(PCAStruct$sdev[1:2]^2)/sum(PCAStruct$sdev^2)

              UsedNodes <- nrow(BasicCircData[[length(BasicCircData)]]$Nodes)

              break()

            }


          }


        }

      }

    }

    FitData <- BasicCircData
  }







  if(GraphType == 'Circle') {

    print("Circle fitting")
    print("Fitting initial circle")

    # Step I - Construct the base circle

    BasicCircData <- computeElasticPrincipalGraph(Data = Data[rownames(Data) %in% NonG0Cell,], NumNodes = 4,
                                                  Method = 'CircleConfiguration', NodeStep = 1)

    PCAStruct <- prcomp(BasicCircData[[length(BasicCircData)]]$Nodes, center = TRUE, scale. = FALSE, retx = TRUE)
    PlanPerc <- sum(PCAStruct$sdev[1:2]^2)/sum(PCAStruct$sdev^2)

    UsedNodes <- nrow(BasicCircData[[length(BasicCircData)]]$Nodes) - 1

    while(UsedNodes < InitStructNodes & PlanPerc > PlanVarLimitIC){

      print("Expanding initial circle")

      # Contiune to add node untill the circle remains planar

      BasicCircData <- append(BasicCircData, computeElasticPrincipalGraph(Data = Data[rownames(Data) %in% NonG0Cell,],
                                                                          NumNodes = UsedNodes + 1,
                                                                          Method = 'CircleConfiguration',
                                                                          NodesPositions = BasicCircData[[length(BasicCircData)]]$Nodes,
                                                                          Edges = BasicCircData[[length(BasicCircData)]]$Edges))


      UsedNodes <- nrow(BasicCircData[[length(BasicCircData)]]$Nodes)

      PCAStruct <- prcomp(BasicCircData[[length(BasicCircData)]]$Nodes, center = TRUE, scale. = FALSE, retx = TRUE)
      PlanPerc <- sum(PCAStruct$sdev[1:2]^2)/sum(PCAStruct$sdev^2)

      print(paste("Initial circle: Nodes = ", UsedNodes, "PercPlan = ", PlanPerc))

    }

    # Step II - Using curves

    while(UsedNodes < nNodes & PlanPerc > PlanVarLimit){

      print("Extending circle")

      BasicCircData <- append(BasicCircData, computeElasticPrincipalGraph(Data = Data,
                                                                          NumNodes = nrow(BasicCircData[[length(BasicCircData)]]$Nodes)+1,
                                                                          Method = 'CurveConfiguration',
                                                                          NodesPositions = BasicCircData[[length(BasicCircData)]]$Nodes,
                                                                          Edges = BasicCircData[[length(BasicCircData)]]$Edges))


      Net <- ConstructGraph(Results = BasicCircData[[length(BasicCircData)]], DirectionMat = NULL, Thr = NULL)

      PCAStruct <- prcomp(BasicCircData[[length(BasicCircData)]]$Nodes, center = TRUE, scale. = FALSE, retx = TRUE)
      PlanPerc <- sum(PCAStruct$sdev[1:2]^2)/sum(PCAStruct$sdev^2)

      UsedNodes <- nrow(BasicCircData[[length(BasicCircData)]]$Nodes)

      print(paste("Final circle: Nodes = ", UsedNodes, "PercPlan = ", PlanPerc))

    }

    FitData <- BasicCircData
  }








  if(GraphType == 'Line') {

    print("Line fitting")
    print("Fitting initial line")

    # Step I - Construct the base circle

    BasicLineData <- computeElasticPrincipalGraph(Data = Data[rownames(Data) %in% NonG0Cell,], NumNodes = 4,
                                                  Method = 'CurveConfiguration', NodeStep = 1)

    PCAStruct <- prcomp(BasicLineData[[length(BasicLineData)]]$Nodes, center = TRUE, scale. = FALSE, retx = TRUE)
    PlanPerc <- sum(PCAStruct$sdev[1:2]^2)/sum(PCAStruct$sdev^2)

    UsedNodes <- nrow(BasicLineData[[length(BasicLineData)]]$Nodes) - 1

    while(UsedNodes < InitStructNodes & PlanPerc > PlanVarLimitIC){

      print("Expanding initial line")

      # Contiune to add node untill the circle remains planar

      BasicLineData <- append(BasicLineData, computeElasticPrincipalGraph(Data = Data[rownames(Data) %in% NonG0Cell,],
                                                                          NumNodes = UsedNodes + 1,
                                                                          Method = 'CircleConfiguration',
                                                                          NodesPositions = BasicLineData[[length(BasicLineData)]]$Nodes,
                                                                          Edges = BasicLineData[[length(BasicLineData)]]$Edges))


      UsedNodes <- nrow(BasicLineData[[length(BasicLineData)]]$Nodes)

      PCAStruct <- prcomp(BasicLineData[[length(BasicLineData)]]$Nodes, center = TRUE, scale. = FALSE, retx = TRUE)
      PlanPerc <- sum(PCAStruct$sdev[1:2]^2)/sum(PCAStruct$sdev^2)

    }

    # Step II - Using curves

    while(UsedNodes < nNodes & PlanPerc > PlanVarLimit){

      print("Extending line")

      BasicLineData <- append(BasicLineData, computeElasticPrincipalGraph(Data = Data,
                                                                          NumNodes = nrow(BasicLineData[[length(BasicLineData)]]$Nodes)+1,
                                                                          Method = 'CurveConfiguration',
                                                                          NodesPositions = BasicLineData[[length(BasicLineData)]]$Nodes,
                                                                          Edges = BasicLineData[[length(BasicLineData)]]$Edges))


      Net <- ConstructGraph(Results = BasicLineData[[length(BasicLineData)]], DirectionMat = NULL, Thr = NULL)

      PCAStruct <- prcomp(BasicLineData[[length(BasicLineData)]]$Nodes, center = TRUE, scale. = FALSE, retx = TRUE)
      PlanPerc <- sum(PCAStruct$sdev[1:2]^2)/sum(PCAStruct$sdev^2)

      UsedNodes <- nrow(BasicLineData[[length(BasicLineData)]]$Nodes)

    }

    FitData <- BasicLineData
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

    Steps[[i]] <- ProjectAndCompute(DataSet = DataSet, GeneSet = UsedGenes[[i]], VarThr = VarThr, nNodes = nNodes,
                                    Log = Log, Categories = Categories, Filter = Filter, OutThr = OutThr, PCAFilter = PCAFilter,
                                    OutThrPCA = OutThrPCA, GraphType = GraphType, PlanVarLimit = PlanVarLimit, PlanVarLimitIC = PlanVarLimitIC,
                                    MinBranDiff = MinBranDiff, InitStructNodes = InitStructNodes, ForceLasso = ForceLasso, EstProlif = EstProlif,
                                    QuaThr = QuaThr, NonG0Cell = Steps[[1]]$NonG0Cell, DipPVThr = DipPVThr, MinProlCells = MinProlCells,
                                    PCACenter = PCACenter, PCAProjCenter = PCAProjCenter, PlotDebug = PlotDebug, PlotIntermediate = PlotIntermediate)

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

    Steps[[i]] <- ProjectAndCompute(DataSet = DataSet, GeneSet = UsedGenes[[i]], VarThr = VarThr, nNodes = nNodes,
                                    Log = Log, Categories = Categories, Filter = Filter, OutThr = OutThr, PCAFilter = PCAFilter,
                                    OutThrPCA = OutThrPCA, GraphType = GraphType, PlanVarLimit = PlanVarLimit, PlanVarLimitIC = PlanVarLimitIC,
                                    MinBranDiff = MinBranDiff, InitStructNodes = InitStructNodes, ForceLasso = ForceLasso, EstProlif = EstProlif,
                                    QuaThr = QuaThr, NonG0Cell = Steps[[1]]$NonG0Cell, DipPVThr = DipPVThr, MinProlCells = MinProlCells,
                                    PCACenter = PCACenter, PCAProjCenter = PCAProjCenter, PlotDebug = PlotDebug, PlotIntermediate = PlotIntermediate)

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








