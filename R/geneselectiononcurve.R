#' Title
#'
#' @param TaxVect
#' @param ExpMat
#' @param Thr
#' @param MadsThr
#' @param Mode
#' @param Net
#' @param PriGraph
#' @param Proj
#'
#' @return
#' @export
#'
#' @examples
SelectGenes <- function(TaxVect,
                        ExpMat,
                        Mode = "CV",
                        Net = NULL,
                        PriGraph = NULL,
                        Proj = NULL,
                        AggFun = min) {

  if(Mode == "CV" | is.null(Net) | is.null(PriGraph) | is.null(Proj)){

    SelCell <- names(TaxVect)

    GeneExprMat.DF.Split <- split(data.frame(ExpMat[SelCell, ]), TaxVect[SelCell])

    GeneExprMat.DF.Split.Mean <- lapply(GeneExprMat.DF.Split, colMeans)
    GeneExprMat.DF.Split.Sd <- lapply(GeneExprMat.DF.Split, function(x){apply(x, 2, sd)})

    GeneExprMat.DF.MeanRemoved <-
      lapply(as.list(1:length(GeneExprMat.DF.Split)), function(i){
        GeneExprMat.DF.Split.Sd[[i]]/GeneExprMat.DF.Split.Mean[[i]]
      })

    GeneExprMat.DF.MeanRemoved.All <- do.call(rbind, GeneExprMat.DF.MeanRemoved)

    return(apply(abs(GeneExprMat.DF.MeanRemoved.All), 2, AggFun, na.rm = TRUE))

  }

  if(Mode == "SmoothOnCircleNodes"){

    NodeOrder <- igraph::subgraph_isomorphisms(target = Net, pattern = igraph::make_lattice(igraph::vcount(Net), circular = TRUE))[[1]] %>%
      names(.) %>%
      strsplit(., split = "_", fixed = TRUE) %>%
      sapply(., "[[", 2)

    OrderedPoints <- OrderOnPath(PrinGraph = PriGraph, Path = as.numeric(c(NodeOrder, NodeOrder[1])), PointProjections = Proj)

    ExtNodeOrder <- rep(NodeOrder, 3)

    names(OrderedPoints$PositionOnPath) <- rownames(ExpMat)
    ExtPT <- c(OrderedPoints$PositionOnPath - sum(OrderedPoints$PathLen),
            OrderedPoints$PositionOnPath,
            OrderedPoints$PositionOnPath + sum(OrderedPoints$PathLen))
    NodePoints <- cumsum(OrderedPoints$PathLen)

    AllFit <- apply(ExpMat, 2, function(x){

      LocDF <- data.frame(Exp = rep(x, 3), PT = ExtPT)

      LOE <- loess(Exp ~ PT, LocDF)
      predict(LOE, data.frame(PT = NodePoints), se = TRUE)

    })

    RetVal <- lapply(AllFit, function(x){x$se.fit/x$fit}) %>%
      sapply(., AggFun)

    return(RetVal)

  }

  if(Mode == "LinearOnCircleNodes"){

    NodeOrder <- igraph::subgraph_isomorphisms(target = Net, pattern = igraph::make_lattice(igraph::vcount(Net), circular = TRUE))[[1]] %>%
      names(.) %>%
      strsplit(., split = "_", fixed = TRUE) %>%
      sapply(., "[[", 2)

    OrderedPoints <- OrderOnPath(PrinGraph = PriGraph, Path = as.numeric(c(NodeOrder, NodeOrder[1])), PointProjections = Proj)

    ExtNodeOrder <- rep(NodeOrder, 3)

    names(OrderedPoints$PositionOnPath) <- rownames(ExpMat)
    ExtPT <- c(OrderedPoints$PositionOnPath - sum(OrderedPoints$PathLen),
               OrderedPoints$PositionOnPath,
               OrderedPoints$PositionOnPath + sum(OrderedPoints$PathLen))
    NodePoints <- cumsum(OrderedPoints$PathLen)
    Borders <- cumsum(OrderedPoints$PathLen)

    Filtered <- (ExtPT <= sum(OrderedPoints$PathLen) & ExtPT >= 0)
    X <- ExtPT[Filtered]

    AllFit <- apply(ExpMat, 2, function(x){

      Y <- rep(x, 3)[Filtered]

      sapply(2:length(Borders), function(i){

        ToUse <- (X <= Borders[i] & X >= Borders[i-1])

        tX <- X[ToUse]
        tY <- Y[ToUse]

        if(sum(ToUse)>3 &
           length(unique(tX))>1 &
           length(unique(tY))>1 ){
          cor.test(tY, tX, method = "spe")$p.value
          } else {
          return(1)
        }

      }) %>% AggFun

    })

    return(AllFit)

  }

}

