
FindAllCorr <- function(Usable, Proc1, Proc2, Start1, Start2){

  CorVect <- NULL

  for(gId in 1:length(Usable)){
    X1 <- Proc1$NodesPT[Start1:length(Proc1$NodesPT)]/max(Proc1$NodesPT)
    X2 <- Proc2$NodesPT[Start2:length(Proc2$NodesPT)]/max(Proc2$NodesPT)

    X1 <- X1 - min(X1)
    X2 <- X2 - min(X2)

    X1 <- X1/max(X1)
    X2 <- X2/max(X2)

    Y1 <- Proc1$NodesExp[Usable[gId],Start1:length(Proc1$NodesPT)]
    Y1 <- Y1 - min(Y1)
    Y1 <- Y1/max(Y1)

    Y2 <- Proc2$NodesExp[Usable[gId],Start2:length(Proc2$NodesPT)]
    Y2 <- Y2 - min(Y2)
    Y2 <- Y2/max(Y2)

    f1 <- approxfun(X1, Y1, method = "linear")
    f2 <- approxfun(X2, Y2, method = "linear")

    CorVect <- c(CorVect, cor(f1(X1), f2(X1)))

    # plot(X1, Y1, col="blue", type = "b", main = paste(gId, Usable[gId]))
    # points(X2, Y2, col = "red", type = "b")

  }

  names(CorVect) <- Usable

  return(CorVect)

}


# SelGene <- Usable[4]
#
# par(mfcol = c(4,4))

Cor.Bue_Sasa <- FindAllCorr(Usable = intersect(rownames(Proc.Exp.Buet$NodesExp), rownames(Proc.Exp.Sasa$NodesExp)),
                            Proc1 = Proc.Exp.Buet, Proc2 = Proc.Exp.Sasa, Start1 = 1, Start2 = 1)

Cor.Bue_Kowa <- FindAllCorr(intersect(rownames(Proc.Exp.Buet$NodesExp), rownames(Proc.Exp.Kowa$NodesExp)),
                           Proc.Exp.Buet, Proc.Exp.Kowa, 1, max(Proc.Exp.Kowa$StageOnNodes$G0)+1)


Cor.Sasa_Kowa <- FindAllCorr(intersect(rownames(Proc.Exp.Sasa$NodesExp), rownames(Proc.Exp.Kowa$NodesExp)),
                            Proc.Exp.Sasa, Proc.Exp.Kowa, 1, max(Proc.Exp.Kowa$StageOnNodes$G0)+1)



hist(Cor.Bue_Sasa)
hist(Cor.Bue_Kowa)
hist(Cor.Sasa_Kowa)


plot(Cor.Bue_Sasa[intersect(names(Cor.Bue_Sasa), names(Cor.Bue_Kowa))],
     Cor.Bue_Kowa[intersect(names(Cor.Bue_Sasa), names(Cor.Bue_Kowa))])


plot(Cor.Bue_Sasa[intersect(names(Cor.Bue_Sasa), names(Cor.Sasa_Kowa))],
     Cor.Sasa_Kowa[intersect(names(Cor.Bue_Sasa), names(Cor.Sasa_Kowa))])

AllGenes <- unique(c(names(Cor.Bue_Sasa), names(Cor.Bue_Kowa), names(Cor.Sasa_Kowa)))

InAll <- AllGenes %in% names(Cor.Bue_Sasa) &
  AllGenes %in% names(Cor.Bue_Kowa) &
  AllGenes %in% names(Cor.Sasa_Kowa)


Combined <- rbind(Cor.Bue_Sasa[AllGenes[InAll]],
      Cor.Bue_Kowa[AllGenes[InAll]],
      Cor.Sasa_Kowa[AllGenes[InAll]])


which(apply(Combined, 2, min) > .7)


plot(sort(apply(Combined, 2, min)))


par(mfcol = c(1,1))

plot(CorVect, apply(Proc.Exp.Buet$NodesExp[Usable,], 1, max) - apply(Proc.Exp.Buet$NodesExp[Usable,], 1, min))

plot(sort(CorVect))

sum(CorVect > .5)/length(CorVect)


hist(CorVect)


gId <- 3

