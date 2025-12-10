library(reshape2)
library(ggplot2)
library(Rcpp)
library(RColorBrewer)
library(stringr)
library(data.table)

library(grid)
library(sp)

## Rcpp
sourceCpp("MatrixZero.cpp")

## Windows endings an begginnings localization
fnVentana <- function(nLargo,nVentana){
  Res <- nLargo%%nVentana
  Div <- nLargo%/%nVentana
  
  ## Result Table
  dtRes <- matrix(NA, nrow=nVentana, ncol=2)
  dtAux <- rep(0, nVentana)
  ## Ponemos los elementos de ajuste al final
  if (Res!=0) {
    dtAux[(length(dtAux)-Res+1):length(dtAux)] <- 1
  }
    
  #print(dtAux)
  
  dtRes[nVentana,2] <- nLargo
  dtRes[nVentana,1] <- nLargo-(Div-1) - dtAux[nVentana]
  
  for (i in (nVentana-1):1) {
    dtRes[i,2] <- dtRes[i+1,1]-1
    dtRes[i,1] <- dtRes[i,2]-(Div-1)-dtAux[i]
  }

  return(dtRes)
}

## Functions for analisis of missing data 
fnNASumUlt <- function(dtX,nLast) {
  return(sum(is.na(dtX[(length(dtX)-nLast+1):length(dtX)])))
}

fnNASumFirst <- function(dtX,nFirst) {
  return(sum(is.na(dtX[1:nFirst])))
}


fnNASum <- function(dtX) {
  return(sum(is.na(dtX)))
}

fnFullSum <- function(dtX) {
  return(sum(!is.na(dtX)))
}


fnLargoNA <- function (a){
  rl <- rle(is.na(c(a)))
  sal <- max(rl$lengths[rl$values])
  salF <- ifelse(is.infinite(sal),0,sal)
}

## Densities in specific points 
fnSumUltD <- function(dtX,nLast) {
  return(sum(dtX[(length(dtX)-nLast):length(dtX)], na.rm=TRUE))
}


## Extreme values labels 
fnAtipicoEtiq <- function(x,Bajo,Alto) {
    res <- ifelse(x<Bajo,-Inf,x)
    res <- ifelse(res>Alto,Inf,res)
    res
}


fnCeroSum <- function(dtX) {
        return(sum((dtX==0)))
}

## Function for counting elements in a segment

## Grouping functions
## The result of the grouping function is a label for each of the data
## There is a dependency in both grouping functions and ordering functions related
## to the number of vertical segments
fnGrouping = function(dtDatosGra, nSections=4, nGroups = 16,lgIndex=FALSE,lgNumGroups=FALSE) {
  dtLimites <- fnVentana(ncol(dtDatosGra[,-c(1:3)]),nSections)
  conShift <- 3
  ## Kmeans Class
  
  if (nSections>=4)
    dtSec4 <- apply(dtDatosGra[,(dtLimites[4,1]+conShift):(dtLimites[4,2]+conShift)],1,fnNASum)
  if (nSections>=3)
    dtSec3 <- apply(dtDatosGra[,(dtLimites[3,1]+conShift):(dtLimites[3,2]+conShift)],1,fnNASum)
  if (nSections>=2) {
    dtSec2 <- apply(dtDatosGra[,(dtLimites[2,1]+conShift):(dtLimites[2,2]+conShift)],1,fnNASum)
    dtSec1 <- apply(dtDatosGra[,(dtLimites[1,1]+conShift):(dtLimites[1,2]+conShift)],1,fnNASum)
  }
  ## Uso del elemento de la agrupacion usando los sectores definidos
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  
  if (nSections==4)
    indexT <- apply(cbind(dtSec1,dtSec2,dtSec3,dtSec4),2,range01)

  if (nSections==3)
    indexT <- apply(cbind(dtSec1,dtSec2,dtSec3),2,range01)
  
  if (nSections==2)
    indexT <- apply(cbind(dtSec1,dtSec2),2,range01)
  
  ## Evaluation
  if (lgNumGroups) {
    scaled_data = indexT
    fviz_nbclust(scaled_data, kmeans, method = "wss") +
      geom_vline(xintercept = 2, linetype = 2)+
      labs(subtitle = "Elbow method",iter.max = 50)
    # Silhouette method
    fviz_nbclust(scaled_data, kmeans, method = "silhouette")+
      labs(subtitle = "Silhouette method",iter.max = 50)
    # Gap statistic
    # nboot = 50 to keep the function speedy. 
    # recommended value: nboot= 500 for your analysis.
    # Use verbose = FALSE to hide computing progression.
    set.seed(123)
    fviz_nbclust(scaled_data, kmeans, nstart = 25,  method = "gap_stat", nboot = 500)+
      labs(subtitle = "Gap statistic method",iter.max = 50)
  }
  
  ## using grouping algorithm    
  nCluster <- nGroups
  groupsKT <- kmeans(indexT,centers = nCluster,nstart = 50, iter.max = 50)
  Grupo <- groupsKT$cluster
  
  ## Ordering algorithm for groups
  dtW <- c(1,4,16,64)
  if (lgIndex) {
    groups <- Grupo
    # lgOrdenGrupo <- FALSE
    # if (lgOrdenGrupo) {
    #   ## Indizado por grupos
    #   dtSumNAT <- apply(dtDatosKG[,-c(1:3)],1,fnFullSum)
    #   indexgroups <- tapply(dtSumNAT,groups,sum,na.rm=TRUE)#/(table(groups)*nnTime)
    #   dtOrdenNuevo <- order(indexgroups, decreasing=c("TRUE","TRUE","TRUE","TRUE"))
    # }
    
    ## Indizado por escalera
    lgEscalera <- TRUE
    if (lgEscalera) {
      if (nSections>=4)
        index4 <- tapply(dtSec4, groups,sum,na.rm=TRUE)/(table(groups)*(dtLimites[4,2] - dtLimites[4,1] +1))
      if (nSections>=3)
        index3 <- tapply(dtSec3, groups,sum,na.rm=TRUE)/(table(groups)*(dtLimites[3,2] - dtLimites[3,1] +1))
      if (nSections>=2) {
        index2 <- tapply(dtSec2, groups,sum,na.rm=TRUE)/(table(groups)*(dtLimites[2,2] - dtLimites[2,1] +1))
        index1 <- tapply(dtSec1, groups,sum,na.rm=TRUE)/(table(groups)*(dtLimites[1,2] - dtLimites[1,1] +1))
      }
      
      if (nSections==4){
        indexgroups <- index4
        indexa <- cbind(index1,index2,index3,index4)
      }
      
      if (nSections==3){
        indexgroups <- index3
        indexa <- cbind(index1,index2,index3)
      }
      
      if (nSections==2){
        indexgroups <- index2
        indexa <- cbind(index1,index2)
      }
      
      dtW = dtW[1:nSections]
      
      dtOrdenNuevo <- order(apply(dtW*t(indexa),2,sum),method="radix", decreasing=c("FALSE"))
    }
    
    dtTabla <- cbind(Nuevo=1:nCluster,Viejo=dtOrdenNuevo ,Etiqueta=as.numeric(names(indexgroups[dtOrdenNuevo])))
    groupsP <- dtTabla[match(groups,dtTabla[,"Etiqueta"]),"Nuevo"]
    groups <- groupsP
    Grupo <- groups
  }
  
  
  return(Grupo)
}

## Ordering functions
## The result of the ordering functions is an index for the data
fnOrdering = function(dtDatosOrd,nSections=4) {
  
  if (!is.null(nSections)) {
    dtLimites <- fnVentana(ncol(dtDatosOrd[,-c(1:3)]),nSections)
    conShift <- 3
  }
  
  ## Creation of a vector for ordering
  if (is.null(nSections)){
    idx <- order(dtDatosOrd[,"Grupo"],
                 apply(dtDatosOrd[,-c(1:3)],1,fnNASum),
                 apply(dtDatosOrd[,-c(1:3)],1,sum,na.rm=TRUE),
                 method="radix",
                 decreasing=c("FALSE","FALSE","TRUE"))
  } else {
  
    ## Statistics by section
    if (nSections>=4)
      dtSec4 <- apply(dtDatosOrd[,(dtLimites[4,1]+conShift):(dtLimites[4,2]+conShift)],1,fnNASum)
    if (nSections>=3)
      dtSec3 <- apply(dtDatosOrd[,(dtLimites[3,1]+conShift):(dtLimites[3,2]+conShift)],1,fnNASum)
    if (nSections>=2) {
      dtSec2 <- apply(dtDatosOrd[,(dtLimites[2,1]+conShift):(dtLimites[2,2]+conShift)],1,fnNASum)
      dtSec1 <- apply(dtDatosOrd[,(dtLimites[1,1]+conShift):(dtLimites[1,2]+conShift)],1,fnNASum)
    }
    
    ## index
    if (nSections==2){
      idx <- order(dtDatosOrd[,"Grupo"], dtSec2, dtSec1,
                   apply(dtDatosOrd[,-c(1:3)],1,fnNASum),
                   apply(dtDatosOrd[,-c(1:3)],1,sum,na.rm=TRUE),
                   method="radix",
                   decreasing=c("FALSE","FALSE","FALSE","FALSE","TRUE"))
    }
    
    if (nSections==3){
      idx <- order(dtDatosOrd[,"Grupo"], dtSec3, dtSec2, dtSec1,
                   apply(dtDatosOrd[,-c(1:3)],1,fnNASum),
                   apply(dtDatosOrd[,-c(1:3)],1,sum,na.rm=TRUE),
                   method="radix",
                   decreasing=c("FALSE","FALSE","FALSE","FALSE","FALSE","TRUE"))
    }
    
    if (nSections==4){
      idx <- order(dtDatosOrd[,"Grupo"], dtSec4, dtSec3, dtSec2, dtSec1,
                 apply(dtDatosOrd[,-c(1:3)],1,fnNASum),
                 apply(dtDatosOrd[,-c(1:3)],1,sum,na.rm=TRUE),
                 method="radix",
                 decreasing=c("FALSE","FALSE","FALSE","FALSE","FALSE","FALSE","TRUE"))
    }
  
  }
  
  return(idx)
  
}

## Summary functions
## The result of the summary functions is a function of the data (sum)
fnSummarySum = function(dtDatosOrd,nSections=4,lgGroup=TRUE) {
  
  if (!is.null(nSections)) {
    dtLimites <- fnVentana(ncol(dtDatosOrd[,-c(1:3)]),nSections)
    conShift <- 3
  }
  
  ## Creation of a summary vector per row
  if (is.null(nSections)){
    idx <- apply(dtDatosOrd[,-c(1:3)],1,sum,na.rm=TRUE)
  } else {
    
    ## Statistics by section
    if (nSections>=4)
      dtSec4 <- apply(dtDatosOrd[,(dtLimites[4,1]+conShift):(dtLimites[4,2]+conShift)],1,sum,na.rm=TRUE)
    if (nSections>=3)
      dtSec3 <- apply(dtDatosOrd[,(dtLimites[3,1]+conShift):(dtLimites[3,2]+conShift)],1,sum,na.rm=TRUE)
    if (nSections>=2) {
      dtSec2 <- apply(dtDatosOrd[,(dtLimites[2,1]+conShift):(dtLimites[2,2]+conShift)],1,sum,na.rm=TRUE)
      dtSec1 <- apply(dtDatosOrd[,(dtLimites[1,1]+conShift):(dtLimites[1,2]+conShift)],1,sum,na.rm=TRUE)
    }
  }
  
  ## Creating a structure for Sum
  if (nSections==2){
    idx <- cbind(Grupo=dtDatosOrd[,"Grupo"], dtSec2, dtSec1)
  }
  
  if (nSections==3){
    idx <- cbind(Grupo=dtDatosOrd[,"Grupo"], dtSec3, dtSec2, dtSec1)
  }
    
  if (nSections==4){
    idx <- cbind(Grupo=dtDatosOrd[,"Grupo"], dtSec4, dtSec3, dtSec2, dtSec1)
  }
  
  idx = as.data.frame(idx)
  
  if (lgGroup) {
    dtSal = aggregate(. ~ Grupo,idx,sum,na.rm=TRUE) 
  } else {
    dtSal = idx
  }
  
  return(dtSal)
  
}




## find square matrices of NA inside a matrix
## m1c is the matrix of data and minArea is a constant
fnPolCuadrados <- function(m1c,minArea) {
  nnren <- nrow(m1c)
  nncol <- ncol(m1c)
  m1c  <- as.matrix(m1c)
  
  ## Shadow matrix for squea algorithm (NA=0 y !NA=1)
  idx <- which(!is.na(m1c))
  m1c[idx] <- 1
  idx <- which(is.na(m1c))
  m1c[idx] <- 0
  
  
  ## Square algorithm
  conArea <- 0.005
  lsCuadradoTotal <- list()
  lsCuadradoArea <- list()
  for (i in 1:(nnren*nncol)) {
    ## i <- 1
    ## Find greatest submatrix (fn_zero_matrix is in MatrixZero.cpp)
    dtCuadrado <- fn_zero_matrix(m1c)
    
    if (dtCuadrado[1]<conArea*nncol*nnren)
      break(i)
    
    dtCoordenadas <- fn_cr_zero(dtCuadrado)
    colnames(dtCoordenadas) <- c("x", "y")
    
    idx  <- dtCoordenadas[1,"x"]:dtCoordenadas[3,"x"]
    idy  <- dtCoordenadas[1,"y"]:dtCoordenadas[2,"y"]
    
    m1c[idx,idy]  <- 1
    
    ## Save submatrix
    dtCoordenadas <- data.frame(dtCoordenadas,id=i)
    lsCuadradoTotal <- c(lsCuadradoTotal,list(dtCoordenadas))
    lsCuadradoArea <- c(lsCuadradoArea,dtCuadrado[1])
  }
  
  dtCuadradosTotal <- rbindlist(lsCuadradoTotal)
  dtCuadradosTotalArea <- c(unlist(lsCuadradoArea))
  
  
  if (!is.null(dtCuadradosTotalArea))
    return(list(Coor=dtCuadradosTotal,Areas=dtCuadradosTotalArea))
  else
    return(NULL)
}

## Find area in squares areas
fnAreaPol <- function(m1, dtCuadradosTotal, dtCuadradosTotalArea) {
  ## Shadow matrix (NA=1, other data=0)
  
  if(is.null(dtCuadradosTotalArea)) {
    return(NULL)
  }
  
  idx <- which(!is.na(m1))
  m1[idx] <- 0
  idx <- which(is.na(m1))
  m1[idx] <- 1
  
  ## Realizar el procedimiento anterior para los cuadrados
  dtCooNACuad <- which(m1==1, arr.ind=TRUE)
  dtCooNoNACuad  <- which(m1==0, arr.ind=TRUE)
  iCuad <- length(table(dtCuadradosTotal$id))
  
  sumTCuad <- rep(0, iCuad)
  for (i in 1:iCuad) {
    sump <- sum(point.in.polygon(dtCooNACuad[,2], dtCooNACuad[,1], unlist(dtCuadradosTotal[dtCuadradosTotal$id==i,"y"]),  unlist(dtCuadradosTotal[dtCuadradosTotal$id==i,"x"]))>0)
    sumTCuad[i] <- sump
  }
  
  sumTNCuad <- rep(0, iCuad)
  for (i in 1:iCuad) {
    sump <- sum(point.in.polygon(dtCooNoNACuad[,2], dtCooNoNACuad[,1], unlist(dtCuadradosTotal[dtCuadradosTotal$id==i,"y"]),  unlist(dtCuadradosTotal[dtCuadradosTotal$id==i,"x"]))>0)
    sumTNCuad[i] <- sump
  }
  
  
  dtPor0Cuad <- sum(sumTCuad)/nrow(dtCooNACuad)
  dtPor0NoCuad  <- sum(sumTNCuad)/nrow(dtCooNoNACuad)
  dtPorCuadSep <- dtCuadradosTotalArea/nrow(dtCooNACuad)
  return(list(Por0=dtPor0Cuad, Por0N=dtPor0NoCuad, PorCuad=dtPorCuadSep ))
  
}


## Conversion of result from values of squares function 
## to x,y coordinates
## dtVector has 5 elements, area, i, d[j], d2[j], d1[j]
fn_cr_zero <- function (dtVector) {
  dtSalida <- matrix(NA,ncol=2,nrow=4)
  i <- dtVector[2]
  d <- dtVector[3]
  d2 <- dtVector[4]
  d1 <- dtVector[5]
  
  ## Coordinates in 0
  xa <- xc <- d1 + 1
  xb <- xd <- d2 - 1
  yd <- yc <- i
  ya <- yb <- i-(i-d-1)
  
  ## Coordinates in 1
  dtSalida[1,] <- c(ya,xa)
  dtSalida[2,] <- c(yb,xb)
  dtSalida[3,] <- c(yd,xd)
  dtSalida[4,] <- c(yc,xc)
  dtSalida <- dtSalida + 1
  return(dtSalida)
}



## Graphical functions
## dtDatosGra matrix of data with a colunm of groups
## lgGroups: mark groups (horizontal) in the table
## conVentana: mark vertical lines 
## lsResCuadF: localization of "squares areas" inside the graph
## AtipRetail: vector with inferior and superior limits for atypical data
## idxOrder: vector with an ordering for plotting
fnGraficaVacio <- function(dtDatosGra, RutaSal,AtipRetail, lgGrupos=FALSE, conVentana=NULL, lsResCuadF=NULL, xlab=NULL,ylab=NULL,idxOrder=NULL, title=NULL,legend=NULL) {
  
  ## Ordenamiento
  if (!is.null(idxOrder))
    dtDatosGra[,"ClienteOrden"] <- idxOrder
  
  ## Cuadrados (verificar los case switch)
  m1 <- dtDatosGra[,4:ncol(dtDatosGra)] 
  conArea <- 0.001
  
  if (!is.null(lsResCuadF))
    lsResCuadFInt = fnPolCuadrados(m1,conArea)
  else 
    lsResCuadFInt = NULL
  
  lgCuadrados <- ifelse(is.null(lsResCuadFInt),FALSE,TRUE)
  lsPorCuad <- fnAreaPol(as.matrix(m1), lsResCuadFInt$Coor, lsResCuadFInt$Areas)
  
  ## Ventanas de datos
  if (!is.null(conVentana)) {
    ## Numero de ventanas
    dtVentana <- fnVentana(ncol(dtDatosGra[,-c(1:3)]),conVentana)
  }
  
  ## Se hacen 5 divisiones para la etiqueta de los datos
  nnVis <- nrow(dtDatosGra)
  nnPaso20 = nnVis/5
  dtEtiquetasPosY <- c(1,seq(from=nnPaso20, to=nnVis, by=nnPaso20))
  dtEtiquetasY <- paste(c(1,seq(from=20, to=100, by=20)),"%")
  
  nnSemanas <- ncol(dtDatosGra) - 3
  dtEtiquetasX <- c(1,seq(from=4, to=nnSemanas, by=4))
  
  dtDatosKGLim <- dtDatosGra
  dtDatosKGLim[,"ClienteOrden"] <- 1:nrow(dtDatosKGLim)
  
  dtGrupos <- table(dtDatosKGLim[,"Grupo"])
  dtGrupos <- c(1, dtGrupos)
  dtGrupos <- cumsum(dtGrupos)
  
  idxElim = which(is.element(colnames(dtDatosKGLim),c("Cliente","Grupo")))
  dtDatosKGLim  <- reshape2::melt(dtDatosKGLim[,-idxElim], id.vars = c("ClienteOrden"), measure.vars = colnames(dtDatosKGLim[,4:ncol(dtDatosKGLim)]))
  colnames(dtDatosKGLim) <- c("Cliente","Dia","Total")
  
  ## Atypical data labels
  dtDatosKGLim[,"Total"] <- fnAtipicoEtiq(dtDatosKGLim[,"Total"],AtipRetail[1],AtipRetail[2])
  
  idxIB <- which(dtDatosKGLim[,"Total"]==-Inf)
  idxIA <- which(dtDatosKGLim[,"Total"]==Inf)
  dtDatosKGLim[idxIB,"Total"] <- NA
  dtDatosKGLim[idxIA,"Total"] <- NA
  
  nnCortes <- 5
  dtDatosKGLim[,"Total"] <- cut(dtDatosKGLim[,"Total"], breaks=nnCortes, dig.lab = 2)
  
  ## Correct label
  dtDatosKGLim$Total <- factor(dtDatosKGLim$Total, levels = c("Atip_Low", levels(dtDatosKGLim$Total) ,"Atip_High"))
  
  dtDatosKGLim[idxIB,"Total"] <- "Atip_Low"
  dtDatosKGLim[idxIA,"Total"] <- "Atip_High"
  
  dtColores <- c("red", brewer.pal(n = nnCortes, name = "YlGn"),"purple")
  
  p <- ggplot(dtDatosKGLim, aes(Dia, Cliente, z= Total)) + geom_tile(aes(fill = Total)) +
    scale_fill_manual(na.value = '#87ceeb', values = dtColores)
    #scale_fill_continuous(na.value = '#FFFFFF',low = "#717171", high = "#000000")
  
  ## Labels for graphing
  conPercentMissing = ifelse(!is.null(lsPorCuad), round(lsPorCuad$Por0,3), 0)
  dtTitulo <- paste(title,"-----" ,"Fraction of Missing Data in square areas-",conPercentMissing,sep="" )
  p <- p + labs(x=xlab,y=ylab,title=dtTitulo) + theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  ## Percentage label needs more work (1000 only)
  #p <- p + labs(fill = "Sales($)") + scale_x_discrete(breaks=dtEtiquetasX) + scale_y_continuous(labels = function(x) paste0(x/10, "%"), sec.axis = dup_axis())
  p <- p + labs(fill = legend) + scale_x_discrete(breaks=dtEtiquetasX) + scale_y_continuous(labels = dtEtiquetasY, breaks= dtEtiquetasPosY, sec.axis = dup_axis())
  
  if (lgGrupos) {
    for (iLinea in 1:length(dtGrupos)){
      p <- p + geom_hline(yintercept=dtGrupos[iLinea], linetype="dashed", color = "black", size=0.25)
    }
  }
  
  if (!is.null(conVentana)) {
      for (iVentana in 1:(nrow(dtVentana)-1)) {
        p <- p + geom_vline(xintercept=dtVentana[iVentana,2]+0.5, linetype="dashed", color = "black", size=0.25)
      }
  }
  
  
  if (lgCuadrados) {
    ##  i <- 1
    dtCuadradosTotal <- lsResCuadFInt$Coor
    iCuad <- length(table(dtCuadradosTotal$id))
    for (i in 1:iCuad) {
      dtPol <- data.frame(Dia=unlist(dtCuadradosTotal[dtCuadradosTotal$id==i,"y"])+c(-0.5,0.5,0.5,-0.5),Cliente=unlist(dtCuadradosTotal[dtCuadradosTotal$id==i,"x"]),Total=0)
      Base <- (dtPol[2,"Dia"]-dtPol[1,"Dia"]) 
      Altura <- (dtPol[3,"Cliente"]-dtPol[1,"Cliente"])
      
      ##conTam <- 1000
      conTam = nnVis
      conAltura <- 31
      if (Altura > conAltura) {
        PosBase <- dtPol[1,"Dia"] + Base/2
        PosAltura <- dtPol[1,"Cliente"]  + Altura/2
        dtEtiquetaCuad <- data.frame(Dia=PosBase, Cliente=PosAltura,Total=0)
        p <- p + geom_polygon(data=dtPol, linetype="longdash", color = "orange3", size=0.8, fill=NA)
        p <- p + geom_label(data=dtEtiquetaCuad,label=paste(Base, paste(round((Altura/conTam)*100,2),"%",sep=""), sep="-"), size=2.2)
      }
    }
  }
  
  p <- p + theme(legend.position="bottom",legend.key.size = unit(0.5,"line")) + theme(legend.background = element_rect(fill="lightblue",
                                                                                    size=0.3, linetype="solid", 
                                                                                    colour ="darkblue"))
  #p <- p + theme(legend.text = element_text(size=8))
  
  p <- p + theme(plot.title = element_text(size=7),
                 legend.title = element_text(size = 5), 
                 legend.text = element_text(size = 7))
  
  Nombre <- RutaSal
  jpeg(Nombre, width = 4800, height = 3000, units = "px", res = 800)
  print(p)
  dev.off()
  
}

