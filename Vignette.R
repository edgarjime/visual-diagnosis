## Random numbers
library(openxlsx)
library(ineq)
library(factoextra)
library(NbClust)
source("Functions.R")

RutaData = getwd()
RutaIma = getwd()


## Online sales
NombreOnline = paste(RutaData, "/OnlineStore.csv", sep="")
dtDatosCliM = read.csv(NombreOnline)
AtipOnline <- quantile(unlist(dtDatosCliM[,-1]), c(0.01,0.99),na.rm=TRUE)

dtDatosKG <- data.frame(Cliente=dtDatosCliM[,"Cliente"], ClienteOrden=1:nrow(dtDatosCliM),
                        Grupo=NA,dtDatosCliM[,-1])
colnames(dtDatosKG) <- c("Cliente","ClienteOrden","Grupo",1:(ncol(dtDatosCliM)-1))

## Ordering using Total by client
dtNombreRetailOrd = paste(RutaIma,"/OnlineOrd.jpeg", sep="")
dtSumaInd = apply(dtDatosKG[,-c(1:3)],1,sum,na.rm=TRUE)
idx = order(dtSumaInd,decreasing=TRUE)
dtDatosKG = dtDatosKG[idx,]
fnGraficaVacio(dtDatosGra=dtDatosKG, RutaSal=dtNombreRetailOrd, AtipRetail=AtipOnline, lgGrupos=FALSE, conVentana=4, lsResCuadF=TRUE,
               xlab="Weeks",ylab="Percentage of clients",title="Online Sales per client and Week")

## Ordering using kmeans and radix sort
idxGrupo = fnGrouping(dtDatosKG,nSections = 4,nGroups = 8,lgIndex = TRUE,lgNumGroups=FALSE)
dtDatosKG$Grupo = idxGrupo
dtDatosKG = dtDatosKG[order(dtDatosKG$Grupo),]
dtNombreOnline = paste(RutaIma,"/Online8.jpeg", sep="")
idx = fnOrdering(dtDatosKG, nSections = 4)
dtDatosKG = dtDatosKG[idx,]
fnGraficaVacio(dtDatosGra=dtDatosKG, RutaSal=dtNombreOnline, AtipRetail=AtipOnline, lgGrupos=TRUE, conVentana=4, 
               lsResCuadF=TRUE, xlab="Weeks",ylab="Percentage of clients",title="Online Sales per client and Week")

## Summary graph
## Sum by periods, then by groups
dtSuma = fnSummarySum(dtDatosKG,nSections=4)
Total = sum(unlist(dtSuma[,-1]))
dtSumaPor = dtSuma[,-1]/Total

nnData = nrow(dtSumaPor)
dtPor= as.matrix(dtSumaPor)
dtPor= dtPor[nnData:1,] 
dtPor= dtPor[,ncol(dtPor):1] 

## Data for graphing
dtDataGroups = cbind(dtSumaPor,table(dtDatosKG$Grupo))
colnames(dtDataGroups) = c(paste("Segment_",4:1,sep=""), "Group","Size")
dtArchivoTabla = paste(RutaIma,"/Online8_Size.csv", sep="")
write.csv(dtDataGroups, dtArchivoTabla)

## Summary graphic routine
dtDataGroups = read.csv(dtArchivoTabla)
dtDataGroups  <- reshape2::melt(dtDataGroups, id.vars = c("Group","Size"), measure.vars = colnames(dtDataGroups[,2:5]))
dtDataGroups$variable = as.numeric(substr(as.character(dtDataGroups$variable),9,9))
dtDataGroups$value = dtDataGroups$value

## Defining sizes for the graph
dtSizes = dtDataGroups$Size[1:8]
dtLabels = paste(1:8, "-", dtSizes, sep="")
dtSizes = dtSizes/max(dtSizes)*10

numG = str_pad(as.character(dtDataGroups$Size),4,side = "left")
dtDataGroups$Group = paste(dtDataGroups$Group," - ",numG,sep="")

## Graphic elements
xlab="Period"
ylab="Percentage against Sales Total"
dtTitulo= "Online Sales - Percentage of sales by Period and Group"

p = ggplot(dtDataGroups, aes(x = variable, y = value,  group=Group, color = Group, size = Size,alpha=0.8))+     
 geom_path() + geom_point() 

p = p + facet_wrap(~ Group, ncol = 4) 

p = p+ scale_y_continuous(labels = scales::percent)
p = p + labs(x=xlab,y=ylab,title=dtTitulo)
p = p  + scale_color_viridis_d() 
p = p + guides(color = guide_legend(title = "Group - # Clients", override.aes = list(lwd=dtSizes)),
              size = "none",
              alpha = "none")
p =p +  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.3, face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "none" 
  ) 

Nombre <- paste(RutaIma,"/Online8_Size_VF.jpeg", sep="")
jpeg(Nombre, width = 4800, height = 3000, units = "px", res = 800)
print(p)
dev.off()






