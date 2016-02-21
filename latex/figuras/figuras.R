require(org.At.tair.db)
#Aplico filtros de evidencia
source("/home/arabinov/At/ecFilter.R")
dropEvidence <- c("IEA", "NAS", "IC", "ND", "RCA", "IEP")

initECfilter("org.At.tair",verbose=TRUE)
setECfilter(dropEvidence,verbose=TRUE)

#Leemos los archivos de datos
data_original <- read.table("/home/arabinov/At/AtGE_Abiostress_gcRMA.txt",header=TRUE, as.is=TRUE)
treatment     <- read.table("/home/arabinov/At/AtGE-Abiostress-sampleList.cvs",header=TRUE, as.is=TRUE)

require(csbl.go)
set.prob.table(filename="/home/arabinov/R/i686-pc-linux-gnu-library/3.2/csbl.go/auxdata/table-3702-no_IEA_NAS_IC_ND_RCA-similarity.csv")
ent<-entities.from.text("/home/arabinov/R/i686-pc-linux-gnu-library/3.2/csbl.go/auxdata/annotationTable.3702.no_IEA_NAS_IC_ND_RCA.csv")

#generamos data con los nombres de fila dados por el probeset
data           <- data_original
rownames(data) <- data_original[, 1]

#Los symbols los traemos del env del chipset y filtramos las que no apunten a nada (NA)
require(ath1121501.db)
symbols             <- mget(rownames(data), ath1121501ACCNUM, ifnotfound = NA)
u_symbols           <- unlist(symbols)
u_symbols_filtrados <- u_symbols[!is.na(u_symbols)]
data                <- (data[names(u_symbols_filtrados), ])

#Mostramos los tratamientos disponibles
rownames(treatment) <- treatment[,1]
table(treatment[,"TYPE"])

#traigo cada categoría de GO con los genes que tiene anotados
cg <- as.list(org.At.tairGO2ALLTAIRS)
# calculo el information content de todos, con GO:0008150 es la cantidad de anotaciones que tiene el nodo raiz BP
#ic <- unlist(lapply(cg, function(x) length(unique(x))))/length(unique(cg[["GO:0008150"]]))

#SQLITE
dbDisconnect(con)
con = dbConnect(RSQLite::SQLite(), dbname="At/db/at.sqlite") #abrimos la base de datos

library(ggplot2)

#Usamos como tratamiento de muestra Cold
tratamiento = "Cold"

#--------------------------------FIGURA 5----------------------------------------------
#Criterios de filtrado de genes
tratamientos <- treatment[treatment[,"TYPE"]==tratamiento & treatment[,"SRC"]=="Roots","ID"]
ipar <- seq(2,length(tratamientos),2)
inon <- ipar-1
atge <- 0.5*(data[,tratamientos[ipar]]+data[,tratamientos[inon]])

gsd   <- apply(atge, 1 ,sd)
plot(ecdf(gsd), main="Distribución acumulada de la desviación estandar de expresión genética", ylab="Probabilidad", xlab="Desviación estandar") #empirical cumulative distribution function de las desviaciones estandar
abline(v=quantile(gsd, 0.9), col="red") #En quantile(gsd, 0.9) tengo un 10% de los genes

#Criterio usual para decidir si hay senal o ruido
plot(density(as.matrix(atge)), main="Kernel density estimation para la distribución de los niveles de expresión genética", xlab="Expresión genética", ylab="Densidad")
abline(v=4, col="red")
#--------------------------------FIGURA 5----------------------------------------------


#--------------------------------FIGURA 4----------------------------------------------
#Perfiles de clusters para un tratamiento
query = paste("SELECT DISTINCT genes.probe, genes.cluster FROM genes INNER JOIN clusters ON genes.cluster = clusters.cluster INNER JOIN tratamientos ON genes.tratamiento_id = tratamientos.tratamiento_id WHERE tratamientos.tratamiento = '", tratamiento, "' AND tratamientos.deepsplit = 1", sep="")
clus<-dbGetQuery(con, query)
tratamientos <- treatment[treatment[,"TYPE"]=="Cold" & treatment[,"SRC"]=="Roots","ID"]
genes <- t(apply(data[clus[, "probe"], tratamientos], 1, function(x){ return( (x-mean(x))/sd(x) ) }))
ipar <- seq(2,length(tratamientos),2)
inon <- ipar-1
genes <- 0.5*(genes[,tratamientos[ipar]]+genes[,tratamientos[inon]])
colnames(genes) <- treatment[(treatment[,"TYPE"] == tratamiento & treatment[,"SRC"]=="Roots"), "TIME"][ipar]
layout(matrix(1:9, ncol=3))
for(i in sort(unique(clus[, "cluster"]))){
  matplot(t(genes[which(clus[, "cluster"]==i), ]), type="b", main=paste("Cluster ", i, " para ", tratamiento, " y ds=1 (", sum(clus[, "cluster"]==i), " genes)", sep=""), xlab = "Tiempo (hs.)", ylab = "Nivel de expresión (normalizado)", xaxt="n")
  axis(1, at = 1:ncol(genes), labels = colnames(genes))
}
layout(matrix(1, ncol=1, nrow=1))
#--------------------------------FIGURA 4----------------------------------------------


#--------------------------------FIGURA 1----------------------------------------------
#Inserta en base de datos la cohesion por tratamiento y por deesplit.
for(deepsplit in c(1, 4)){
  for(tratamiento in unique(treatment[, "TYPE"])){
    print(tratamiento)
    cohesion_por_tamano(con, tratamiento, data, treatment, deepsplit = deepsplit)
  }
}

#Muestra Cohesion deepsplit 1
deepsplit = 1
query = paste("SELECT correlaciones.correlacion, correlaciones.cuartil_1, correlaciones.cuartil_2, correlaciones.cuartil_3, clusters.total AS tamanio FROM correlaciones INNER JOIN clusters ON correlaciones.cluster_id = clusters.cluster_id INNER JOIN tratamientos ON clusters.tratamiento_id = tratamientos.tratamiento_id WHERE tratamientos.deepsplit = ", deepsplit, "  ORDER BY total ASC", sep="")
correlaciones<-dbGetQuery(con, query)

q <- ggplot(correlaciones, aes(tamanio, correlacion))
r <- geom_ribbon(data = correlaciones, aes(ymin=cuartil_1, ymax=cuartil_3), alpha=I(1/10))
l <- geom_line(data = correlaciones, aes(y=cuartil_2))
s <- stat_bin2d(binwidth = c(10, .03))
q+r+l+s+ggtitle(paste("Correlacion media por tamaño para todos los tratamientos con deepsplit =", deepsplit))

#Muestra Cohesion deepsplit 4
deepsplit = 4
query = paste("SELECT correlaciones.correlacion, correlaciones.cuartil_1, correlaciones.cuartil_2, correlaciones.cuartil_3, clusters.total AS tamanio FROM correlaciones INNER JOIN clusters ON correlaciones.cluster_id = clusters.cluster_id INNER JOIN tratamientos ON clusters.tratamiento_id = tratamientos.tratamiento_id WHERE tratamientos.deepsplit = ", deepsplit, "  ORDER BY total ASC", sep="")
correlaciones<-dbGetQuery(con, query)

q <- ggplot(correlaciones, aes(tamanio, correlacion))
r <- geom_ribbon(data = correlaciones, aes(ymin=cuartil_1, ymax=cuartil_3), alpha=I(1/10))
l <- geom_line(data = correlaciones, aes(y=cuartil_2))
s <- stat_bin2d(binwidth = c(5, .03))
q+r+l+s+ggtitle(paste("Correlacion media por tamaño para todos los tratamientos con deepsplit =", deepsplit))
#--------------------------------FIGURA 1----------------------------------------------

#--------------------------------FIGURA 2----------------------------------------------
#Histograma de tamaños de cluster por método
#Hist de dtc con ds1
query = paste("SELECT clusters.total FROM clusters WHERE clusters.tratamiento_id >=1 AND  clusters.tratamiento_id <=11", sep="")
clus<-dbGetQuery(con, query)
hist(as.integer(as.matrix(clus["total"])), breaks=20, main="Histograma de tamaños de cluster para DTC con DS=1", xlab="Tamaño de cluster", ylab="Frecuencia")

#Hist de dtc con ds2
query = paste("SELECT clusters.total FROM clusters WHERE clusters.tratamiento_id >=12 AND  clusters.tratamiento_id <=22", sep="")
clus<-dbGetQuery(con, query)
hist(as.integer(as.matrix(clus["total"])), breaks=20, main="Histograma de tamaños de cluster para DTC con DS=2", xlab="Tamaño de cluster", ylab="Frecuencia")
#--------------------------------FIGURA 2----------------------------------------------

#--------------------------------FIGURA 3----------------------------------------------
#Interacting densities
query = paste("SELECT interacting_densities.* FROM interacting_densities INNER JOIN tratamientos ON interacting_densities.tratamiento_id = tratamientos.tratamiento_id WHERE interacting_densities.categoria_interacting_density_id = (SELECT categoria_interacting_density_id FROM categorias_interacting_densities WHERE categoria_interacting_density LIKE 'DTC%_DS1') AND interacting_densities.tamano_total > 1", sep="")
ID<-dbGetQuery(con, query)
#cid<-t(apply(ID, 1, function(x) {calc_id(x["cantidades"], as.integer(x["tamano_total"]), as.integer(x["tamano_prefiltrado"]))}))
cID<-t(apply(ID, 1, function(x) {calc_ID(x["cantidades"], as.integer(x["tamano_total"]), as.integer(x["tamano_prefiltrado"]))}))
ID<-cbind(ID, cID)
colnames(ID) <- c("interacting_density_id", "cantidades", "tamano_total", "tamano_prefiltrado", "termino", "tratamiento_id", "categoria_interacting_density_id", "ID", "ID_prefiltrado")
B<-binning(ID[, c("ID", "tamano_total")], mean)
ID_DTCDS1 <- B

query = paste("SELECT interacting_densities.* FROM interacting_densities INNER JOIN tratamientos ON interacting_densities.tratamiento_id = tratamientos.tratamiento_id WHERE interacting_densities.categoria_interacting_density_id = (SELECT categoria_interacting_density_id FROM categorias_interacting_densities WHERE categoria_interacting_density LIKE 'DTC%_DS4') AND interacting_densities.tamano_total > 1", sep="")
ID<-dbGetQuery(con, query)
#cid<-t(apply(ID, 1, function(x) {calc_id(x["cantidades"], as.integer(x["tamano_total"]), as.integer(x["tamano_prefiltrado"]))}))
cID<-t(apply(ID, 1, function(x) {calc_ID(x["cantidades"], as.integer(x["tamano_total"]), as.integer(x["tamano_prefiltrado"]))}))
ID<-cbind(ID, cID)
colnames(ID) <- c("interacting_density_id", "cantidades", "tamano_total", "tamano_prefiltrado", "termino", "tratamiento_id", "categoria_interacting_density_id", "ID", "ID_prefiltrado")
B<-binning(ID[, c("ID", "tamano_total")], mean)
ID_DTCDS4 <- B

query = paste("SELECT interacting_densities.* FROM interacting_densities INNER JOIN tratamientos ON interacting_densities.tratamiento_id = tratamientos.tratamiento_id WHERE interacting_densities.categoria_interacting_density_id = (SELECT categoria_interacting_density_id FROM categorias_interacting_densities WHERE categoria_interacting_density LIKE 'KEGG%') AND interacting_densities.tamano_total > 1", sep="")
ID<-dbGetQuery(con, query)
#cid<-t(apply(ID, 1, function(x) {calc_id(x["cantidades"], as.integer(x["tamano_total"]), as.integer(x["tamano_prefiltrado"]))}))
cID<-t(apply(ID, 1, function(x) {calc_ID(x["cantidades"], as.integer(x["tamano_total"]), as.integer(x["tamano_prefiltrado"]))}))
ID<-cbind(ID, cID)
colnames(ID) <- c("interacting_density_id", "cantidades", "tamano_total", "tamano_prefiltrado", "termino", "tratamiento_id", "categoria_interacting_density_id", "ID", "ID_prefiltrado")
B<-binning(ID[, c("ID", "tamano_total")], mean)
ID_PIN <- B

query = paste("SELECT interacting_densities.* FROM interacting_densities INNER JOIN tratamientos ON interacting_densities.tratamiento_id = tratamientos.tratamiento_id WHERE interacting_densities.categoria_interacting_density_id = (SELECT categoria_interacting_density_id FROM categorias_interacting_densities WHERE categoria_interacting_density LIKE 'PIN%') AND interacting_densities.tamano_total > 1", sep="")
ID<-dbGetQuery(con, query)
#cid<-t(apply(ID, 1, function(x) {calc_id(x["cantidades"], as.integer(x["tamano_total"]), as.integer(x["tamano_prefiltrado"]))}))
cID<-t(apply(ID, 1, function(x) {calc_ID(x["cantidades"], as.integer(x["tamano_total"]), as.integer(x["tamano_prefiltrado"]))}))
ID<-cbind(ID, cID)
colnames(ID) <- c("interacting_density_id", "cantidades", "tamano_total", "tamano_prefiltrado", "termino", "tratamiento_id", "categoria_interacting_density_id", "ID", "ID_prefiltrado")
B<-binning(ID[, c("ID", "tamano_total")], mean)
ID_KEGG <- B


query = paste("SELECT interacting_densities.* FROM interacting_densities INNER JOIN tratamientos ON interacting_densities.tratamiento_id = tratamientos.tratamiento_id WHERE interacting_densities.categoria_interacting_density_id = (SELECT categoria_interacting_density_id FROM categorias_interacting_densities WHERE categoria_interacting_density LIKE '%random%') AND interacting_densities.tamano_total > 1", sep="")
ID<-dbGetQuery(con, query)
#cid<-t(apply(ID, 1, function(x) {calc_id(x["cantidades"], as.integer(x["tamano_total"]), as.integer(x["tamano_prefiltrado"]))}))
cID<-t(apply(ID, 1, function(x) {calc_ID(x["cantidades"], as.integer(x["tamano_total"]), as.integer(x["tamano_prefiltrado"]))}))
ID<-cbind(ID, cID)
colnames(ID) <- c("interacting_density_id", "cantidades", "tamano_total", "tamano_prefiltrado", "termino", "tratamiento_id", "categoria_interacting_density_id", "ID", "ID_prefiltrado")
B<-binning(ID[, c("ID", "tamano_total")], mean)
ID_RANDOM <- B

options(scipen=5)
plot(ID_DTCDS1[, 2], ID_DTCDS1[, 1], log="xy", main="Densidad de interacción", type="b", xaxt = "n", ylim=c(0.00015,0.105), col="red", xlab="Cantidad de anotaciones por término", ylab="Densidad de interacción", pch=1)
axis(side = 1, at = ID_DTCDS1[, 2], labels = ID_DTCDS1[, 2], tck=-.02)
points(ID_DTCDS4[, 2], ID_DTCDS4[, 1], type="b", col="blue", pch=2)
points(ID_PIN[, 2], ID_PIN[, 1], type="b", col="green", pch=3)
points(ID_KEGG[, 2], ID_KEGG[, 1], type="b", col="cyan", pch=4)
points(ID_RANDOM[, 2], ID_RANDOM[, 1], type="b", col="black", pch=5)
legend(list(x = 120,y = 0.13), title = "Algoritmos de agrupamiento", legend = c("Dynamic Tree Cut (Deepsplit 1)", "Dynamic Tree Cut (Deepsplit 4)", "Infomap (Protein Interaction Network)", "Infomap (Methabolyc Pathways)", "Control aleatorio"), col = c("red", "blue", "green", "cyan", "black"), pch = 1:5, lty = 1, merge = TRUE)
#--------------------------------FIGURA 3----------------------------------------------
