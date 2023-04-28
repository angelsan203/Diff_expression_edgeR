packages <- c("Biocmnanager","pheatmap", "edgeR", "RColorBrewer","colorspace","")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  BiocManager::install(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = TRUE))

directory <- "~/Documentos1/GF-BS/Reto_GFBS1"
setwd(directory);	


raw.table<-read.delim(file.choose(), stringsAsFactors=F);	#Table with counts
conditions<-read.delim(file.choose(), row.names=1);			#Experimental conditions

#Visualize if data was loaded correctly
dim(raw.table);	
raw.table[1:2,];	

#OK, so we need a table with just the counts
counts<-raw.table[,-c(1:6)];	#So we'll take everything except the first 6 columns
rownames(counts)<-raw.table[,"Geneid"];	#And add the IDs as rownames

#So, what does it look like?
counts[1:2,];

#Names don't tell me much
#But that's why we have a file with the experimental conditions
conditions[,"name"]; 
#let's make sure the columns match the order of the rows in conditions
conditions<-conditions[colnames(counts),]; 
#and re-name the columns
colnames(counts)<-conditions[colnames(counts),"name"]; 

#let's take a look at the counts
summary(counts)

#Which are the most expressed genes?
apply(counts,2,which.max)

rownames(counts)[unique(apply(counts,2,which.max))]
##Don't forget these are transcripts, with a biological meaning
#See if it makes sense: https://www.kegg.jp/brite/ang00001+ANI_1_760064


#Cómo se ven los datos?
boxplot(counts+1, log = "y", col = conditions$color, las = 2, main = "Distribución de muestras con posibles outliers", ylab = "counts", pch = 20, cex.axis = 0.7);


###vamos a empezar por un pequeño análisis explorativo
pheatmap(cor(log(counts+1,2)), main = "Heatmap explorativo con posibles outliers", color = colorRampPalette(rev(brewer.pal(n = 7, name = "Spectral")))(100)); #Hay cosas que se ven raras
#¿por qué? Porque son datos reales. Por eso. 

#Pero si se acuerdan, los análisis explorativos no los hacemos con counts
#Hay que transformar primero; con edgeR, primero hay que poner los datos en el formato adecuado: DGEList
?DGEList

gene.list<-DGEList(counts=counts, group=conditions$group)

#Then, we filter out genes with very low counts. We do this because of 2 things: 
#	(1)	To remove likely artifacts
#	(2)	Genes with low counts mess up the statistics (the model does not handle discreetnes well, curiously enough)
keep<-filterByExpr(gene.list,min.count=5,min.prop=0.5,group=conditions$group);

#then we just remove the genes that do not meet the criteria, and update the library sizes
gene.list<-gene.list[keep,,keep.lib.sizes=FALSE];

#Now we calculate the normalization factors (remember all the equations in the slides)
#Why do we need to calculate the normalization factor?
gene.list<-calcNormFactors(gene.list)

##Let's se what we get
boxplot(cpm(gene.list, log=TRUE), col = conditions$color, las = 2, main = "Distribución de muestras", ylab = "log cpm (TMM normalized)", cex.axis = 0.7);
boxplot(counts+1, log = "y", col = conditions$color, las = 2, main = "Distribución de muestras con factor de normalización", ylab = "counts", cex.axis = 0.7, pch = 20);


##Now, we explore the data. 
#MDS (Multi-dimensional scaling) is something like a PCA - If you ask me to elaborate I will take some time
#Pero las distancias son (roughly) leading log-fold-changes
#The leading log-fold-change is the average (root-mean-square) of the largest absolute log-foldchanges between each pair of samples
plotMDS(gene.list)

#This is the closest to a PCA
plotMDS(gene.list,gene.selection="common", pch=conditions$treatment, col = conditions$color);

plotMDS(gene.list, col = conditions$color, main = "Escalado multi-dimensional con posibles outliers");

pheatmap(cor(cpm(gene.list, log = TRUE)), color = colorRampPalette(rev(brewer.pal(n = 7, name = "Spectral")))(100)); #Definitivamente las muestras parecen ser outliers

#### Hay 2 opciones: 
#(1)	Seguir el análisis así y esperar que la estadística se arregle
#(2)	Quitar las muestras problemáticas como "outliers"

###	Para esto hay que leer un poco de los métodos que estamos utilizando; no necesariemente entender todas las matemáticas, pero entender el principio

###		edgeR lidia muy mal con outliers, y sí puede manejar diseños no balanceados (aunque no sea ideal)
#####	No siempre es el caso
#####	¿en qué punto borramos las muestras?

borrar<-c(1,13,16);

counts[1:2,borrar]
conditions[borrar,]

counts<-counts[,-borrar];
conditions<-conditions[-borrar,];

gene.list<-DGEList(counts=counts, group=conditions$group);
keep<-filterByExpr(gene.list,min.count=5,min.prop=0.6,group=conditions$group);	###LAS CONDICIONES CAMBIAN; ¿porqué?

gene.list<-gene.list[keep,,keep.lib.sizes=FALSE];
gene.list<-calcNormFactors(gene.list)

boxplot(cpm(gene.list, log=TRUE), col = conditions$color, las = 2, main = "Distribución de muestras", ylab = "log cpm (TMM normalized)");
boxplot(counts+1, log = "y", col = conditions$color, las = 2, main = "Distribución de muestras normalizadas", ylab = "counts", cex.axis = 0.7, pch = 20);

plotMDS(gene.list,gene.selection="common", col = conditions$color, main = "Escalado multi-dimensional normalizado");	#Meh. Mejor, pero no demasiado.

pheatmap(cor(cpm(gene.list, log = TRUE)), color = colorRampPalette(rev(brewer.pal(n = 7, name = "Spectral")))(100), main = "Heatmap explorativo con datos normalizados"); #Definitivamente las muestras parecen ser outliers

#Le alimentamos al modelo el diseño experimental
design<-model.matrix(~0+conditions$group)	#Esto elimina el intercepto

#¿qué es la matriz design?
t(design)

##Y estimamos la dispersión
gene.list<-estimateDisp(gene.list,design);
##this estimates the quantile-adjusted conditional maximum likelihood (qCML) dispersion for A vs. B experiments
##And for multiple factors, it uses Cox-Reid profile-adjusted likelihood (CR) - approximate conditional likelihoods
##In our case ¿which one is it calculating?

#Esta función genera seudo-conteos, que son conteos normalizados: 
#los números de cuentas que hubieran resultado si todas las librerías fueran iguales. 
#Sobre estas seudo-cuentas es que edgeR hace los cálculos

##Now finally we calculate the DE genes
#by Quasi Likelihood F-test	- There's also the option of doing this with an exact test that is similar to Fisher's exact test
#you can find out which one works better. 
QL.fit<-glmQLFit(gene.list,design)
#these actually fit a log-linear model
#which is equivalent to take the log of both dependent and independent variables (like we did in the slides - logs are awesome!)

#Remember the case with the ducks?

glmTest<-glmQLFTest(QL.fit, poisson.bound=FALSE);	#Do note that it is the log (base 2) fold-change  
topTags(glmTest)

#Great! We have detected some changes that are significantly different...
#...
#...
#from what?

#We need to specify pairwise comparisons
colnames(design)
glmTest<-glmQLFTest(QL.fit, contrast=c(1,0,0,0,0,0,-1), poisson.bound=FALSE);
topTags(glmTest)

##Not the most elegant solution, but it helps avoid the manual work:
valores<-c(-1,0,1);
contrastes<-expand.grid(rep(list(valores),ncol(design)));
contrastes<-contrastes[(rowSums(contrastes)==0)&(rowSums(abs(contrastes))==2),];
colnames(contrastes)<-colnames(design);

#This is definitely not elegant, and a good case for using pipes
rownames(contrastes)<-gsub("conditions[$]group","",apply(contrastes,1,function(x) {paste0(names(x)[x==1],".vs.",names(x)[x==-1])}))

#but now we can specify better which contrasts to make
glmTest<-glmQLFTest(QL.fit, contrast=as.numeric(contrastes["CTRL.0.vs.CTRL.2h",]), poisson.bound=FALSE);
topTags(glmTest1)

plotMD(glmTest1)

#Now, we DO have a threshold above which we consider the effect to be meaningful
#La socio formadora nos dijo un umbral de fold-change; ¿se acuerdan cuál era?
thr<- 1.2;


##We can make a heatmap

log.cpm<-cpm(gene.list, log = TRUE);

#Also, we can just export the whole thing to an Excel compatible file
full.table<-topTags(glmTest, n = "Inf")$table;

nrow(contrastes)

#actually 
nrow(contrastes)/2

#which ones will give us valuable, biologically relevant information?

##let's go back to the drawing board 
##(Literally, not figuratively)


#QLF works well to control type I error ("false positives") 
#But when the dispersions are very large and the counts are very small, the QL framework doesn't work that well
#It also doesn't work well with 2 replicates (when you have no replicates, it doesn't work at all)
#edgeR also offers the LRT (likelihood-ratio test) 
#Which may not control that well type I error, but is more flexible when you want to control type II error ("false negatives")
#We can also try that if we have time
glmTest1<-glmQLFTest(QL.fit, contrast=as.numeric(contrastes["CTRL.2h.vs.CTRL.8h",]), poisson.bound=FALSE);
topTags(glmTest1) 
GLM.fit<-glmFit(gene.list,design)
library("RColorBrewer")
library("VennDiagram")

chosen.qlf1 <-decideTests(glmTest1,adjust.method="BH", p.value=0.05,lfc=log(thr,2));
summary(chosen.qlf1);

up.regulated.qlf<-rownames(chosen.qlf1)[chosen.qlf1==1];
down.regulated.qlf<-rownames(chosen.qlf1)[chosen.qlf1==-1];


lrtTest1<-glmLRT(GLM.fit,contrast=as.numeric(contrastes["CTRL.2h.vs.CTRL.8h",]));
chosen.lrt1<-decideTests(lrtTest1,adjust.method="BH", p.value=0.05,lfc=log(thr,2));
summary(chosen.lrt1)

up.regulated.lrt1<-rownames(chosen.lrt1)[chosen.lrt1==1];
down.regulated.lrt1<-rownames(chosen.lrt1)[chosen.lrt1==-1];

pheatmap(t(log.cpm[((chosen.lrt1!=0) & (chosen.qlf1!=0)),(conditions$group=="CTRL.2h")|(conditions$group=="CTRL.8h")]),
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "Spectral")))(100), main = "Heatmap de expresión Control 2h vs Control 8h ",
         scale="column", show_colnames = F, cutree_rows = 2, clustering_method = "ward.D2",clustering_distance_rows = "manhattan")

venn.diagram(x=list(QLM=up.regulated.qlf, LRT=up.regulated.lrt1), 
             fill=c("Gray", "Cornflowerblue"),	main = "Upregulated Genes",cat.col=c("Black","DarkBlue"),filename="DEG_Upregulated_LRTvsQLF.tiff");
venn.diagram(x=list(QLM=down.regulated.qlf, LRT=down.regulated.lrt1),
             fill=c("Red", "Yellow"),	main = "Downregulated Genes",cat.col=c("DarkRed","Orange"), filename="DEG_Downregulated_LRTvsQLF.tiff");

int.up <- intersect(up.regulated.lrt1,up.regulated.qlf)
int.down <- intersect(down.regulated.lrt1,down.regulated.qlf)
write.csv(int.up, file = "intersect_upregulated.csv")
write.csv(int.down, file = "intersect_downregulated.csv")

full.qlf.table <- topTags(glmTest1, n="inf")$table
full.gene.list <- c(int.down,int.up)
write.csv(full.qlf.table, file = "QLF_table_CTRL2-8h.csv")
write.csv(full.gene.list, file = "CTRL2-8h_list.csv")

glmTest2<-glmQLFTest(QL.fit, contrast=as.numeric(contrastes["T4kV.2h.vs.T4kV.8h",]), poisson.bound=FALSE);
lrtTest2<-glmLRT(GLM.fit,contrast=as.numeric(contrastes["T4kV.2h.vs.T4kV.8h",]));

chosen.qlf2 <-decideTests(glmTest2,adjust.method="BH", p.value=0.05,lfc=log(thr,2));
chosen.lrt2<-decideTests(lrtTest2,adjust.method="BH", p.value=0.05,lfc=log(thr,2));

pheatmap(t(log.cpm[((chosen.lrt2!=0) & (chosen.qlf2!=0)),(conditions$group=="T4kV.2h")|(conditions$group=="T4kV.8h")]),
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "Spectral")))(100), main = "Heatmap de expresión tratamiento 4kV 2h vs 8h ", 
         scale="column", show_colnames = F, cutree_rows = 2, clustering_method = "ward.D2",clustering_distance_rows = "manhattan")


up.regulated.qlf2<-rownames(chosen.qlf2)[chosen.qlf2==1];
down.regulated.qlf2<-rownames(chosen.qlf2)[chosen.qlf2==-1];

up.regulated.lrt2<-rownames(chosen.lrt2)[chosen.lrt2==1];
down.regulated.lrt2<-rownames(chosen.lrt2)[chosen.lrt2==-1];

int.up1 <- intersect(up.regulated.lrt2,up.regulated.qlf2)
int.down1 <- intersect(down.regulated.lrt2,down.regulated.qlf2)
write.csv(int.up1, file = "intersect_upregulated_2.csv")
write.csv(int.down1, file = "intersect_downregulated_2.csv")

full.gene.list1 <- c(int.down1,int.up1)

write.csv(full.gene.list1, file = "T4kv_gene_list.csv")
full.qlf.table1 <- topTags(glmTest2, n="inf")$table
write.csv(full.qlf.table1, file = "T4kV_full_table.csv")

plotMD(glmTest1, main = "Expresión diferencial de Control 2h vs 8h",
       pch = 19, hl.col = c("purple", "orange"))
plotMD(glmTest2, main = "Expresión diferencial de Tratamiento 4kV 2h vs 8h", 
       pch = 19, hl.col = c("purple","orange"))

plot(full.qlf.table$logFC,-1*log(full.qlf.table$FDR,10), ylab="-log10(FDR)",xlab="log2(Fold Change)", 
       pch = 19, cex=0.3, main = "Volcano plot de genes diferencialmente expresados CTRL 2h vs CTRL 8h");
points(full.qlf.table[up.regulated.qlf,"logFC"],-1*log(full.qlf.table[up.regulated.qlf,"FDR"],10), 
       pch = 19, cex=0.3,col="hotpink");
points(full.qlf.table[down.regulated.qlf,"logFC"],-1*log(full.qlf.table[down.regulated.qlf,"FDR"],10), 
       pch = 19, cex=0.3,col="purple");
plot(full.qlf.table1$logFC,-1*log(full.qlf.table1$FDR,10), ylab="-log10(FDR)",xlab="log2(Fold Change)",
       pch = 19, cex=0.3, main = "Volcano plot de genes diferencialmente expresados T4kV 2h vs T4kV 8h");
points(full.qlf.table1[up.regulated.qlf2,"logFC"],-1*log(full.qlf.table1[up.regulated.qlf2,"FDR"],10), 
       pch = 19, cex=0.3,col="hotpink");
points(full.qlf.table1[down.regulated.qlf2,"logFC"],-1*log(full.qlf.table1[down.regulated.qlf2,"FDR"],10), 
       pch = 19, cex=0.3,col="purple");

