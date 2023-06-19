## De .csv a phyloseq
TAXA <- read_delim("PF_GF/domain;phylum;class;order;family;ge.txt", 
                   +     delim = ";", escape_double = FALSE, trim_ws = TRUE)
o<-sequence(3759)
t<-rep("otu", 3759) 
otu<-paste(t,o, sep = "_")
otu
taxa<-cbind(TAXA,otu) #TAXA es la tabla que envié antes solo que le pusi ese nombre

##hola esto es una prueba
otu_mat <- as.matrix(taxa) # convitiendo tabla en taxa
TAX = tax_table(otu_mat)

otus<-read.csv("PF_GF/raw_data/tabla2.csv")
otus1<-otus[,-1] ## quitando primera COLUMNA
otu_mat <- as.matrix(otus1)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)# tomará cada fila como un otu
muestras<-data.frame(c(Sample, samplenames))

samples = sample_names(sample)

library(phyloseq)
datos <- phyloseq(OTU, TAX, samples)
plot_bar(datos, fill = "DOMAIN")
