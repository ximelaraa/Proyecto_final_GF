## De .csv a phyloseq
TAXA <- read_delim("PF_GF/domain;phylum;class;order;family;ge.txt", 
                   +     delim = ";", escape_double = FALSE, trim_ws = TRUE)
o<-sequence(3759)
t<-rep("otu", 3759) 
otu<-paste(t,o, sep = "_")
otu
taxa<-cbind(TAXA,otu) #TAXA es la tabla que envié antes solo que le pusi ese nombre

##hola esto es una prueba
otu_mat <- as.matrix(taxa)
