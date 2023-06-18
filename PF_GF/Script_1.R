## De .csv a phyloseq
o<-sequence(3759)
t<-rep("otu", 3759) 
otu<-paste(t,o, sep = "_")
otu
taxa<-cbind(TAXA,otu) #TAXA es la tabla que envié antes solo que le pusi ese nombre

##hola esto es una prueba