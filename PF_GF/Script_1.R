### De .csv a phyloseq ----

library(phyloseq)
library(readxl)
library(ggplot2)
library(tibble)
library(vegan)
TAXA <- read_delim("PF_GF/domain;phylum;class;order;family;ge.txt", #lee la tabla donde viene la taxonomía
        +   delim = ";", escape_double = FALSE, trim_ws = TRUE)
o<-sequence(3759) 
t<-rep("otu", 3759) 
otu<-paste(t,o, sep = "_") #creando la columna de número de otu
otu
taxa<-cbind(TAXA,otu) #añadiendo la comluna de otu a la tabla

otu_mat <- as.matrix(taxa) # convitiendo tabla en taxa
TAX = tax_table(otu_mat)


otus<-read.csv("PF_GF/raw_data/tabla2.csv")


otus1<-otus[,-1] ## quitando primera COLUMNA
otu_mat <- as.matrix(otus1)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)# tomará cada fila como un otu
samples_df <- read_excel("PF_GF/raw_data/samplenames1.xlsx")


samples_df <- samples_df %>% #definiendo los nombres de las filas de la columna otu
  tibble::column_to_rownames("Sample")
samples = sample_data(samples_df)


datos <- phyloseq(OTU, TAX,samples)
datos
 #visualización de datos
sample_names(datos)
rank_names(datos)
sample_variables(datos)

### Gráficos con el objeto phyloseq ----
plot_bar(datos, fill = "DOMAIN")+ #sample vs abundance por dominio
  geom_bar(aes(color=DOMAIN, fill=DOMAIN), stat="identity", position="stack")

plot_bar(datos, fill="subject", facet_grid=~DOMAIN) +
  geom_bar(aes(color=subject, fill=subject), stat="identity", position="stack")

plot_bar(datos, fill="zona", facet_grid=~DOMAIN) +
  geom_bar(aes(color=zona, fill=zona), stat="identity", position="stack")

plot_bar(datos, fill="tratamiento", facet_grid=~DOMAIN) +
  geom_bar(aes(color=tratamiento, fill=tratamiento), stat="identity", position="stack")

plot_bar(datos, fill="dia", facet_grid=~DOMAIN) +
  geom_bar(aes(color=dia, fill=dia), stat="identity", position="stack")

plot_bar(datos,"dia", facet_grid =~DOMAIN) +
  geom_bar(aes(color=DOMAIN), stat="identity", position="stack")

plot_bar(datos,"dia",fill="tratamiento", facet_grid =~DOMAIN) +
  geom_bar(aes(color=tratamiento, fill=tratamiento), stat="identity", position="stack")

plot_bar(datos,"tratamiento", fill="subject", facet_grid =~DOMAIN)

#estimadores de la diversidad alfa
alpha_meas = c("Shannon", "Simpson", "InvSimpson")
(p <- plot_richness(datos,"dia","tratamiento",measures = alpha_meas))


### Creación del filtrado que deje únicamente los datos más relevantes ----
tabla_abundancias <- otu_table(datos)
suma_abundancias <- rowSums(tabla_abundancias)
umbral <- 0.00001
filtro_taxa <- names(suma_abundancias[suma_abundancias > umbral])
datos_filtrados <- prune_taxa(filtro_taxa, datos)
datos_filtrados #objeto phyloseq con el filtrado

#normalización de los datos filtrados
samples_filtro <- sample_names(datos_filtrados)
nuevos_datos <- transform_sample_counts(datos_filtrados, function(x) x / sum(x))
sum(otu_table(nuevos_datos)[,1]) #Si la suma es igual a 1 entonces se ha normalizado


sample_names(nuevos_datos)
rank_names(nuevos_datos)
sample_variables(nuevos_datos)
### Gráficos con el filtrado ----
plot_bar(nuevos_datos, fill = "DOMAIN") + 
  geom_bar(aes(color=DOMAIN, fill=DOMAIN), stat="identity", position="stack")

plot_bar(nuevos_datos, fill="subject", facet_grid=~DOMAIN) +
  geom_bar(aes(color=subject, fill=subject), stat="identity", position="stack")

plot_bar(nuevos_datos, fill="zona", facet_grid=~DOMAIN) +
  geom_bar(aes(color=zona, fill=zona), stat="identity", position="stack")

plot_bar(nuevos_datos, fill="tratamiento", facet_grid=~DOMAIN) +
  geom_bar(aes(color=tratamiento, fill=tratamiento), stat="identity", position="stack")

plot_bar(nuevos_datos, fill="dia", facet_grid=~DOMAIN) +
  geom_bar(aes(color=dia, fill=dia), stat="identity", position="stack")

plot_bar(nuevos_datos,"dia",fill="tratamiento", facet_grid =~DOMAIN) +
  geom_bar(aes(color=tratamiento, fill=tratamiento), stat="identity", position="stack")

plot_bar(nuevos_datos,"tratamiento", fill="subject", facet_grid =~DOMAIN)
#estimadores de la diversidad alfa
alpha_meas = c("Shannon", "Simpson", "InvSimpson")
(p <- plot_richness(nuevos_datos,"dia","tratamiento",measures = alpha_meas))

#Gráficos exploratorios

barplot(sort(taxa_sums(datos), TRUE)/nsamples(datos), las=2) + title(main = "Todos los datos") 
barplot(sort(taxa_sums(datos), TRUE)[1:30]/nsamples(datos), las=2)  + title(main = "Primeros 30")

#Demasiados datos irrelevantes

barplot(sort(taxa_sums(nuevos_datos), TRUE)/nsamples(nuevos_datos), las=2) + title(main = "Todos los datos filtrados")
barplot(sort(taxa_sums(nuevos_datos), TRUE)[1:30]/nsamples(nuevos_datos), las=2)  + title(main = "Primeros 30")
plot_bar(nuevos_datos,"DOMAIN",fill="tratamiento" ,facet_grid =~dia)
plot_bar(nuevos_datos,"DOMAIN",fill="tratamiento" ,facet_grid =~subject)
plot_bar(nuevos_datos,"DOMAIN",fill="dia" ,facet_grid =~zona)


###diversidad beta metodo bray curtis
bray <- phyloseq::distance(datos_filtrados, method = "bray")
bray <- as.matrix(bray)
boxplot(bray)


ord = ordinate(datos_filtrados, method="PCoA", distance = "bray")

plot_ordination(datos_filtrados, ord, color = "dia", shape="subject") + 
  geom_point(size=4) + 
  stat_ellipse(aes(group=tratamiento))+
  scale_shape_manual(values = c(0,1,2,3,
                                5,6,7,8,9,10,11,12,13,14))

plot_ordination(datos_filtrados, ord, type="taxa", color="DOMAIN", 
               title="OTUs", label="GENUS") + 
  facet_wrap(~DOMAIN, 3)