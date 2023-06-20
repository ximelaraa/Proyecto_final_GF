### De .csv a phyloseq ----

library(phyloseq)
library(readxl)
library(ggplot2)
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
plot_bar(datos, fill = "DOMAIN")
plot_bar(datos, fill="subject", facet_grid=~DOMAIN)
plot_bar(datos, fill="zona", facet_grid=~DOMAIN)
plot_bar(datos, fill="tratamiento", facet_grid=~DOMAIN)
plot_bar(datos, fill="dia", facet_grid=~DOMAIN)
plot_bar(datos,"dia", facet_grid =~DOMAIN)
plot_bar(datos,"dia",fill="tratamiento", facet_grid =~DOMAIN)
plot_bar(datos,"tratamiento", fill="subject", facet_grid =~DOMAIN)
#estimadores de la diversidad alfa
alpha_meas = c("Shannon", "Simpson", "InvSimpson")
(p <- plot_richness(datos,"dia","tratamiento",measures = alpha_meas))

### Creación del filtrado que deje únicamente los datos más relevantes ----
tabla_abundancias <- otu_table(datos)
suma_abundancias <- rowSums(tabla_abundancias)
umbral <- 0.00001
filtro_taxa <- names(suma_abundancias[suma_abundancias < umbral])
datos_filtrados <- prune_taxa(filtro_taxa, datos)
datos_filtrados #objeto phyloseq con el filtrado 

### Gráficos con el filtrado ----
plot_bar(datos_filtrados, fill = "DOMAIN")
plot_bar(datos_filtrados, fill="subject", facet_grid=~DOMAIN)
plot_bar(datos_filtrados, fill="zona", facet_grid=~DOMAIN)
plot_bar(datos_filtrados, fill="tratamiento", facet_grid=~DOMAIN)
plot_bar(datos_filtrados, fill="dia", facet_grid=~DOMAIN)
plot_bar(datos_filtrados,"dia", facet_grid =~DOMAIN)
plot_bar(datos_filtrados,"dia",fill="tratamiento", facet_grid =~DOMAIN)
plot_bar(datos_filtrados,"tratamiento", fill="subject", facet_grid =~DOMAIN)
#estimadores de la diversidad alfa
alpha_meas = c("Shannon", "Simpson", "InvSimpson")
(p <- plot_richness(datos_filtrados,"dia","tratamiento",measures = alpha_meas))

