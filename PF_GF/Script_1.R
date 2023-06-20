## De .csv a phyloseq

library(phyloseq)
library(readxl)
library(ggplot2)
TAXA <- read_delim("PF_GF/domain;phylum;class;order;family;ge.txt", 
        +   delim = ";", escape_double = FALSE, trim_ws = TRUE)
o<-sequence(3759)
t<-rep("otu", 3759) 
otu<-paste(t,o, sep = "_")
otu
taxa<-cbind(TAXA,otu) #TAXA es la tabla que envié antes solo que le pusi ese nombre

otu_mat <- as.matrix(taxa) # convitiendo tabla en taxa
TAX = tax_table(otu_mat)

otus<-read.csv("PF_GF/raw_data/tabla2.csv")

otus1<-otus[,-1] ## quitando primera COLUMNA
otu_mat <- as.matrix(otus1)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)# tomará cada fila como un otu
samples_df <- read_excel("PF_GF/raw_data/samplenames1.xlsx")


samples_df <- samples_df %>% 
  tibble::column_to_rownames("Sample")
samples = sample_data(samples_df)


datos <- phyloseq(OTU, TAX,samples)
datos


plot_bar(datos, fill = "DOMAIN")

plot_bar(datos, fill="subject", facet_grid=~DOMAIN)
plot_bar(datos, fill="zona", facet_grid=~DOMAIN)
plot_bar(datos, fill="tratamiento", facet_grid=~DOMAIN)
plot_bar(datos, fill="dia", facet_grid=~DOMAIN)
plot_bar(datos,"dia", facet_grid =~DOMAIN)
plot_bar(datos,"dia",fill="tratamiento", facet_grid =~DOMAIN)
plot_bar(datos,"tratamiento", fill="subject", facet_grid =~DOMAIN)

alpha_meas = c("Shannon", "Simpson", "InvSimpson")
(p <- plot_richness(datos,"dia","tratamiento",measures = alpha_meas))
