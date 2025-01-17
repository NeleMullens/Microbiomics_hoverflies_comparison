#Figure 1 ----------------------------------------------------------------------

#QGIS map

#Figure 2 a --------------------------------------------------------------------
#boxplot I. aegyptius Shannon females vs males

divdat <- read.csv("outputs/Alpha_diversity/Diversity_indices_hoverflies_Ischiodon.csv")
my_comparisons = list( c("female", "male"))

# H <- ggboxplot(
#   divdat, 
#   x = "Sex", 
#   y= "Shannon", 
#   ylab = "Shannon",
#   fill = "Sex",
#   xlab = "",  
#   legend = "none",
#   # palette = c("red", "blue")) +
# palette = c("brown2", "blueviolet")) +
  # scale_fill_manual(values = c("red", "blue")) +
  # scale_fill_brewer(palette="Dark2") +
#   geom_jitter(width = 0.3, height = 0, size=3, shape = 1) +
#   theme(text = element_text(size = 15)) +
#   stat_compare_means(comparisons = my_comparisons, label.y = c(5.5),tip.length=0.0,
#                      symnum.args = list(cutpoints = c(0.0001, 0.01, 0.05, 1),
#                                         symbols = c("**", "*", "ns")))


H <- ggboxplot(
  divdat, 
  x = "Sex", 
  y= "Shannon", 
  ylab = "Shannon",
  #fill = "Sex",
  xlab = "",  
  legend = "none",
  color = "Sex", palette =c("#E7B800", "#00AFBB"), #"#FC4E07", "#00AFBB", "#E7B800" 
  add = "jitter", size = 0.8) + #, shape = "Sex"
  theme(text = element_text(size = 15)) +
  stat_compare_means(comparisons = my_comparisons, label.y = c(5.5),tip.length=0.0,
                     symnum.args = list(cutpoints = c(0.0001, 0.01, 0.05, 1),
                                        symbols = c("**", "*", "ns")))

H

#width = 0.2, height = 0, alpha = 0.7,


#Figure 2 b --------------------------------------------------------------------
#boxplot I. aegyptius Shannon females vs males

IS <- ggboxplot(
  divdat, 
  x = "Sex", 
  y= "inverse_simpson", 
  ylab = "Inverse Simpson",
  #fill = "Sex",
  xlab = "",
  legend = "none", 
  #palette = c("brown2", "blueviolet")) +
  color = "Sex", palette =c("#E7B800", "#00AFBB"), #"#FC4E07", "#00AFBB", "#E7B800" 
  add = "jitter", size = 0.8) + #, shape = "Sex"
  theme(text = element_text(size = 15)) +
  stat_compare_means(comparisons = my_comparisons, label.y = c(95, 95),tip.length=0.0,
                     symnum.args = list(cutpoints = c(0.0001, 0.01, 0.05, 1),
                                        symbols = c("***", "*", "ns")))
IS

ggarrange(H, IS, labels = c("A", "B"), ncol = 2, nrow = 1, widths = c(30, 35))

#Figure 3 --------------------------------------------------------------------
#PCoA I. aegyptius Shannon females vs males

PCA_I <- read.csv("outputs/PCA_Ischiodon_sex.csv", row.names = 1)
colnames(PCA_I) <- c("pcoa1", "pcoa2", "Sex")
PCA_I %>% as_tibble(rownames="samples")

pca_sex <- PCA_I %>% 
  ggplot(aes(x=pcoa1, y=pcoa2, color = Sex)) +
  geom_point(aes(col = Sex)) +
  geom_point(size=3) +
  labs(x="PCo 1 (11.8%)", y="PCo 2 (8.2%)") +
  #scale_color_manual(values = c("brown2", "blueviolet")) +
  scale_color_manual(values = c("#E7B800", "#00AFBB")) +
  stat_ellipse(level = 0.95) +
  theme_minimal() +
  theme(axis.title=element_text(size=15),
        axis.text=element_text(size=12),
        legend.title = element_text(size = 15, face = "bold"),
        legend.text=element_text(size=15)) #+
  #theme_bw()

pca_sex

#theme(legend.position = "top") +
#scale_color_manual(values = c("darkolivegreen4", "darkgoldenrod3"))  


#ggarrange(pca_sex, labels = c("C"), ncol = 1, nrow = 1)

ggarrange(H, IS, pca_sex, labels = c("A", "B", "C"), ncol = 3, nrow = 1, widths = c(20, 20, 40))

ggsave("outputs/Figure_2_thicker_lines_pollinators.tiff", device='tiff', dpi=300)


#Figure 3 ----------------------------------------------------------------------

library(microeco)

otu_table <- read.csv("Outputs/Discription_microbiome/Ischiodon/Discription_microbiome_input_Ischiodon.csv", row.names = 1)
taxonomy_table <- read.csv("outputs/Discription_microbiome/Ischiodon/Discription_microbiome_taxonomy_Ischiodon.csv", row.names=1)
sample_info <- read.csv("outputs/Discription_microbiome/Ischiodon/Discription_microbiome_metadata_Ischiodon.csv", row.names = 1)


# create an object of microtable
dataset <- microtable$new(sample_table = sample_info, otu_table = otu_table, tax_table = taxonomy_table)

phy <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 8)
P1 <- phy$plot_bar(
  others_color = "grey70", 
  facet = "Sex", 
  xtext_keep = FALSE, 
  legend_text_italic = FALSE) +
  theme(strip.text = element_text(size = 14),
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 12))

P1

gen <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 8)
P2 <- gen$plot_bar(
  others_color = "grey70", 
  facet = "Sex", 
  xtext_keep = FALSE, 
  legend_text_italic = FALSE) +
  theme(strip.text = element_text(size = 14),
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 12))

P2

ggarrange(P1, P2, labels = c("A", "B"), ncol = 2, nrow = 1, widths = c(0.8, 1))
plot_grid(phy, gen, labels = c('A', 'B'), label_size = 12, ncol = 2)

ggsave("outputs/Figure_3_pollinators.tiff", device='tiff', dpi=300)

fam <- trans_abund$new(dataset = dataset, taxrank = "Family", ntaxa = 8)
P3 <- fam$plot_bar(
  others_color = "grey70", 
  facet = "Sex", 
  xtext_keep = FALSE, 
  legend_text_italic = FALSE) +
  theme(strip.text = element_text(size = 14),
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 12))

P3
