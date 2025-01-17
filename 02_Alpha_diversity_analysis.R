# Alpha diversity --------------------------------------------------------------

## Calculate diversity indices -------------------------------------------------
# import normalized data

freq_data <- read.csv("outputs/Frequency_table.csv", row.names = 1)

## Shannon
data_shannon <- diversity(t(freq_data), index = "shannon")
head(data_shannon)
write.table(data_shannon,
            "outputs/Alpha_diversity/alpha_diversity.xlsx",
            sheetName = "Shannon",
            col.names = TRUE,
            row.names = TRUE,
            append = TRUE)

## Inverse Simpson
data_invsimpson <- diversity(t(freq_data), index = "invsimpson")
head(data_invsimpson)
write.table(data_invsimpson,
            "outputs/Alpha_diversity/alpha_diversity.xlsx",
            sheetName = "Inverse_Simpson",
            col.names = TRUE,
            row.names = TRUE,
            append = TRUE)

## ACE
AceIndex <- estimateR(t(freq_data))
write.table(AceIndex,
            file = "outputs/Alpha_diversity/alpha_diversity.xlsx",
            sheetName = "ACE",
            col.names = TRUE,
            row.names = TRUE,
            append = TRUE)

## PD
library(phytools)
library("phyloseq")

phyl_tree <- ape::read.nexus("outputs/Hoverflies_microbiome.nex")
phyl_tree_rooted <- midpoint.root(phyl_tree)

Faith <- pd(t(freq_data), phyl_tree_rooted)
write.table(Faith, "outputs/Alpha_diversity/alpha_diversity.xlsx",
            sheetName = "PD", col.names = TRUE, row.names = TRUE, append = TRUE)

# Combine diversity indices in one file:
# outputs/Alpha_diversity/Diversity_indices_hoverflies.csv

# First look:
library(ggpubr)
divdat <- read.csv("outputs/Alpha_diversity/Diversity_indices_hoverflies.csv")

ggboxplot(divdat,
          x = "Species",
          y = "ACE",
          color = "black",
          fill = "Species")
ggboxplot(divdat,
          x = "Species",
          y = "Shannon",
          color = "black",
          fill = "Species")
ggboxplot(divdat,
          x = "Species",
          y = "inverse_simpson",
          color = "black",
          fill = "Species")
ggboxplot(divdat,
          x = "Species",
          y = "PD",
          color = "black",
          fill = "Species")


# Linear mixed models ----------------------------------------------------------
#compare differences in alpha diversity within and between species


# Toxomerus x Apis x Paragus----------------------------------------------------
divdat_PAT <- read.csv("outputs/Alpha_diversity/Diversity_indices_hoverflies_PAT.csv")

## Shannon ---------------------------------------------------------------------
# go from complex model to simplified model, drop all non-significant
# interactions and random factor
lm.H.PAT <- lmer(Shannon ~ Species +
                   Management +
                   Species:Management +
                   (1 | Field),
                 data = divdat_PAT)

# isSingular, remove random factor and switch to lm
lm.H.PAT <- lm(Shannon ~ Species +
                 Management +
                 Species:Management,
               data = divdat_PAT)

summary(lm.H.PAT)

# simplify model: remove insignificant interactions.
####### final model #######
lm.H.PAT <- lm(Shannon ~ Species +
                 Management,
               data = divdat_PAT)


# where do we see significance?
summary(lm.H.PAT)

# significances between the three "species" groups (male and female P.
# borbonicus and A. mellifera)
emmeans(lm.H.PAT, pairwise ~ Species, lmer.df = "Satterthwaite")

# which direction is the difference between species?
emmeans(lm.H.PAT, ~Species)



## inverse simpson -------------------------------------------------------------
# go from complex model to simplified model, drop all non-significant
# interactions and random factor
lm.IS.PAT <- lmer(inverse_simpson ~ Species +
                    Management +
                    Species:Management +
                    (1 | Field),
                  data = divdat_PAT)

# isSingular, remove random factor and switch to lm
lm.IS.PAT <- lm(inverse_simpson ~ Species +
                  Management +
                  Species:Management,
                data = divdat_PAT)

summary(lm.IS.PAT)

# simplify model: remove insignificant interactions.
####### final model #######
lm.IS.PAT <- lm(inverse_simpson ~ Species +
                  Management,
                data = divdat_PAT)

# where do we see significance?
summary(lm.IS.PAT)

# significances between the three "species" groups (male and female P.
# borbonicus and A. mellifera)
emmeans(lm.IS.PAT, pairwise ~ Species, lmer.df = "Satterthwaite")

# which direction is the difference between species?
emmeans(lm.IS.PAT, ~Species)



## ACE -------------------------------------------------------------------------
# go from complex model to simplified model, drop all non-significant
# interactions and random factor
lm.ACE.PAT <- lmer(S.ACE ~ Species +
                     Management +
                     Species:Management +
                     (1 | Field),
                   data = divdat_PAT)


#error: Model may not have converged with 1 eigenvalue close to zero: 4.3e-14
#standardize data:
divdat_PAT_st <- divdat_PAT %>% mutate_at(c("S.ACE"), ~(scale(.) %>% as.vector))

lm.ACE.PAT <- lmer(S.ACE ~ Species +
                     Management +
                     Species:Management +
                     (1 | Field),
                   data = divdat_PAT_st)

# isSingular, remove random factor and switch to lm
lm.ACE.PAT <- lm(S.ACE ~ Species +
                   Management +
                   Species:Management,
                 data = divdat_PAT_st)

summary(lm.ACE.PAT)

# simplify model: remove insignificant interactions.
####### final model #######
lm.ACE.PAT <- lm(S.ACE ~ Species +
                   Management,
                 data = divdat_PAT_st)

# where do we see significance?
summary(lm.ACE.PAT)

# significances between the three "species" groups (male and female P.
# borbonicus and A. mellifera)
emmeans(lm.ACE.PAT, pairwise ~ Species, lmer.df = "Satterthwaite")

# which direction is the difference between species?
emmeans(lm.ACE.PAT, ~Species)




## PD --------------------------------------------------------------------------
# go from complex model to simplified model, drop all non-significant
# interactions and random factor
lm.PD.PAT <- lmer(PD ~ Species +
                    Management +
                    Species:Management +
                    (1 | Field),
                  data = divdat_PAT)

# isSingular, remove random factor and switch to lm
lm.PD.PAT <- lm(PD ~ Species +
                  Management +
                  Species:Management,
                data = divdat_PAT)

summary(lm.PD.PAT)

# simplify model: remove insignificant interactions.
####### final model #######
lm.PD.PAT <- lm(PD ~ Species +
                  Management,
                data = divdat_PAT)

# where do we see significance?
summary(lm.PD.PAT)

# significances between the three "species" groups (male and female P.
# borbonicus and A. mellifera)
emmeans(lm.PD.PAT, pairwise ~ Species, lmer.df = "Satterthwaite")

# which direction is the difference between species?
emmeans(lm.PD.PAT, ~Species)




# Ischiodon x Apis x Paragus ---------------------------------------------------
divdat_PAIS <- read.csv("outputs/Alpha_diversity/Diversity_indices_hoverflies_PAIS.csv")


## Shannon ---------------------------------------------------------------------
#go from complex model to simplified model, drop all non-significant
# interactions and random factor
####### final model #######
lm.H.PAIS <- lmer(Shannon ~ Species +
                    (1 | Field),
                  data = divdat_PAIS)

summary(lm.H.PAIS)

# where do we see significance?
summary(lm.H.PAIS)

# significances between the three "species" groups (male and female P.
# borbonicus and A. mellifera)
emmeans(lm.H.PAIS, pairwise ~ Species, lmer.df = "Satterthwaite")

# which direction is the difference between species?
emmeans(lm.H.PAIS, ~Species)


## Inverse Simpson -------------------------------------------------------------

# go from complex model to simplified model, drop all non-significant
# interactions and random factor
####### final model #######
lm.IS.PAIS <- lmer(inverse_simpson ~ Species +
                     (1 | Field),
                   data = divdat_PAIS)

summary(lm.IS.PAIS)

# where do we see significance?
summary(lm.IS.PAIS)

# significances between the three "species" groups (male and female P.
# borbonicus and A. mellifera)
emmeans(lm.IS.PAIS, pairwise ~ Species, lmer.df = "Satterthwaite")

# which direction is the difference between species?
emmeans(lm.IS.PAIS, ~Species)


## ACE -------------------------------------------------------------------------


# go from complex model to simplified model, drop all non-significant
# interactions and random factor

lm.ACE.PAIS <- lmer(S.ACE ~ Species +
                      (1 | Field),
                    data = divdat_PAIS)

#error: Model may not have converged with 1 eigenvalue close to zero: 1.3e-14
#standardize data:
divdat_PAIS_st <- divdat_PAIS %>% mutate_at(c("S.ACE"), ~(scale(.) %>% as.vector))

####### final model #######
lm.ACE.PAIS <- lmer(S.ACE ~ Species +
                      (1 | Field),
                    data = divdat_PAIS_st)


# where do we see significance?
summary(lm.ACE.PAIS)

# significances between the three "species" groups (male and female P.
# borbonicus and A. mellifera)
emmeans(lm.ACE.PAIS, pairwise ~ Species, lmer.df = "Satterthwaite")

# which direction is the difference between species?
emmeans(lm.ACE.PAIS, ~Species)


## PD --------------------------------------------------------------------------

# go from complex model to simplified model, drop all non-significant
# interactions and random factor
####### final model #######
lm.PD.PAIS <- lmer(PD ~ Species +
                     (1 | Field),
                   data = divdat_PAIS)

# isSingular, remove random factor and switch to lm
lm.PD.PAIS <- lm(PD ~ Species,
                 data = divdat_PAIS)


# where do we see significance?
summary(lm.PD.PAIS)

# significances between the three "species" groups (male and female P.
# borbonicus and A. mellifera)
emmeans(lm.PD.PAIS, pairwise ~ Species, lmer.df = "Satterthwaite")

# which direction is the difference between species?
emmeans(lm.PD.PAIS, ~Species)
