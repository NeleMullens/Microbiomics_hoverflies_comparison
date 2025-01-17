# ALDEx2 analysis -----------------------

### Aldex honeybees ---------------------
Genera <- read.csv("outputs/Aldex/Genera/Aldex_Genera_Honeybees.csv", row.names=1) 

Reduce_genera <- Genera %>% 
  group_by(Genus) %>% 
  summarise(across(everything(), ~ sum(.x)))

write.xlsx(Reduce_genera, file = "outputs/Aldex/Genera/Input_Aldex_genera_honeybees.xlsx" , sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE)

aldexmeta <- read.csv("outputs/Aldex2/Genera/Metadata_apis.csv", row.names = 1)
filt_data <- read.csv("outputs/Aldex2/Genera/Input_genera_Apis.csv", row.names = 1)

# can only compare 2 conditions at the time, not for all combined
Aldex_analysis_ma <- aldex(filt_data, conditions = aldexmeta$Management, mc.samples = 1000, test="t", effect=TRUE, denom="all", verbose=TRUE)
Aldex_analysis_al <- aldex(filt_data, conditions = aldexmeta$Altitude, mc.samples = 1000, test="t", effect=TRUE, denom="all", verbose=TRUE)

write.xlsx(Aldex_analysis_ma, file = "outputs/Aldex/Genera/Output_Aldex_Genera_A.xlsx" , sheetName = "Management", 
           col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx(Aldex_analysis_al, file = "outputs/Aldex/Genera/Output_Aldex_Genera_A.xlsx" , sheetName = "Altitude", 
           col.names = TRUE, row.names = TRUE, append = TRUE)


##################
###### TEST ######
##################
# repeat, but with ASV's
aldexmeta <- read.csv("outputs/Aldex2/Metadata_Apis.csv", row.names = 1)
filt_data <- read.csv("outputs/Aldex2/ASV/Complete_filtered_data_Apis.csv", row.names = 1)

# can only compare 2 conditions at the time, not for all combined
Aldex_analysis_ma <- aldex(filt_data, conditions = aldexmeta$Management, mc.samples = 1000, test="t", effect=TRUE, denom="all", verbose=TRUE)
Aldex_analysis_al <- aldex(filt_data, conditions = aldexmeta$Altitude, mc.samples = 1000, test="t", effect=TRUE, denom="all", verbose=TRUE)

write.xlsx(Aldex_analysis_ma, file = "outputs/Aldex2/ASV/Output_Aldex_ASV_A.xlsx" , sheetName = "Management", 
           col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx(Aldex_analysis_al, file = "outputs/Aldex2/ASV/Output_Aldex_ASV_A.xlsx" , sheetName = "Altitude", 
           col.names = TRUE, row.names = TRUE, append = TRUE)

# paragus 

aldexmeta <- read.csv("outputs/Aldex2/Metadata_Paragus.csv", row.names = 1)
filt_data <- read.csv("outputs/Aldex2/ASV/Complete_filtered_data_Paragus.csv", row.names = 1)

# can only compare 2 conditions at the time, not for all combined
Aldex_analysis_ma <- aldex(filt_data, conditions = aldexmeta$Management, mc.samples = 1000, test="t", effect=TRUE, denom="all", verbose=TRUE)
Aldex_analysis_al <- aldex(filt_data, conditions = aldexmeta$Altitude, mc.samples = 1000, test="t", effect=TRUE, denom="all", verbose=TRUE)

write.xlsx(Aldex_analysis_ma, file = "outputs/Aldex2/ASV/Output_Aldex_ASV_P.xlsx" , sheetName = "Management", 
           col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx(Aldex_analysis_al, file = "outputs/Aldex2/ASV/Output_Aldex_ASV_P.xlsx" , sheetName = "Altitude", 
           col.names = TRUE, row.names = TRUE, append = TRUE)


# repeat genera but with the different sexes. 
aldexmeta <- read.csv("outputs/Aldex2/Metadata_Paragus_m.csv", row.names = 1)
filt_data <- read.csv("outputs/Aldex2/Genera/Input_genera_male_Paragus.csv", row.names = 1)

# can only compare 2 conditions at the time, not for all combined
#Aldex_analysis_ma <- aldex(filt_data, conditions = aldexmeta$Management, mc.samples = 1000, test="t", effect=TRUE, denom="all", verbose=TRUE)
Aldex_analysis_al <- aldex(filt_data, conditions = aldexmeta$Altitude, mc.samples = 1000, test="t", effect=TRUE, denom="all", verbose=TRUE)

write.xlsx(Aldex_analysis_al, file = "outputs/Aldex2/Genera/Output_Aldex_SexAltitude_P.xlsx" , sheetName = "Males", 
           col.names = TRUE, row.names = TRUE, append = TRUE)


aldexmeta <- read.csv("outputs/Aldex2/Metadata_Paragus_f.csv", row.names = 1)
filt_data <- read.csv("outputs/Aldex2/Genera/Input_genera_female_Paragus.csv", row.names = 1)

# can only compare 2 conditions at the time, not for all combined
#Aldex_analysis_ma <- aldex(filt_data, conditions = aldexmeta$Management, mc.samples = 1000, test="t", effect=TRUE, denom="all", verbose=TRUE)
Aldex_analysis_al <- aldex(filt_data, conditions = aldexmeta$Altitude, mc.samples = 1000, test="t", effect=TRUE, denom="all", verbose=TRUE)

write.xlsx(Aldex_analysis_al, file = "outputs/Aldex2/Genera/Output_Aldex_SexAltitude_P.xlsx" , sheetName = "Females", 
           col.names = TRUE, row.names = TRUE, append = TRUE)

### Aldex paragus ---------------------
Genera <- read.csv("outputs/Aldex/Genera/Aldex_Phyla_Paragus.csv", row.names=1) 

Reduce_genera <- Genera %>% 
  group_by(Genus) %>% 
  summarise(across(everything(), ~ sum(.x)))

write.xlsx(Reduce_genera, file = "outputs/Aldex/Genera/Input_Aldex_genera_Paragus.xlsx" , sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE)

aldexmeta <- read.csv("outputs/Aldex/Genera/Metadata_Paragus.csv", row.names = 1)
filt_data <- read.csv("outputs/Aldex/Genera/Input_Aldex_genera_Paragus.csv", row.names = 1)

Aldex_analysis_sex <- aldex(filt_data, conditions = aldexmeta$Sex, mc.samples = 1000, test="t", effect=TRUE, denom="all", verbose=TRUE)
Aldex_analysis_ma <- aldex(filt_data, conditions = aldexmeta$Management, mc.samples = 1000, test="t", effect=TRUE, denom="all", verbose=TRUE)
Aldex_analysis_al <- aldex(filt_data, conditions = aldexmeta$Altitude, mc.samples = 1000, test="t", effect=TRUE, denom="all", verbose=TRUE)

write.xlsx(Aldex_analysis_ma, file = "outputs/Aldex/Genera/Output_Aldex_Genera_P.xlsx" , sheetName = "Management", 
           col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx(Aldex_analysis_al, file = "outputs/Aldex/Genera/Output_Aldex_Genera_P.xlsx" , sheetName = "Altitude", 
           col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(Aldex_analysis_sex, file = "outputs/Aldex/Genera/Output_Aldex_Genera_P.xlsx" , sheetName = "Sex", 
           col.names = TRUE, row.names = TRUE, append = TRUE)



### Aldex toxomerus ---------------------
Genera <- read.csv("outputs/Aldex/Genera/Aldex_Phyla_Toxomerus.csv", row.names=1) 

Reduce_genera <- Genera %>% 
  group_by(Genus) %>% 
  summarise(across(everything(), ~ sum(.x)))

write.xlsx(Reduce_genera, file = "outputs/Aldex/Genera/Input_Aldex_genera_Toxomerus.xlsx" , sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE)

aldexmeta <- read.csv("outputs/Aldex/Genera/Metadata_Toxomerus.csv", row.names = 1)
filt_data <- read.csv("outputs/Aldex/Genera/Input_Aldex_genera_Toxomerus.csv", row.names = 1)

Aldex_analysis_sex <- aldex(filt_data, conditions = aldexmeta$Sex, mc.samples = 1000, test="t", effect=TRUE, denom="all", verbose=TRUE)
Aldex_analysis_ma <- aldex(filt_data, conditions = aldexmeta$Management, mc.samples = 1000, test="t", effect=TRUE, denom="all", verbose=TRUE)

write.xlsx(Aldex_analysis_ma, file = "outputs/Aldex/Genera/Output_Aldex_Genera_T.xlsx" , sheetName = "Management", 
           col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx(Aldex_analysis_sex, file = "outputs/Aldex/Genera/Output_Aldex_Genera_T.xlsx" , sheetName = "Sex", 
           col.names = TRUE, row.names = TRUE, append = TRUE)



### Aldex Ischiodon ---------------------
Genera <- read.csv("outputs/Aldex/Genera/Aldex_Phyla_Ischiodon.csv", row.names=1) 

Reduce_genera <- Genera %>% 
  group_by(Genus) %>% 
  summarise(across(everything(), ~ sum(.x)))

write.xlsx(Reduce_genera, file = "outputs/Aldex/Genera/Input_Aldex_genera_Ischiodon.xlsx" , sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE)

aldexmeta <- read.csv("outputs/Aldex/Genera/Metadata_Ischiodon.csv", row.names = 1)
filt_data <- read.csv("outputs/Aldex/Genera/Input_Aldex_genera_Ischiodon.csv", row.names = 1)

Aldex_analysis_sex <- aldex(filt_data, conditions = aldexmeta$Sex, mc.samples = 1000, test="t", effect=TRUE, denom="all", verbose=TRUE)

write.xlsx(Aldex_analysis_sex, file = "outputs/Aldex/Genera/Output_Aldex_Genera_I.xlsx" , sheetName = "Sex", 
           col.names = TRUE, row.names = TRUE, append = FALSE)



### Aldex paragus honeybees ---------------------
Genera <- read.csv("outputs/Aldex/Genera/Aldex_Phyla_Honeybees_Paragus.csv", row.names=1) 

Reduce_genera <- Genera %>% 
  group_by(Genus) %>% 
  summarise(across(everything(), ~ sum(.x)))

write.xlsx(Reduce_genera, file = "outputs/Aldex/Genera/Input_Aldex_genera_Honeybees_Paragus.xlsx" , sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE)

aldexmeta <- read.csv("outputs/Aldex/Genera/Metadata_honeybees_Paragus.csv", row.names = 1)
filt_data <- read.csv("outputs/Aldex/Genera/Input_Aldex_genera_Honeybees_Paragus.csv", row.names = 1)

Aldex_analysis_sp <- aldex(filt_data, conditions = aldexmeta$Species, mc.samples = 1000, test="t", effect=TRUE, denom="all", verbose=TRUE)
Aldex_analysis_ma <- aldex(filt_data, conditions = aldexmeta$Management, mc.samples = 1000, test="t", effect=TRUE, denom="all", verbose=TRUE)
Aldex_analysis_al <- aldex(filt_data, conditions = aldexmeta$Altitude, mc.samples = 1000, test="t", effect=TRUE, denom="all", verbose=TRUE)

write.xlsx(Aldex_analysis_ma, file = "outputs/Aldex/Genera/Output_Aldex_Genera_P_A.xlsx" , sheetName = "Management", 
           col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx(Aldex_analysis_al, file = "outputs/Aldex/Genera/Output_Aldex_Genera_P_A.xlsx" , sheetName = "Altitude", 
           col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(Aldex_analysis_sp, file = "outputs/Aldex/Genera/Output_Aldex_Genera_P_A.xlsx" , sheetName = "Species", 
           col.names = TRUE, row.names = TRUE, append = TRUE)



### Aldex paragus toxomerus ---------------------
Genera <- read.csv("outputs/Aldex/Genera/Aldex_Genera_Paragus_Toxomerus.csv", row.names=1) 

Reduce_genera <- Genera %>% 
  group_by(Genus) %>% 
  summarise(across(everything(), ~ sum(.x)))

write.xlsx(Reduce_genera, file = "outputs/Aldex/Genera/Input_Aldex_genera_Paragus_Toxomerus.xlsx" , sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE)

aldexmeta <- read.csv("outputs/Aldex/Genera/Metadata_Paragus_toxomerus.csv", row.names = 1)
filt_data <- read.csv("outputs/Aldex/Genera/Input_Aldex_genera_Paragus_Toxomerus.csv", row.names = 1)


Aldex_analysis_sp <- aldex(filt_data, conditions = aldexmeta$Species, mc.samples = 1000, test="t", effect=TRUE, denom="all", verbose=TRUE)
Aldex_analysis_ma <- aldex(filt_data, conditions = aldexmeta$Management, mc.samples = 1000, test="t", effect=TRUE, denom="all", verbose=TRUE)
Aldex_analysis_sex <- aldex(filt_data, conditions = aldexmeta$Sex, mc.samples = 1000, test="t", effect=TRUE, denom="all", verbose=TRUE)

write.xlsx(Aldex_analysis_ma, file = "outputs/Aldex/Genera/Output_Aldex_Genera_P_T.xlsx" , sheetName = "Management", 
           col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx(Aldex_analysis_sp, file = "outputs/Aldex/Genera/Output_Aldex_Genera_P_T.xlsx" , sheetName = "Species", 
           col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(Aldex_analysis_sex, file = "outputs/Aldex/Genera/Output_Aldex_Genera_P_T.xlsx" , sheetName = "Sex", 
           col.names = TRUE, row.names = TRUE, append = TRUE)



### Aldex paragus ischiodon ---------------------
Genera <- read.csv("outputs/Aldex/Genera/Aldex_Genera_Paragus_Ischiodon.csv", row.names=1) 

Reduce_genera <- Genera %>% 
  group_by(Genus) %>% 
  summarise(across(everything(), ~ sum(.x)))

write.xlsx(Reduce_genera, file = "outputs/Aldex/Genera/Input_Aldex_genera_Paragus_Ischiodon.xlsx" , sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE)

aldexmeta <- read.csv("outputs/Aldex/Genera/Metadata_Paragus_Ischiodon.csv", row.names = 1)
filt_data <- read.csv("outputs/Aldex/Genera/Input_Aldex_genera_Paragus_Ischiodon.csv", row.names = 1)

Aldex_analysis_sp <- aldex(filt_data, conditions = aldexmeta$Species, mc.samples = 1000, test="t", effect=TRUE, denom="all", verbose=TRUE)
Aldex_analysis_sex <- aldex(filt_data, conditions = aldexmeta$Sex, mc.samples = 1000, test="t", effect=TRUE, denom="all", verbose=TRUE)

write.xlsx(Aldex_analysis_sp, file = "outputs/Aldex/Genera/Output_Aldex_Genera_P_I.xlsx" , sheetName = "Species", 
           col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx(Aldex_analysis_sex, file = "outputs/Aldex/Genera/Output_Aldex_Genera_P_I.xlsx" , sheetName = "Sex", 
           col.names = TRUE, row.names = TRUE, append = TRUE)