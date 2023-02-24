#PP analysis done for Bergholz et al Nature 2023#

#Input to this code is mcmicro output www.mcmicro.org

#Load packages
library(doBy)
library(tidyverse)
library(ggridges)
library(summarytools)
library(naivestates)
library(ggpubr)
library(ggsci)
library(miscTools)

if( !require(devtools) ) install.packages("devtools")
devtools::install_github( "labsyspharm/naivestates" )

#Set working dir

setwd('~/R/cycif_analysis/data/')

#Import mcmicro quantification files from /data 

fns <- c( PPA = "unmicst-PP_J6Feb19_PPA.csv",
          PPB = "unmicst-PP_J6Feb19_PPB.csv",
          PP = "unmicst-PP_J6Feb19_PP.csv")

# Load each file and concatenate into a single data frame
# Sample names will be appear in column SampleID
# Generate unique IDs by combining sampleIDs and cellIDs
# Apply log10 normalization

X <- map( fns, read_csv ) %>%
  bind_rows( .id = "SampleID" ) 

X <- X %>%
  set_names(~ str_replace_all(.,"_cellMask", "")) %>%
  nest( uniqueID = c(SampleID, CellID) )

# Vector with markers of interest names MUST match colnames for tibble

markers <- c('DAPI_1', 'PTEN', 'CD8', 'CD3e', 'CD4', 'CD45', 'CKc11', 'CKae', 'CD11b', 'F480', 'Vimentin',
             'pSTAT3', 'FOXP3', 'NFKB', 'cCaspase3', 'Ki67', 'IRF3', 'pAKT')

X <- X%>%
  mutate_at( markers, ~log10(.x+1) )

Y <- X %>% unnest(uniqueID)

Y%>%
  ggplot(aes(x = CD45, y = SampleID))+
  stat_density_ridges(quantile_lines = TRUE, quantiles = 0.5, bandwidth = 0.1, alpha = 0.5)+
  theme_ridges()+
  labs(x = 'CD45')

# Ask GMMfit to fit all markers
GMM <- GMMfit( X, uniqueID, markers )

#Retrieve marker data after fit

# Value  - original expression
# AdjVal - normalized expression
# Prob   - probability of marker expression

pSTAT3 <- GMM %>% filter(Marker == 'pSTAT3') %>% pluck('Values', 1) %>% unnest(uniqueID)

# Compare raw values across samples

ggplot(CD45, aes(x =  Value, y = SampleID))+
 stat_density_ridges(quantile_lines = TRUE, quantiles = 0.5, bandwidth = 0.01, scale = 0.9, alpha = 0.5)+
 theme_ridges()+
  labs(x = 'CD45')

median(CD45$AdjVal)

Y%>%
  filter(SampleID == 'PPA')%>%
  filter(CD45 < 0.5)

#map cells using raw threshold value. Value is chosen from corresponding adjval 
Y%>%
  filter(CD45 < 2) %>%
  filter(SampleID == 'PPA')%>%
ggscatter("X_centroid", "Y_centroid", orientation = "reverse",
          size = 0.25, color = 'CD45')

# Review model data after fit to determine thresholds

# lambda - relative fraction of cells assigned to each Gaussian
# mu     - mean of each Gaussian (i.e., location of the peak) in normalized expression space
# sigma  - st.d. of each Gaussian in normalized expression space
# lo/hi  - the interval in the original expression space that the data was normalized to

GMM %>% filter(Marker == 'CD45') %>% pluck('GMM', 1) 

# Use non-normalized dataset to subset CD45 < median

statesCD45neg <- filter(X, CD45 < 2)
as_tibble(statesCD45neg)

states <- statesCD45neg %>% unnest(uniqueID)
states$SampleID <- factor(states$SampleID, levels = c("PP", "PPA", "PPB"), ordered = TRUE)

length(states$DAPI_1)
length(Y$DAPI_1)

# Plot log10 markers for CD45 <2. Related to panel E

ggplot(states, aes(x = pSTAT3, y = fct_rev(as_factor(SampleID)), fill = SampleID))+
  stat_density_ridges(quantile_lines = TRUE, 
                      quantiles = 0.75, 
                      bandwidth = 0.1, 
                      alpha = 0.5)+
  theme_ridges()+
  theme(legend.position = 'none')+
  ylab("Sample")+
  scale_fill_manual(values = c("blue", "red", "green"))

# Thresholded pSTAT3

cells <- states%>%
  group_by(SampleID)%>%
  filter(pSTAT3 > 2) %>%
  tally()

#Count and subset total cells

total <- states%>%
  group_by(SampleID)%>%
  tally()

sum(total$n)

sum_pSTAT3 <- left_join(cells, total, by = 'SampleID')

cohort <- rename(sum_pSTAT3, pSTAT3 = n.x, total = n.y)

# calculate percent difference in markers by treatment
order <- c("PP", "PPA", "PPB")

pct_cohort <- cohort%>%
    mutate(Percentage=paste0(round(pSTAT3/total*100,2)))%>%
  arrange(order)

pct_cohort$Percentage <- as.numeric(pct_cohort$Percentage)


ggbarplot(pct_cohort, x = 'SampleID', y = 'Percentage', fill = 'SampleID',
          xlab = 'sample',
          ylab = 'pSTAT3 cancer cells (%)',
          legend = 'none')+
  theme(axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title = element_text(size = 22),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 24))+
  scale_fill_manual(values = c('blue', 'red', 'green'))

#END