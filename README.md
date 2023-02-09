missing_data_pca_tutorial
================
Marc A. Beer
2/8/2023

<br><br><br>

# Nonrandom missing data in population genetics analysis (and how to correct it)

<br> <font size="5"> The randomness (or lack thereof) of missing data is
an important consideration for analysis. Population genetic datasets
frequently contain missing genotypes due to variation in sequencing
quality across genetic loci and individual samples. When these missing
data are random with respect to key characteristics of the data (e.g.,
with respect to samples’ geographic regions of origin), we expect that
they will lead to little distortion to the patterns we study; we will
see a demonstration of this in just a moment. However, when data are
missing not at random (MNAR), genetic relationships among samples can
indeed become distorted and adversely affect inference. A more thorough
treatment of random and non-random missing data can be found at a
wonderful series of webpages by Stef van Buuren
(<https://stefvanbuuren.name/fimd/sec-MCAR.html>) </font> <br>

<font size="5"> You might be wondering why genetic data MNAR might
occur. What might lead to a correlation between missing data and some
characteristics of our samples? One example is the improper
randomization of samples into multiplexed sequencing libraries. If we
sample individuals from six populations and treat each population’s
samples separately (i.e., when preparing sequencing libraries and/or
carrying out sequencing itself), nonrandomness in missing data or other
nonrandom aberrations may result. In this case, some sequencing runs may
go poorly, leading to poor sequencing depth and missing data that is
correlated with samples’ populations of origin; this misfortune could
have been avoided by randomizing samples into sequencing libraries with
respect to population of origin. Another case that is salient to
reduced-representation sequencing is inter-population restriction site
polymorphism. Some populations may lack restriction sites that are
present in other populations. The former populations would have
completely missing data at loci that are otherwise genotyped in other
populations. </font> <br>

<font size="5"> Okay, the potential for genetic data MNAR are clearly
present, but besides some nebulous warning that it can cause distrotion,
should we really care? In this walkthrough, I’ll show you how data MNAR
can distort inference from a common population genetic analysis:
principal components analysis (PCA). Then, I’ll show you how to remedy
the problem using a related analysis, discriminant analysis of principal
components (DAPC). </font> <br>

<font size="5"> Let’s start by loading in some packages. One of these is
adegenet, which contains numerous functions for carrying out population
genetics analyses, including PCA and DAPC. Handily, it also includes
genetic datasets that we can use as examples. </font>

``` r
library(adegenet)
library(dplyr)
library(knitr)
library(ggplot2)
library(patchwork)
library(rmarkdown)
```

<br>

## 1. Load in a population genetic dataset

<br> <font size="5"> The package adegenet contains several genetic
datasets, and we will use one as an example. Let’s load in a dataset
containing genotypes for 30 microsatellite loci for 600 diploid
individuals from six populations simulated under an island model. Much
of this walkthrough will involve the locus/allele matrix stored in the
tab subset of the dataset. If you are less familiar with genind objects,
take a few minutes to look through some of the data object’s
subcomponents. Let’s visualize the locus/allele matrix below. Rows are
individual samples and each column refers to a locus and one of its
alleles. E.g., column names loc-1.03 and loc-1.19 refer to two different
alleles at locus 1; entries in the matrix indicate how many copies of
given allele are present (ranging from 0 to 2). Your own datasets may
look different: while some of the microsatellite loci in this example
dataset have many alleles, the SNPs typical of modern genomic datasets
often have only two (and thus biallelic loci uniformly will only have
two columns per locus). </font> <br>

``` r
#####
#load in example data distributed with adegenet
##in genind format, which is useful for our purposes
data(dapcIllus)

#we'll use only the first dataset
data <- dapcIllus$a

#isolating the genotype matrix
data_tab <- data$tab

#rows are individuals while columns record the presence/absence of a given allele for a codominant locus
#notice that column names refer to a locus (e.g., loc-1 or loc-2) followed by the allele of interest
#(e.g., loc-1.03 and loc-1.19 refer to two different alleles at locus 1)
#for biallelic loci, such as SNPs, there will uniformly be only two columns per locus
kable(data_tab[1:5,1:10])
```

| loc-1.03 | loc-1.19 | loc-2.01 | loc-2.04 | loc-2.40 | loc-2.41 | loc-2.44 | loc-3.08 | loc-3.10 | loc-3.38 |
|---------:|---------:|---------:|---------:|---------:|---------:|---------:|---------:|---------:|---------:|
|        1 |        1 |        0 |        0 |        1 |        0 |        1 |        0 |        2 |        0 |
|        1 |        1 |        0 |        0 |        2 |        0 |        0 |        0 |        2 |        0 |
|        1 |        1 |        0 |        0 |        1 |        0 |        1 |        0 |        2 |        0 |
|        2 |        0 |        0 |        0 |        0 |        0 |        2 |        0 |        0 |        0 |
|        1 |        1 |        0 |        0 |        1 |        0 |        1 |        0 |        1 |        0 |

<br>

<font size="5"> Note that this is a complete dataset - there are no
missing data. We will introduce missing data both randomly and
nonrandomly to demonstrate their effects on results from principal
components analysis (PCA). To get a baseline understanding of the
complete dataset, let’s run PCA on it now. </font>

``` r
#run PCA on the complete dataset
data_pca <- dudi.pca(data, center=TRUE, scale=TRUE, scannf=FALSE, nf=2)

#summarize PC coordinates by population
pc_pop <- data.frame(pop=data$pop, data_pca$li)
pc_pop_summ <- pc_pop %>%
  group_by(pop) %>%
  summarise(pc1_mean=mean(Axis1), pc2_mean=mean(Axis2))

#plot results
p_fulldata <- ggplot()+
  geom_vline(xintercept=0)+
  geom_hline(yintercept=0)+
  geom_point(data_pca$li, mapping=aes(x=Axis1, y=Axis2, colour=data$pop), size=3, alpha=0.5)+
  geom_point(pc_pop_summ, mapping=aes(x=pc1_mean, y=pc2_mean, fill=pop), colour="black", shape=24, size=4)+
    labs(x="PC1", y="PC2", colour="Population", fill="Population \n mean", title = "Full dataset (0% missing data)")+
  scale_y_continuous(limits=c(-10, 10), breaks=seq(from=-10, to=10, by=2))+
  scale_x_continuous(limits=c(-10, 10), breaks=seq(from=-10, to=10, by=2))+
  theme(aspect.ratio=1)+
  coord_fixed()+
  theme_bw()
```

<br>

<font size="5"> The results of PCA are visualized below. Individuals
from different populations (coloured points) clearly cluster together.
Noteworthy is the relatively large distance between populations P2 and
P4 along PC2 (the vertical axis). The relationship between these two
populations will be used to understand how random and nonrandom missing
data can impact our inference </font>
<img src="missing_data_pca_tutorial_files/figure-gfm/unnamed-chunk-4-1.png" style="display: block; margin: auto;" />
<br>

## 2. Data missing at random

<br>

<font size="5"> Before investigating nonrandom missing data, let’s start
with the best-case scenario of random missing data. </font> <br>

<font size="5"> We will convert the genotypes of 20% (6) of the loci to
missing data for 50% (150) of the individuals, leading to a missing data
rate of 10% in terms of a genotype matrix; note that this percentage
should be roughly the same in terms of the tab matrix in the genind
object, with some variation since loci can have different numbers of
alleles. The code is shown and annotated below, but briefly, we first
randomly select 20% of the loci, and then for each of those loci, we
randomly select 50% of the individuals to have missing data. We will
repeat this process selecting 80% of the loci and 50% of the
individuals, leading to a missing data rate of 40%. Thus, we have
randomly introduced two different magnitudes of missing data in a way
that is uncorrelated with properties of our samples (e.g., their
populations of origin). If we wanted, we could introduce missing data
more uniformly across loci, but the above procedure is sufficient for
our purposes. </font>

``` r
###
#simulate random missing data (10%)

#save a copy of the data_tab to modify
data_random <- data
data_tab_random <- data_tab

#we will select 20% of the loci to add missing data to
set.seed(100)
loc_select <- sample(x=levels(data$loc.fac),
              size=round(0.20*length(levels(data$loc.fac)), digits=0),
              replace=FALSE
              )

#for each locus, we will randomly select 50% of individuals (i.e., 150 individuals) to introduce missing data to
for (i in 1:length(loc_select)){
  
  #randomly sample from all the individuals in the dataset
  set.seed(100+i)
  ind <- sort(sample(x=rownames(data$tab),
                size=round(0.50*length(rownames(data$tab)), digits=0),
                replace=FALSE
          ),
          decreasing=FALSE)
  
  #for the current locus, replace genotypes of randomly selected individuals with NA values
  data_tab_random[ind, which(data$loc.fac==loc_select[i])] <- NA
}

###
#simulate random missing data (40%)

#save a copy of the data_tab to modify
data_random_040 <- data
data_tab_random_040 <- data_tab

#we will select 80% of the loci to add missing data to
set.seed(100)
loc_select <- sample(x=levels(data$loc.fac),
                     size=round(0.8*length(levels(data$loc.fac)), digits=0),
                     replace=FALSE
)

#for each locus, we will randomly select 50% of individuals (i.e., 300 individuals) to introduce missing data to
for (i in 1:length(loc_select)){
  
  #randomly sample from all the individuals in the dataset
  set.seed(100+i)
  ind <- sort(sample(x=rownames(data$tab),
                     size=round(0.5*length(rownames(data$tab)), digits=0),
                     replace=FALSE
  ),
  decreasing=FALSE)
  
  #for the current locus, replace genotypes of randomly selected individuals with NA values
  data_tab_random_040[ind, which(data$loc.fac==loc_select[i])] <- NA
}
```

<br>

<font size="5"> Let’s conduct PCA on the modified datasets. Note that
PCA requires complete data, so we will first fill in the now-missing
data using mean value imputation. The code below both imputes missing
data and carries out PCA. </font>

``` r
#re-do PCA to see effects

###
#10% missing data

#insert tab with missing data back into the genind object
data_random$tab <- data_tab_random

#note that there cannot be missing entries in PCA, so we will impute the now-missing data
data_random_imputed <- scaleGen(x=data_random, center=TRUE, scale=TRUE, NA.method="mean")

#run the PCA
data_random_pca <- dudi.pca(data_random_imputed, center=TRUE, scale=TRUE, scannf=FALSE, nf=2)

#summarize PC coordinates by population
pc_random_pop <- data.frame(pop=data$pop, data_random_pca$li)
pc_random_pop_summ <- pc_random_pop %>%
  group_by(pop) %>%
  summarise(pc1_mean=mean(Axis1), pc2_mean=mean(Axis2))

p_random <- ggplot()+
  geom_vline(xintercept=0)+
  geom_hline(yintercept=0)+
  geom_point(data_random_pca$li, mapping=aes(x=Axis1, y=Axis2, colour=data$pop), size=3, alpha=0.5)+
  geom_point(pc_random_pop_summ, mapping=aes(x=pc1_mean, y=pc2_mean, fill=pop), colour="black", shape=24, size=4)+
  labs(x="PC1", y="PC2", colour="Population", fill="Population \n mean", title = "10% random missing data")+
  scale_y_continuous(limits=c(-10, 10), breaks=seq(from=-10, to=10, by=2))+
  scale_x_continuous(limits=c(-10, 10), breaks=seq(from=-10, to=10, by=2))+
  theme(aspect.ratio=1)+
  coord_fixed()+
  theme_bw()

###
#40% missing data

#insert tab with missing data back into the genind object
data_random_040$tab <- data_tab_random_040

#note that there cannot be missing entries in PCA, so we will impute the now-missing data
data_random_040_imputed <- scaleGen(x=data_random_040, center=TRUE, scale=TRUE, NA.method="mean")

#run the PCA
data_random_040_pca <- dudi.pca(data_random_040_imputed, center=TRUE, scale=TRUE, scannf=FALSE, nf=2)

#summarize PC coordinates by population
pc_random_040_pop <- data.frame(pop=data$pop, data_random_040_pca$li)
pc_random_040_pop_summ <- pc_random_040_pop %>%
  group_by(pop) %>%
  summarise(pc1_mean=mean(Axis1), pc2_mean=mean(Axis2))

p_random_040 <- ggplot()+
  geom_vline(xintercept=0)+
  geom_hline(yintercept=0)+
  geom_point(data_random_040_pca$li, mapping=aes(x=Axis1, y=Axis2, colour=data$pop), size=3, alpha=0.5)+
  geom_point(pc_random_040_pop_summ, mapping=aes(x=pc1_mean, y=pc2_mean, fill=pop), colour="black", shape=24, size=4)+
    labs(x="PC1", y="PC2", colour="Population", fill="Population \n mean", title = "40% random missing data")+
  scale_y_continuous(limits=c(-10, 10), breaks=seq(from=-10, to=10, by=2))+
  scale_x_continuous(limits=c(-10, 10), breaks=seq(from=-10, to=10, by=2))+
  theme(aspect.ratio=1)+
  coord_fixed()+
  theme_bw()
```

<br>

<font size="5"> The results of PCA on the two imputed datasets are shown
below (panels B and C), adjacent to the original PCA we conducted on the
complete, unimputed dataset (Panel A). At a respectable missing data
rate of 10%, the broad patterns are largely undistorted. Even the
missing data rate of 40% captures most of the patterns found in the
complete dataset. </font> <br>

<font size="5"> However, you might notice that individual PC scores, as
well as population means, are shifted towards the graph origin;
distances between population means have accordingly decreased. This
“collapse” towards the origin is a symptom of mean value imputation, and
the obscuration of the true patterns becomes more extreme with
increasing missing data. More sophisticated imputation methods (such as
one based on snmf in the R package LEA) may combat this distortion.
Indeed, it would be fair to say that the observed distortion is a
consequence of both missing data and imputation. In general, however,
the effects of random missing data are rather mild because individuals
of all populations are affected similarly. </font>

<img src="missing_data_pca_tutorial_files/figure-gfm/unnamed-chunk-7-1.png" style="display: block; margin: auto;" />
<br>
