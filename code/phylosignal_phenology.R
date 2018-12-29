rm(list = ls())
# First read in the arguments listed at the command line
args=(commandArgs(TRUE))
library(picante)
library(ape)
library(phylobase)
library(adephylo)
library(phytools)

m <- read.csv("data/species_model_coefficients.csv", stringsAsFactors = FALSE)

m$species <- gsub(" ", "_", as.character(m$species))
row.names(m) <- m[,1]
m1 <- m[,8, drop=FALSE]


t1 <- read.nexus("data/Protea_phylogeny.tre")

subphy <- match.phylo.data(t1, m1)$phy
subdat <- match.phylo.data(t1, m1)$data

# 1. Moran's I
# Combine both the phylogenetic and trait data into one file
phylotraits <- phylo4d(subphy, subdat)

# Run a Moran test using some Monte Carlo simulations (default is 999 randomizations)
moran.test <- abouheif.moran(phylotraits, method="Abouheif", nrepet = 999)
moran.test

# 2. Abouheif's Cmean
# Combine both the phylogenetic and trait data into one file
phylotraits <- phylo4d(subphy, subdat)
abouheif.test <- abouheif.moran(phylotraits, method="oriAbouheif")
abouheif.test

# 3. Pagel's λ
# First, you need to define which trait you want to test and give names to each value according to species
trait <- subdat[,1]
names(trait) <- rownames(subdat)
phylosig(subphy, trait, method="lambda", test=TRUE, nsim=1000)

# 4. Blomberg's K
# First, you need to define which trait you want to test and give names to each value according to species
trait <- subdat[,1]
names(trait) <- rownames(subdat)

# We do the test with 999 randomizations:
phylosig(subphy, trait, method="K", test=TRUE, nsim=1000)



