library("ape")
library("U.PhyloMaker")
library("maps")
library("phytools")
library("readxl")


#-------------------------------------
# first time phylogeny construction
#-------------------------------------

data <- readxl::read_xlsx(path = "../data/chapter2/FRED/subsets/final.xlsx", sheet = "final")
data$binominal <- gsub(data$binominal, pattern = ' ', replacement = '_') # replace the spaces with _, so the binominal names match the phylogeny
# sp.list needs to have the following columns => species, genus, family, species.relative, genus.relative
species_list <- data.frame(species = data$binominal, genus = data$F01286, family = data$F01289, species.relative = NA, genus.relative = NA)
genus_list <- data.frame(genus = data$F01286, family = data$F01289) # genus and family relationships
genus_list <- genus_list[!duplicated(genus_list), ] # remove duplicated rows

megatree <- ape::read.tree("../data/chapter2/uphylomaker/GBOTB.extended.TPL.tre") # TPL megatree

runtime <- Sys.time()
phylogeny <- U.PhyloMaker::phylo.maker(sp.list = species_list, tree = megatree, gen.list = genus_list)
runtime <- Sys.time() - runtime

# serialize the new phylogenetic tree
ape::write.tree(phy = phylogeny$phylo, file = "../data/chapter2/uphylomaker/FRED4_1301.tre")


#---------------------------------------------------------------
# name match the finalized root trait data with the phylogeny
#---------------------------------------------------------------

tree <- ape::read.tree("../data/chapter2/uphylomaker/FRED4_1301.tre")
stopifnot(all(tree$tip_label %in% data$binominal))

data <- data[match(tree$tip.label, data$binominal), ] # name match the trait data with the phylogeny
stopifnot(all(data$binominal == tree$tip.label))

# serialize the name matched table
write.table(x = data, file = "../data/chapter2/FRED/subsets/name_matched_FRED4_1301.csv", sep = ',', row.names = FALSE)