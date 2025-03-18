library(vegan)

d <- read.csv("dist.csv", header = FALSE)
m <- read.csv("meta.csv", header = TRUE)
m$loc <- factor(m$loc)

# adonis2(d ~ batch, m, permutations = 1000000)
# adonis2(d ~ part, m, permutations = 1000000)
# adonis2(d ~ loc, m, permutations = 1000000)
adonis2(d ~ batch + part + loc, m, permutations = 1000000, by = "margin")
