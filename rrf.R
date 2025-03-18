library(vegan)

otu_table <- read.csv("rarefaction_table.csv", row.names = 1)
nrow(otu_table)
ncol(otu_table)

# otu_table <- head(otu_table, 3)
png("Rplots.png", width = 2400, height = 1800, res = 300)

rarecurve(otu_table,
    step = 1000, cex = 0.6, label = FALSE, col = "blue",
    main = "Rarefaction Curves",
    xlab = "Sequencing Depth (Reads)", ylab = "Number of OTUs"
)

# Add a legend
# legend("bottomright",
#     legend = rownames(otu_table), col = "blue",
#     lty = 1, cex = 0.5
# )
