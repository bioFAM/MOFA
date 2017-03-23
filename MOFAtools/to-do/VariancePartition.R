# (1) Stacked barplot

# Collect dimensionalities
N <- object@Dimensions[["N"]]
M <- object@Dimensions[["M"]]
K <- object@Dimensions[["K"]]

# Collect relevant expectations
W <- getExpectations(object,"SW","EW")
S <- getExpectations(object,"SW","ES")
SW <- getExpectations(object,"SW","ESW")
Z <- getExpectations(object,"Z","E")
Y <- getExpectations(object,"Y","E")


# (2) For a particular gene plot how much variation is explained by a given factor and how much is noise

view <- "expr_mRNA"
Y <- object@TrainData[[view]]
W <- getExpectations(object,"SW","E")[[view]]; rownames(W) <- colnames(Y)
Z <- getExpectations(object,"Z","E"); rownames(Z) <- rownames(Y)

genes <- c("ENSG00000003147","ENSG00000005189","ENSG00000005812")

# ypred <- t(as.matrix(W[gene,])) %*% t(Z)

fvar <- sapply(1:K, function(k) apply(as.matrix(Z[,k]) %*% t(as.matrix(W[genes,k])), 2, var) / apply(Y[,genes],2,var))

tmp <- as.data.frame(fvar)
tmp$gene <- rownames(tmp)
tmp <- tidyr::gather(tmp, key="K", value="fvar", -gene)

p <- ggplot2::ggplot(tmp, ggplot2::aes(x=gene,y=fvar,fill=K)) +
  geom_bar(stat="identity", position="stack") +
  ylab("Variance explained") + xlab("") +
  # ylim(c(0,0.5)) + 
  coord_flip() +
  theme(
    # plot.margin=margin(10,10,10,10),
    axis.text.x=element_text(size=rel(1.4), color='black', margin=margin(7,0,0,0)),
    # axis.text.y=element_text(size=rel(1.3), color='black'),
    axis.text.y=element_text(size=rel(1.3), color='black'),
    axis.title.x=element_text(size=rel(1.3), margin=margin(10,0,0,0)),
    axis.title.y=element_blank(),
    axis.line = element_line(colour="black", size=0.8),
    axis.ticks.x = element_line(colour="black", size=0.5),
    axis.ticks.y = element_blank(),
    legend.position='right',
    legend.title=element_text(size=rel(1.1)),
    legend.text=element_text(size=rel(1.0)),
    legend.key=element_rect(fill='transparent'),
    panel.border=element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()
  )
print(p)



