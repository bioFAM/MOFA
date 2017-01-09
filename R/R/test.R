
setwd("/Users/ricard/git/scGFA/R/R")


source("loadModel.R")
file = "/Users/ricard/git/britta/scGFA/expr/model1.hdf5"
model <- loadModel(file)

source("proportion_variance.R")
# prvar_mk <- CalculateProportionResidualVariance(model, plot=T)

# View vs Factor plot
source("view_vs_factor.R")
PlotViewVsFactor(model, color=colorRampPalette(c("grey100", "grey0"))(100), title="", cluster_cols=F, cluster_rows=F, show_rownames=T, show_colnames=T,
                 legend=T, treeheight_row=20, treeheight_col=20, fontsize_row=20, fontsize_col=20, cellheight=NA, cellwidth=NA, outfile=NA)
  

# Boxplot of weights of gene sets in a given latent varaible and view
abs_value = TRUE
m=1
k=1
x=""
gene_sets=list(c(1,2,3), c(4,5,6), c(7,8,9))
xlabel <- ""
ylabel <- ""


if is.null(names(gene_sets))
  names(gene_sets) <- str_c("GeneSet_",1:length(gene_sets))

tmp <- list()
for (gs in names(gene_sets)) {
  
  if (is.numeric(gene_sets[[gs]])) {
    stopifnot(gene_sets[[gs]] %in% 1:D[m])
  } else if (all(is.character(gene_sets[[gs]]))) {
    stopifnot(all(gene_sets[[gs]] %in% getfeatureNames(model)[[m]]))
  } else {
    stop("Features have to be specified either as character or numeric vector")
  }

  tmp[[gs]] <- model@Expectations$SW[[m]]$E[gene_sets[[gs]],k] %>% as.data.frame %>% 
    `colnames<-`(str_c("K=",k)) %>% `rownames<-`(features) %>% tibble::rownames_to_column("features") %>%
    mutate(geneset=gs) %>% gather(k,value,-features,-geneset)
}
df <- do.call(rbind,tmp)  
df <- data.table::rbindlist(tmp)
  
if (abs_value)
  df$value <- abs(df$value)

ggplot(df, aes(x=geneset, y=value, fill=geneset)) +
  # ggtitle("") +
  geom_boxplot() +
  xlab(xlabel) + ylab(ylabel) +
  theme(
    axis.title.y = element_text(colour="black", size=16, vjust=1.5),
    axis.title.x = element_text(colour="black", size=16, vjust=1.5),
    axis.text.x = element_text(colour="black",size=rel(1.3)),
    axis.text.y = element_text(colour="black",size=rel(1.3)),
    axis.line = element_line(colour="black"),
    axis.ticks.y = element_line(colour="black", size=rel(0.8)),
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position="none"
  )





  
# Plot top weights
abs_value = TRUE
m=1
k=1
xlabel <- ""
ylabel <- ""
n <- 30

weights <-  model@Expectations$SW[[m]]$E[,k]
# names(weights) <- getfeatureNames(model)[[m]]
names(weights) <- str_c("feature",1:length(weights))

if (abs_value==T)
  weights <- abs(weights)

df <- data.frame(value=sort(weights,decreasing=T)[1:n]) %>% tibble::rownames_to_column("feature") %>%
  mutate(feature=factor(feature, levels=feature))

ggplot(df, aes(x=feature, y=value, color=value)) +
  geom_point(size=3) +
  geom_segment(aes(xend=feature, yend=0)) +
  scale_colour_gradient(low="grey", high="black") +
  coord_flip() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=rel(1.2), hjust=1, color='black', margin=margin(0,-8,0,0)),
        axis.text.x = element_text(size=rel(1.5), color='black'),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_line(),
        legend.position='none',
        legend.text=element_text(size=rel(1.5),angle=90, hjust=1, color="black"),
        legend.key=element_rect(fill='transparent'),
        panel.background = element_blank(),
        aspect.ratio = .7)

# ggplot(df, aes(x=feature, y=value)) +
#   ggtitle("") +
#   geom_bar(stat='identity', position="dodge", width=0.5) +
#   scale_colour_gradient(low="grey", high="black") +
#   xlab("") + ylab("") +
#   coord_flip() +
#   guides(fill=F) +
#   theme(
#     axis.title.y = element_text(colour="black", size=16, vjust=1.5),
#     axis.title.x = element_text(colour="black", size=16, vjust=1.5),
#     axis.text.x = element_text(colour="black",size=rel(1.3)),
#     axis.text.y = element_text(colour="black",size=rel(1.3)),
#     axis.line = element_line(colour="black"),
#     axis.ticks.x = element_line(colour="black", size=rel(0.8)),
#     axis.ticks.y = element_blank(),
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_blank()
#   )


# Training statistics
source("training_statistics.R")
elbo_trainCurve(model)
activeK_trainCurve(model)

# Scatterplot
source("scatterplot.R")
scatterPlot(model, idx=2, idy=4)


# Corplot
source("corPlot.R")
FactorsCorrPlot(model)

  
# Y <- model@TrainData
e <- model@Expectations
p <- model@Parameters
data <- model@Data


e$Alpha$bernoulli$E
colMeans(abs(e$SW$bernoulli$EW))
e$Theta$bernoulli$E
abs(e$SW$bernoulli$EW)[,1]
e$SW$bernoulli$EW[e$SW$bernoulli$ES[,1] > 0.75,1]

i <- which(e$SW$bernoulli$ES[,1] > 0.75)
apply(e$Y$bernoulli$E,2,var)[i]
mean(apply(e$Y$bernoulli$E,2,var))


getTrainOpts(model)
getTrainStats(model)

m = "foo"
e$Alpha[[m]]$E
e$Theta[[m]]$E
k=1
e$Alpha[[m]]$E[k]
e$Theta[[m]]$E[k]

p$Alpha[[m]]$a[k]
p$Alpha[[m]]$b[k]

p$Theta[[m]]$a[k]
p$Theta[[m]]$b[k]

mean(e$SW[[m]]$ES[,k])
mean(e$SW[[m]]$ES[,k] > 0.75)

e$SW[[m]]$ES[,k]
plot(e$Z$E[,k])

