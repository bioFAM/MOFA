
## Correlation plot ##
if (opts$corplot) {
  suppressMessages(library(corrplot))
  h <- cor(Z,Z)
  # pdf(file=sprintf("%s/cor.pdf",outdir))
  corrplot(h, method='color', order="original", diag=T, type="upper",
           tl.col="black", tl.cex=1.8, cl.cex=0.9)
  # dev.off()
}



## Faced scatterplot (correlation between latent variables) ##
if (opts$scatterplot) {
  
  # Center and scale latent variable matrix
  Z <- model$Z
  Z <- scale(Z,center=T,scale=T)
  
  # Convert latent variable matrix into dataframe
  z_order <- 1:nrow(alpha)
  Z <- Z[,as.numeric(z_order)]; colnames(Z) <- z_order
  Z_long <- as.data.frame(t(Z)) %>% tbl_df %>% mutate(Z=factor(1:n())) %>%
    gather(sample, value, -Z)
  
  # Concate PCs
  joined <- Z_long %>% inner_join(Z_long, by='sample') %>%
    dplyr::rename(Z1=Z.x, Z2=Z.y, value1=value.x, value2=value.y)
  
  # Keep the following components
  # keep <- 1:10
  # d <- joined %>% filter(Z1 %in% keep, Z2 %in% keep) %>% droplevels
  
  # Plot
  #png(filename=file.path(plots_outdir,"scatterplot.png"), height=1800, width=1800)
  postscript(file=sprintf("%s/scatterplot.eps",outdir))
  p <- ggplot(joined, aes(x=value1, y=value2)) +
    ggtitle("") +
    stat_smooth(method=lm, color='black') +
    geom_point(aes(color=sample), size=0.5, alpha=0.5) +
    xlab('') + ylab('') +
    facet_grid(Z1~Z2) +
    guides(color=F)
  print(p)
  dev.off()
}
