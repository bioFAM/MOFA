# Scatterplot of the latent variables

df <- data.frame(
  sample=rownames(Z),
  met=Z[,7],
  pluri=Z[,2]
)

df$Phenotype <- c("a","b")[as.numeric(df$pluri>0.5 & df$met < 0)+1]
df[df$sample %in% c("D05","E01"),"Phenotype"] <- "c"

ggplot(df,aes(x=met,y=pluri,color=Phenotype)) + 
  geom_point(size=2.5) +
  geom_text(aes(label=sample),hjust=0, vjust=0) +
  xlab("Latent variable 7\n(Global methylation axis)") +
  ylab("Latent variable 2\n(Pluripotency axis)") +
  theme(
    axis.text=element_text(size=rel(1.2), color='black'),
    axis.title=element_text(size=rel(1.5)),
    axis.title.y=element_text(margin=margin(0,15,0,0)),
    axis.title.x=element_text(margin=margin(15,,0,0)),
    axis.line = element_line(colour="black", size=1),
    axis.ticks = element_line(colour="black", size=0.5),
    legend.position='none',
    panel.border=element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()
  )
