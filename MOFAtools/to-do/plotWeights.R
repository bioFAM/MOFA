
# Plot top weights for a given view and latent variable

m = 1
k = 1
n <- 20

SW <- getExpectations(model,"SW","E")[[m]][,k]
alpha <- getExpectations(model,"Alpha","E")[k,]

df <- data.frame(name=colnames(model@TrainData[[m]]), weight=SW)
df <- df[with(df, order(-weight)), ]
df <- head(df,n)
df$name <- factor(df$name, levels=df$name)

p <- ggplot2::ggplot(df,aes(x=name, y=weight, color=weight)) +
  geom_point(size=5) +
  geom_segment(aes(xend=name, yend=0)) +
  scale_colour_gradient(low="grey", high="black") +
  coord_flip() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=rel(1.3), hjust=1, color='black', margin=margin(0,-8,0,0)),
        axis.text.x = element_text(size=rel(1.5), color='black'),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_line(),
        legend.position='none',
        legend.text=element_text(size=rel(1.5),angle=90, hjust=1, color="black"),
        legend.key=element_rect(fill='transparent'),
        panel.background = element_blank(),
        aspect.ratio = .7)
print(p)

