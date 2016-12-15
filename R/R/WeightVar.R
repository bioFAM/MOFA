
library(dplyr)
library(tidyr)
library(ggplot2)

## Plot weights (upper) and variance (lower) ##

# Plot average weight and average variance

m = 1
k = 1
SW <- abs(model@Expectations$SW[[m]]$ESW[,k])
# S <- model@Expectations$SW[[m]]$ES[,k]
# W <- model@Expectations$SW[[m]]$EW[,k]
model@Expectations$Alpha[[m]]$E

n <- 20
df <- data.frame(name=colnames(model@Data[[m]]), weight=SW) %>% arrange(desc(weight)) %>% head(n) %>%
  mutate(name=factor(name, levels=name))

p <- ggplot(df,aes(x=name, y=weight, color=weight)) +
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


