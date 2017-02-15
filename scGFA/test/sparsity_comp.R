library(ggplot2)
library(reshape2)

# sparsities 
inferred_sparsity_no_opt = read.table('/tmp/test/model/sparsity_inf_no_opt.txt', header=TRUE)
inferred_sparsity_opt = read.table('/tmp/test/model/sparsity_inf_opt.txt', header=TRUE)
true_sparsity = read.table('/tmp/test/model/sparsity_truth.txt', header=TRUE)

# alpha
alpha_no_opt = read.table('/tmp/test/model/alpha_no_opt.txt', header=TRUE)
alpha_opt = read.table('/tmp/test/model/alpha_opt.txt', header=TRUE)

alpha_no_opt_N =  apply(alpha_no_opt[,-c(1,2)], 1, function(x) sum(x != 0))
alpha_no_opt_N= cbind(alpha_no_opt[,c(1,2)], alpha_no_opt_N)
alpha_no_opt_mean =  apply(alpha_no_opt[,-c(1,2)], 1, sum)/alpha_no_opt_N[,3]
alpha_no_opt_mean= cbind(alpha_no_opt[,c(1,2)], alpha_no_opt_mean)

alpha_opt_N =  apply(alpha_opt[,-c(1,2)], 1, function(x) sum(x != 0))
alpha_opt_N= cbind(alpha_opt[,c(1,2)], alpha_opt_N)
alpha_opt_mean =  apply(alpha_opt[,-c(1,2)], 1, sum)/alpha_opt_N[,3]
alpha_opt_mean= cbind(alpha_opt[,c(1,2)], alpha_opt_mean)

sparsity_opt = merge(inferred_sparsity_opt, alpha_opt_N, by=c('run', 'factor'))
sparsity_opt = merge(sparsity_opt, alpha_opt_mean, by=c('run', 'factor'))
sparsity_opt = merge(true_sparsity, sparsity_opt, by=c('run', 'factor'))
sparsity_opt =  cbind(sparsity_opt, 'optimisation')
colnames(sparsity_opt) = cbind('run', 'factor', 'truth', 'inferred_sparsity', 'alive_factors', 'alpha_average','method')

sparsity_no_opt = merge(inferred_sparsity_no_opt, alpha_no_opt_N, by=c('run', 'factor'))
sparsity_no_opt = merge(sparsity_no_opt, alpha_no_opt_mean, by=c('run', 'factor'))
sparsity_no_opt = merge(true_sparsity, sparsity_no_opt, by=c('run', 'factor'))
sparsity_no_opt = cbind(sparsity_no_opt, 'no optimisation')
colnames(sparsity_no_opt) = cbind('run', 'factor', 'truth', 'inferred_sparsity', 'alive_factors', 'alpha_average', 'method')

sparsities = rbind(sparsity_opt, sparsity_no_opt)


ggplot(sparsities, aes(x=truth, y = inferred_sparsity))+
  geom_point(aes(shape=method, color=alive_factors))+
  scale_colour_gradient2(low = 'red', high='blue', mid = 'green', midpoint = 6)

ggplot(sparsities, aes(x=truth, y = inferred_sparsity))+
  geom_point(aes(shape=method, color=alpha_average))+
  scale_colour_gradient2(low = 'red', high='blue', mid='green')

ggplot(sparsities, aes(x=alive_factors))+
  geom_histogram(aes(fill=method), position = 'dodge')
