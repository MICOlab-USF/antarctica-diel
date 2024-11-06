# library(devtools)
# install_github("zdk123/SpiecEasi")
library(SpiecEasi)

#-------------------------------------------------------------------------------
# Lets simulate some multivariate data under zero-inflated negative binomial model,
# based on (high depth/count) round 1 of the American gut project, with a sparse network. The basic steps are
# 
# 1. load the data and normalize counts to to common scale (min depth)
# 2. fit count margins to the model
# 3. generate a synthetic network
# 4. generate some synthetic data
# 5. clr transformation
# 6. inverse covariance estimation along a lambda (sparsity) path
# 7. stability selection using the StARS criterion
# 8. evaluate performance
# 
# Obviously, for real data, skip 1-4.

load("Cellvibrio_test.RData",verbose = TRUE)
rownames(df.count) <- df.count$kegg.vec


Cellvibrio <- t(df.count[sample(nrow(df.count),24),-ncol(df.count)])
# Cellvibrio <- as.matrix(df.count[,-ncol(df.count)])
rm(df.count)

depths <- rowSums(Cellvibrio)
Cellvibrio.n  <- t(apply(Cellvibrio, 1, norm_to_total))
Cellvibrio.cs <- round(Cellvibrio.n * min(depths))

d <- ncol(Cellvibrio.cs)
n <- nrow(Cellvibrio.cs)
e <- d

#-------------------------------------------------------------------------------
# Synthesize the data

set.seed(10010)
graph <- make_graph('cluster', d, e)
Prec  <- graph2prec(graph)
Cor   <- cov2cor(prec2cov(Prec))

X <- synth_comm_from_counts(Cellvibrio.cs, mar=2, distr='zinegbin', Sigma=Cor, n=n)

#-------------------------------------------------------------------------------
# the main SPIEC-EASI pipeline: Data transformation, sparse inverse covariance estimation and model selection

se <- spiec.easi(X, method='mb', lambda.min.ratio=1e-2, nlambda=15)
# Applying data transformations...
# Selecting model with pulsar using stars...
# Fitting final estimate with mb...
# done

#-------------------------------------------------------------------------------
# examine ROC over lambda path and PR over the stars index for the selected graph

huge::huge.roc(se$est$path, graph, verbose=FALSE)
stars.pr(getOptMerge(se), graph, verbose=FALSE)
# stars selected final network under: getRefit(se)

#-------------------------------------------------------------------------------
# The above example does not cover all possible options and parameters. For example, 
# other generative network models are available, the lambda.min.ratio 
# (the scaling factor that determines the minimum sparsity/lambda parameter) 
# shown here might not be right for your dataset, and its possible that you’ll 
# want more repetitions (number of subsamples) for StARS.

#-------------------------------------------------------------------------------
# Analysis of American Gut data
# Now let’s apply SpiecEasi directly to the American Gut data. Don’t forget that
# the normalization is performed internally in the spiec.easi function. Also,
# we should use a larger number of stars repetitions for real data. We can pass
# in arguments to the inner stars selection function as a list via the parameter
# pulsar.params. If you have more than one processor available, you can also
# supply a number to ncores. Also, let’s compare results from the MB and glasso
# methods as well as SparCC (correlation).

se.mb.amgut <- spiec.easi(Cellvibrio, method='mb', lambda.min.ratio=1e-2,
                          nlambda=20, pulsar.params=list(rep.num=50))
se.gl.amgut <- spiec.easi(Cellvibrio, method='glasso', lambda.min.ratio=1e-2,
                          nlambda=20, pulsar.params=list(rep.num=50))
sparcc.amgut <- sparcc(Cellvibrio)
## Define arbitrary threshold for SparCC correlation matrix for the graph
sparcc.graph <- abs(sparcc.amgut$Cor) >= 0.3
diag(sparcc.graph) <- 0
library(Matrix)
sparcc.graph <- Matrix(sparcc.graph, sparse=TRUE)
## Create igraph objects
ig.mb     <- adj2igraph(getRefit(se.mb.amgut))
ig.gl     <- adj2igraph(getRefit(se.gl.amgut))
ig.sparcc <- adj2igraph(sparcc.graph)

#-------------------------------------------------------------------------------
# Visualize using igraph plotting:

library(igraph)
## set size of vertex proportional to clr-mean
vsize    <- rowMeans(clr(Cellvibrio, 1))+6
am.coord <- layout_with_fr(ig.mb)

par(mfrow=c(1,3))
# plot(ig.mb, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="MB")

plot(ig.mb, vertex.size=vsize, vertex.label.dist=1.1,
     vertex.label.cex=0.5,
     vertex.label.color="black", main="R1+R2 glasso", layout=am.coord)

plot(ig.gl, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="glasso")
plot(ig.sparcc, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="sparcc")

#-------------------------------------------------------------------------------
# We can evaluate the weights on edges networks using the terms from the
# underlying model. SparCC correlations can be used directly, while SpiecEasi
# networks need to be massaged a bit. Note that since SPIEC-EASI is based on
# penalized estimators, the edge weights are not directly comparable to SparCC
# (or Pearson/Spearman correlation coefficients)

library(Matrix)
secor  <- cov2cor(getOptCov(se.gl.amgut))
sebeta <- symBeta(getOptBeta(se.mb.amgut), mode='maxabs')
elist.gl     <- summary(triu(secor*getRefit(se.gl.amgut), k=1))
elist.mb     <- summary(sebeta)
elist.sparcc <- summary(sparcc.graph*sparcc.amgut$Cor)

hist(elist.sparcc[,3], main='', xlab='edge weights')
hist(elist.mb[,3], add=TRUE, col='forestgreen')
hist(elist.gl[,3], add=TRUE, col='red')

#-------------------------------------------------------------------------------
# Lets look at the degree statistics from the networks inferred by each method.S

dd.gl     <- degree_distribution(ig.gl)
dd.mb     <- degree_distribution(ig.mb)
dd.sparcc <- degree_distribution(ig.sparcc)

plot(0:(length(dd.sparcc)-1), dd.sparcc, ylim=c(0,1), type='b',
     ylab="Frequency", xlab="Degree", main="Degree Distributions")
points(0:(length(dd.gl)-1), dd.gl, col="red" , type='b')
points(0:(length(dd.mb)-1), dd.mb, col="forestgreen", type='b')
legend("topright", c("MB", "glasso", "sparcc"),
       col=c("forestgreen", "red", "black"), pch=1, lty=1)

#-------------------------------------------------------------------------------
# Working with phyloseq
# SpiecEasi includes some convience wrappers to work directly with phyloseq objects.

library(phyloseq)
## Load round 2 of American gut project
data('amgut2.filt.phy')
se.mb.amgut2 <- spiec.easi(amgut2.filt.phy, method='mb', lambda.min.ratio=1e-2,
                           nlambda=20, pulsar.params=list(rep.num=50))
ig2.mb <- adj2igraph(getRefit(se.mb.amgut2),  vertex.attr=list(name=taxa_names(amgut2.filt.phy)))
plot_network(ig2.mb, amgut2.filt.phy, type='taxa', color="Rank3")
