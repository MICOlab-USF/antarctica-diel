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

data(amgut1.filt)
depths <- rowSums(amgut1.filt)
amgut1.filt.n  <- t(apply(amgut1.filt, 1, norm_to_total))
amgut1.filt.cs <- round(amgut1.filt.n * min(depths))

d <- ncol(amgut1.filt.cs)
n <- nrow(amgut1.filt.cs)
e <- d

#-------------------------------------------------------------------------------
# Synthesize the data

set.seed(10010)
graph <- make_graph('cluster', d, e)
Prec  <- graph2prec(graph)
Cor   <- cov2cor(prec2cov(Prec))

X <- synth_comm_from_counts(amgut1.filt.cs, mar=2, distr='zinegbin', Sigma=Cor, n=n)

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