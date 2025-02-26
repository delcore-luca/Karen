% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{get.cdn}
\alias{get.cdn}
\title{Get the cell differentiation network from a fitted Kalman Reaction Network.}
\usage{
get.cdn(res.fit, edges.lab = FALSE, AIC = FALSE, cell.cols = NULL)
}
\arguments{
\item{res.fit}{A list returned by get.fit() containing the information of a fitted Kalman Reaction Network.}

\item{edges.lab}{(logical) Defaults to FALSE, in which case the labels (weights) will not be printed on the network edges.}

\item{AIC}{(logical) Defaults to FALSE, in which case the Akaike Information Criterion is not reported.}

\item{cell.cols}{Color legend for the cell types. Defaults to NULL, in which case no color legend for the cell types is provided.}
}
\description{
This function returns the cell differentiation network from a Kalman Reaction Network previously fitted on a clonal tracking dataset.
}
\examples{
\donttest{
cat("\nInstall/load packages")
inst.pkgs <- installed.packages() ## installed packages

l.pkgs <- c("expm",
            "Matrix",
            "parallel",
            "gaussquad",
            "splines",
            "scales",
            "mvtnorm",
            "tmvtnorm",
            "MASS",
            "igraph",
            "stringr",
            "Karen")
lapply(l.pkgs, function(pkg){
  if(!(pkg \%in\% rownames(inst.pkgs))){
    install.packages(pkg)
  }
})
## load packages
lapply(l.pkgs, function(pkg){library(pkg,character.only=TRUE)})

rm(list = ls())

rcts <- c("HSC->P1", ## reactions
          "HSC->P2",
          "P1->T",
          "P1->B",
          "P1->NK",
          "P2->G",
          "P2->M",
          "T->0",
          "B->0",
          "NK->0",
          "G->0",
          "M->0"
          ,"HSC->1"
          ,"P1->1"
          ,"P2->1"
)

cnstr <- c("theta\\\\[\\\\'HSC->P1\\\\'\\\\]=(theta\\\\[\\\\'P1->T\\\\'\\\\] + theta\\\\[\\\\'P1->B\\\\'\\\\] + theta\\\\[\\\\'P1->NK\\\\'\\\\])",
           "theta\\\\[\\\\'HSC->P2\\\\'\\\\]=(theta\\\\[\\\\'P2->G\\\\'\\\\] + theta\\\\[\\\\'P2->M\\\\'\\\\])") ## reaction constraints
latsts <- c("HSC", "P1", "P2") ## latent cell types

ctps <- unique(setdiff(c(sapply(rcts, function(r){ ## all cell types
  as.vector(unlist(strsplit(r, split = "->", fixed = T)))
}, simplify = "array")), c("0", "1")))

########## TRUE PARAMETERS ##########
th.true <- c(0.65, 0.9, 0.925, 0.975, 0.55, 3.5, 3.1, 4, 3.7, 4.1, 0.25, 0.225, 0.275) ## dynamic parameters
names(th.true) <- tail(rcts, -length(cnstr))
s2.true <- 1e-8 ## additonal noise
r0.true <- .1 ## intercept noise parameter
r1.true <- .5 ## slope noise parameter
phi.true <- c(th.true, r0.true, r1.true) ## whole vector parameter
names(phi.true) <- c(names(th.true), "r0", "r1")

########## SIMULATION PARAMETERS ##########
S <- 1000 ## trajectories length
nCL <- 3 ## number of clones
X0 <- rep(0, length(ctps)) ## initial condition
names(X0) <- ctps
X0["HSC"] <- 100
ntps <- 30 ## number of time-points
f_NA <- .75 ## fraction of observed data

###########################
## SIMULATE TRAJECTORIES ##
###########################

XY <- get.sim.trajectories(rct.lst = rcts,
                           constr.lst = cnstr,
                           latSts.lst = latsts,
                           ct.lst = ctps,
                           th = th.true,
                           S = S,
                           nCL = nCL,
                           X0 = X0,
                           s2 = s2.true,
                           r0 = r0.true,
                           r1 = r1.true,
                           f = f_NA,
                           ntps = ntps,
                           trunc = FALSE)

#####################################
## Fitting Karen on simulated data ##
#####################################

nProc <- 1 # number of cores
cat(paste("\tLoading CPU cluster...\n", sep = ""))
cat(paste("Cluster type: ", "PSOCK\n", sep = ""))
cpu <- Sys.getenv("SLURM_CPUS_ON_NODE", nProc) ## define cluster CPUs
hosts <- rep("localhost",cpu)
cl <- makeCluster(hosts, type = "PSOCK") ## make the cluster
rm(nProc)

## mean vector of the initial condition:
m_0 <- replicate(nCL, X0, simplify = "array")
colnames(m_0) <- 1:nCL
## covariance matrix of the initial condition:
P_0 <- Diagonal(length(ctps) * nCL, 1e-5)
rownames(P_0) <- colnames(P_0) <- rep(1:nCL, each = length(ctps))
## Fit Karen on the simulated data:
res.fit <- get.fit(rct.lst = rcts,
                   constr.lst = cnstr,
                   latSts.lst = latsts,
                   ct.lst = ctps,
                   Y = XY$Y[,setdiff(ctps, latsts),],
                   m0 = m_0,
                   P0 = P_0,
                   cl = cl,
                   list(nLQR = 3,
                        lmm = 25,
                        pgtol = 0,
                        relErrfct = 1e-9,
                        tol = 1e-9,
                        maxit = 1000,
                        maxitEM = 10,
                        trace = 1,
                        FORCEP = FALSE))

stopCluster(cl) ## stop the cluster

#########################
## Visualizing results ##
#########################

library(devtools)
source_url("https://raw.githubusercontent.com/jevansbio/igraphhack/master/igraphplot2.R")
environment(plot.igraph2) <- asNamespace('igraph')
environment(igraph.Arrows2) <- asNamespace('igraph')

## Cell differentiation network
legend_image <- as.raster(matrix(colorRampPalette(c("lightgray", "red", "black"))(99), ncol=1))
pdf(file = paste(getwd(), "f", f_NA, "_diffNet.pdf", sep = ""), width = 5, height = 5)
layout(mat = matrix(c(1,1,1,2), ncol = 1))
par(mar = c(0,0,3,0))
get.cdn(res.fit = res.fit,
        edges.lab = F)
plot(c(0,1),c(-1,1),type = 'n', axes = F,xlab = '', ylab = '')
text(x=seq(0,1,l=5), y = -.2, labels = seq(0,1,l=5), cex = 2, font = 2)
rasterImage(t(legend_image), 0, 0, 1, 1)
dev.off()
}
}
