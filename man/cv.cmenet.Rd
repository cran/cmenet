\name{cv.cmenet}
\alias{cv.cmenet}
\docType{package}
\title{
  Bi-level selection of conditional main effects
}
\description{
  The main function in this package. \code{cv.cmenet} performs variable selection of conditional main effects (CMEs) via a bi-level penalization framework, with penalty parameters tuned via cross-validation.
}
\usage{
  cv.cmenet(xme, xcme, y,
            nfolds = 10, var.names = NULL,
            nlambda.sib=20, nlambda.cou=20, lambda.min.ratio=1e-6,
            ngamma=10, max.gamma=150, ntau=20, max.tau=0.01,
            tau.min.ratio=ifelse(length(y)/(ncol(xme)+ncol(xcme))<0.25,0.002,0.01),
            it.max=10, it.warm = 3, warm.str="lasso")
}
\arguments{
  \item{xme}{An \eqn{n} x \eqn{p} binary model matrix for MEs.}
  \item{xcme}{An \eqn{n} x (4*choose(\eqn{p},2)) model matrix for CMEs.}
  \item{y}{An \eqn{n}-length response vector.}
  \item{nfolds}{Number of folds for cross-validation.}
  \item{var.names}{A (\eqn{p}+4*choose(\eqn{p},2))-length string vector for variable names.}
  \item{nlambda.sib}{Number of values for \code{lambda.sib}.}
  \item{nlambda.cou}{Number of values for \code{lambda.cou}.}
  \item{lambda.min.ratio}{Smallest value for \code{lambda.sib} and \code{lambda.cou}, as a fraction of \code{lambda.max} (the smallest penalty for which all coefficients are zero).}
  \item{ngamma}{Number of values for \code{gamma}.}
  \item{max.gamma}{Maximum value for \code{gamma}.}
  \item{ntau}{Number of values for \code{tau}.}
  \item{max.tau}{Maximum value for \code{tau}.}
  \item{tau.min.ratio}{Smallest value for \code{tau}, as a fraction of \code{max.tau}.}
  \item{it.max}{Number of iterations after warm start.}
  \item{it.warm}{Number of iterations in warm start.}
  \item{warm.str}{A string (currently \code{"lasso"}, \code{"hierNet"} or NULL) indicating which method (if any) should be used for initializing active variables in optimization.}
}
\value{
  \item{cme.fit}{Fitted \code{cmenet} object.}
  \item{params}{Fitted penalty parameters (\code{lambda.sib}, \code{lambda.cou}, \code{gamma} and \code{tau}).}
  \item{select.names}{Selected ME and CME variables.}
  \item{select.idx}{Indices for \code{select.names}.}
}
\references{
Mak and Wu (2018). cmenet: a new method for bi-level variable selection of conditional main effects. \emph{Journal of the American Statistical Association}, to appear.
}
\examples{
  \donttest{
    library(MASS)
    n <- 50 #number of observations
    p <- 50 #number of main effects

    ## Simulate model matrix for MEs and CMEs
    set.seed(1000)
    rho <- 0.5 #correlation
    ones <- matrix(1,p,p)
    covmtx <- rho*ones+(1-rho)*diag(p)
    latmtx <- mvrnorm(n,p,mu=rep(0,p),Sigma=covmtx) #equicorrelated cov. matrix
    memtx <- (latmtx>=0)-(latmtx<0) #simulate model matrix for MEs
    model.mtx <- full.model.mtx(memtx)$model.mtx #generate model matrix for MEs and CMEs

    ## Set true model and generate response
    num.act <- 2 # two siblings active
    num.grp <- 4 # ... within four active groups
    ind <- c()
    for (ii in 1:num.grp){
      eff <- sample(seq(2*(p-1)),num.act)
      ind <- c(ind, p + eff + (ii-1)*(2*(p-1)))
    }
    colnames(model.mtx)[ind] # active CMEs

    des.mtx <- model.mtx[,ind]
    inter <- 12 #intercept
    xbtrue <- inter + rowSums(des.mtx)
    y <- xbtrue + rnorm(n,sd=1) #response
    xme <- model.mtx[,1:p]
    xcme <- model.mtx[,(p+1):ncol(model.mtx)]

    #---------------------------------------------------------------
    # Selection of MEs and CMEs:
    #---------------------------------------------------------------

    ## cmenet (parameters tuned via cross-validation)
    cv.cme <- cv.cmenet(xme, xcme, y, var.names=colnames(model.mtx))
    fit.cme <- cv.cme$cme.fit
    sel.cme <- cv.cme$select.idx
    colnames(model.mtx)[ind] #true model
    colnames(model.mtx)[sel.cme] #selected effects from cmenet
    colnames(model.mtx)[setdiff(sel.cme,ind)] #selected effects not in true model
    colnames(model.mtx)[setdiff(ind,sel.cme)] #true effects not in selected model

    ## lasso
    library(glmnet)
    cv.las <- cv.glmnet(cbind(xme,xcme),y)
    fit.las <- glmnet(cbind(xme,xcme),y)
    sel.las <- which(fit.las$beta[,which(cv.las$lambda==cv.las$lambda.min)]!=0)
    colnames(model.mtx)[ind] #true model
    colnames(model.mtx)[sel.las] #selected effects from lasso
    colnames(model.mtx)[setdiff(sel.las,ind)] #selected effects not in true model
    colnames(model.mtx)[setdiff(ind,sel.las)] #true effects not in selected model

    # ## sparsenet (load from CRAN archive)
    # library(sparsenet)
    # cv.sn <- cv.sparsenet(cbind(xme,xcme),y)
    # fit.sn <- sparsenet(cbind(xme,xcme),y)
    # sel.sn <- which(fit.sn$coefficients[[cv.sn$which.min[2]]]$beta[,cv.sn$which.min[1]]!=0)
    # colnames(model.mtx)[ind] #true model
    # colnames(model.mtx)[sel.sn] #selected effects from sparsenet
    # colnames(model.mtx)[setdiff(sel.sn,ind)] #selected effects not in true model
    # colnames(model.mtx)[setdiff(ind,sel.sn)] #true effects not in selected model

    #---------------------------------------------------------------
    ## Comparison:
    #---------------------------------------------------------------

    ## (a) Misspecifications
    length(setdiff(sel.cme,ind)) + length(setdiff(ind,sel.cme)) # cmenet:    14
    length(setdiff(sel.las,ind)) + length(setdiff(ind,sel.las)) # lasso:     32
    # length(setdiff(sel.sn,ind)) + length(setdiff(ind,sel.sn))   # sparsenet: 31

    ## (b) False discovery rate
    length(setdiff(sel.cme,ind))/length(sel.cme) # cmenet:    0.65
    length(setdiff(sel.las,ind))/length(sel.las) # lasso:     0.82
    # length(setdiff(sel.sn,ind))/length(sel.sn)   # sparsenet: 0.81

    ## (c) MSPE
    set.seed(1000)
    ntst <- 50
    latmtx <- mvrnorm(ntst,p,mu=rep(0,p),Sigma=covmtx)
    memtx <- (latmtx>=0)-(latmtx<0)
    tst.mtx <- full.model.mtx(memtx)$model.mtx
    xbtst <- inter + rowSums(tst.mtx[,ind])
    ytst <- xbtst + rnorm(ntst,sd=1)

    pred.cme <- predictcme(fit.cme,newx=tst.mtx)
    pred.cme <- pred.cme[,which(cv.cme$lambda.sib==cv.cme$params[1]),
                which(cv.cme$lambda.cou==cv.cme$params[2])]
    pred.las <- predict(fit.las,newx=tst.mtx)[,which(cv.las$lambda==cv.las$lambda.min)]
    # pred.sn <- predict(fit.sn,newx=tst.mtx)
    # pred.sn <- pred.sn[[which(fit.sn$gamma==cv.sn$parms.min[1])]][,
    #             which(fit.sn$lambda==cv.sn$parms.min[2])]
    mean( (ytst-pred.cme)^2 ) # cmenet:    1.67
    mean( (ytst-pred.las)^2 ) # lasso:     1.85
    # mean( (ytst-pred.sn)^2 )  # sparsenet: 1.81
  }
}
