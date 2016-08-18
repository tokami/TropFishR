# add description
# change name?
# This functions allows to estimate inital values for optim functinos for 4 of 7 rtypes
# This is very convient as guessing inital values can be hard, the method is very sensitive to these values
# and for each rtypes different intital values are required
# this is based on the gillnet functions from Millar (gillnet select)
# add these options for trawl selectivities

find_init_Millar <- function(param, rtype="norm.loc", rel.power=NULL){
  # Adapted R code from Russell Millar (https://www.stat.auckland.ac.nz/~millar/selectware/)

  res <- param
  meshSizes <- res$meshSizes
  classes <- res$midLengths
  CatchPerNet_mat <- res$CatchPerNet_mat
  lens <- rep(classes,ncol(CatchPerNet_mat))
  msizes <- rep(meshSizes,rep(length(classes),ncol(CatchPerNet_mat)))
  dat <- as.vector(CatchPerNet_mat)

  var1 <- lens * msizes
  var2 <- msizes^2
  var3 <- (lens / msizes)^2
  var4 <- lens / msizes
  var5 <- -log(msizes)
  var6 <- log(msizes / msizes[1])
  var7 <- var6 * log(lens) - 0.5 * var6 * var6
  var8 <- lens * lens
  var9 <- msizes / lens

  if(is.null(rel.power)){
    os <- 0
  }else{
    os <- rep(log(rel.power), rep(length(classes), ncol(CatchPerNet_mat)))
  }

  switch(rtype,

         "norm.loc"={
           if(is.null(rel.power))
             fit <- glm(dat ~ -1 + var1 + var2 + as.factor(lens),family=poisson)
           else
             fit <- glm(dat ~ -1 + var1 + var2 + as.factor(lens) + offset(os),family=poisson)
           x <- coef(fit)[c("var1","var2")]
           varx <- summary(fit)$cov.unscaled[1:2,1:2]
           k <- -2*x[2]/x[1]
           sigma <- sqrt(-2*x[2]/(x[1]^2))
           pars <- c(k*msizes[1],sigma) #c("mode(mesh1)","std_dev(all meshes)")
           },

         "norm.sca"={
           if(missing(rel.power))
             fit <- glm(dat ~ -1 + var3 + var4 + as.factor(lens),family=poisson)
           else
             fit <- glm(dat ~ -1 + var3 + var4 + as.factor(lens) + offset(os),family=poisson)
           x <- coef(fit)[c("var3","var4")]
           varx <- summary(fit)$cov.unscaled[1:2,1:2]
           k1 <- -x[2]/(2*x[1])
           k2 <- -1/(2*x[1])
           pars <- c(k1*msizes[1],sqrt(k2*msizes[1]^2)) #c("mode(mesh1)","std_dev(mesh1)")
           },

         "gamma"={
           if(is.null(rel.power))
             fit <- glm(dat ~ -1 + var4 + var5 + as.factor(lens), family = poisson)
           else
             fit <- glm(dat ~ -1 + var4 + var5 + as.factor(lens) + offset(os), family = poisson)
           x <- coef(fit)[c("var4","var5")]
           varx <- summary(fit)$cov.unscaled[1:2,1:2]
           alpha <- x[2] + 1
           k <- -1 / x[1]
           pars <- c(x[1],x[2],alpha,k,(alpha-1)*k*msizes[1],sqrt(alpha*(k*msizes[1])^2)) #c("mode(mesh1)","std_dev(mesh1)")
           },

         "lognorm"={
           if(is.null(rel.power))
             fit <- glm(dat ~ -1 + var6 + var7 + as.factor(lens),family=poisson)
           else
             fit <- glm(dat ~ -1 + var6 + var7 + as.factor(lens) + offset(os),family=poisson)
           x <- coef(fit)[c("var6","var7")]
           varx <- summary(fit)$cov.unscaled[1:2,1:2]
           mu1 <- -(x[1]-1)/x[2]
           sigma <- sqrt(1/x[2])
           pars <- c(exp(mu1-sigma^2),
                     sqrt(exp(2*mu1+sigma^2)*(exp(sigma^2)-1))) #c("mode(mesh1)","std_dev(mesh1)")
           },

         stop(paste("\n For ",rtype, " it is not yet possible to estimate initial parameters.",
                    "\n Please enter initial parameters (x0).")))
  return(pars)
}

