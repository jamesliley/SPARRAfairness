% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SPARRAfairness_functions.R
\name{counterfactual_yhat}
\alias{counterfactual_yhat}
\title{counterfactual_yhat}
\usage{
counterfactual_yhat(dat, X, x = NULL, G, g, gdash, excl = NULL, n = NULL)
}
\arguments{
\item{dat}{data frame containing variables in X, Yhat, G, excl. Variables in U are assumed to be colnames(dat)\(X U Yhat U G U excl)}

\item{X}{set of variables and values on which to 'condition'; essentially we assume that any causal pathway from G to Yhat is through X}

\item{x}{values of variables X on which to condition; can be null in which case we use marginal distribution of X}

\item{G}{grouping variable, usually a sensitive attribute}

\item{g}{conditioned value of g}

\item{gdash}{counterfactual value of g}

\item{excl}{variable names to exclude from U}

\item{n}{number of samples; if NULL return all}
}
\value{
indices representing sample(s) from counterfactual  Yhat_{g' <- G} |X=x,G=g
}
\description{
Estimation of counterfactual quantities by resampling.
}
\details{
Counterfactual fairness is with respect to the causal graph:

G <----------U
|         /   |
|     /       |
| /           |
VV            V
X -------->  Yhat

where

G=group (usually sensitive attribute);  Yhat=outcome; X=set of variables
through which G can act on Yhat, U=set of background variables;

We want the counterfactual Yhat_{g' <- G} |X=x,G=g (or alternatively
Yhat_{g' <- G} |G=g).

This can be interpreted as:
\verb{the distribution of values of Yhat amongst indivdiduals whose values of U  are distributed as though they were in group G=g (and, optionally, had  values X=x, but whose value of G is g'}

Essentially, comparison of the counterfactual quantity above to the
conditional Yhat|G=g isolates the difference in Yhat due to the effect of
G on Yhat through X, removing any effect due to different distributions of
U due to different values of G.

To estimate Y'=Yhat_{g' <- G} | G=g, we need to
\enumerate{
\item Compute U'~(U|G=g)
\item Compute the distribution X' as X'~(X|U~U',G=g')
\item Sample Y'~(Yhat|X~X',U~U')
}

To estimate Y'=Yhat_{g' <- G} |X=x, G=g, we need to
\enumerate{
\item Compute U'~(U|G=g,X=x)
\item Compute the distribution X' as X'~(X|U~U',G=g')
\item Sample Y'~(Yhat|X~X',U~U')
}

This function approximates this samplying procedure as follows
\enumerate{
\item Look at individuals with G=g (and optionally X=x)
\item Find the values of U for these individuals
\item Find a second set of individuals with the same values of U but for whom G=g'
\item Return the indices of these individuals
}

The values of Yhat for these individuals constitute a sample from the
desired counterfactual.
}
\examples{

set.seed(23173)
N=10000

# Background variables sampler
background_U=function(n) runif(n) # U~U(0,1)

# Structural equations
struct_G=function(u,n) rbinom(n,1,prob=u) # G|U=u ~ Bern(u)
struct_X=function(u,g,n) rbinom(n,1,prob=u*(0.5 + 0.5*g)) # X|U=u,G=g ~ Bern(u(1+g)/2)
struct_Yhat=function(u,x,n) (runif(n,0,x) + runif(n,0,u))/2 # Yhat|X,N ~ (U(0,X) + U(0,U))/2

# To see that the counterfactual 'isolates' the difference in Yhat due to the 
#  causal pathway from G to Yhat through X, change the definition of struct_G to
#  
#  struct_G=function(u,n) rbinom(n,1,prob=1/2) # G|U=u ~ Bern(1/2)
#
# so the posterior of U|G=g does not depend on g. Note that, with this definition, the 
#  counterfactual Yhat_{G<01}|G=1 coincides with the conditional Yhat|G=0, since
#  the counterfactual G<-1 is equivalent to just conditioning on G=1.
#
# By contrast, if we change struct_G back to its original definition, but 
#  change the definition of struct_Yhat to
#  
#  struct_Yhat=function(u,x,n) (runif(n,0,1) + runif(n,0,u))/2 # Yhat|X,N ~ (U(0,1) + U(0,U))/2
#
# so Yhat depends on G only through the change in posterior of U from changing g, 
#  the counterfactual Yhat_{G<01}|G=1 coincides with the conditional Yhat|G=1.


# Sample from complete causal model
U=background_U(N)
G=struct_G(U,N)
X=struct_X(U,G,N)
Yhat=struct_Yhat(U,X,N)
dat=data.frame(U,G,X,Yhat)


# True counterfactual Yhat_{G <- 0}|G=1
w1=which(dat$G==1)
n1=length(w1)
UG1=dat$U[w1] # This is U|G=1
XG1=struct_X(UG1,rep(0,n1),n1)
YhatG1=struct_Yhat(UG1,XG1,n1)

# Estimated counterfactual Yhat_{G <- 0}|G=1
ind_G1=counterfactual_yhat(dat,X="X",G="G",g = 1, gdash = 0)
YhatG1_resample=dat$Yhat[ind_G1]



# True counterfactual Yhat_{G <- 0}|G=1,X=1
w11=which(dat$G==1 & dat$X==1)
n11=length(w11)
UG1X1=dat$U[w11] # This is U|G=1,X=1
XG1X1=struct_X(UG1X1,rep(0,n11),n11)
YhatG1X1=struct_Yhat(UG1X1,XG1X1,n11)

# Estimated counterfactual Yhat_{G <- 0}|G=1
ind_G1X1=counterfactual_yhat(dat,X="X",G="G",g = 1, gdash = 0,x=1)
YhatG1X1_resample=dat$Yhat[ind_G1X1]




# Compare CDFs
x=seq(0,1,length=1000)
par(mfrow=c(1,2))

plot(0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="Value",
     ylab=expression(paste("Prop. ",hat('Y')," < x")))
lines(x,ecdf(dat$Yhat)(x),col="black") # Unconditional CDF of Yhat
lines(x,ecdf(dat$Yhat[which(dat$G==1)])(x),col="red") # Yhat|G=1
lines(x,ecdf(dat$Yhat[which(dat$G==0)])(x),col="blue") # Yhat|G=0

# True counterfactual Yhat_{G <- 0}|G=1
lines(x,ecdf(YhatG1)(x),col="blue",lty=2) 

# Estimated counterfactual Yhat_{G <- 0}|G=1
lines(x,ecdf(YhatG1_resample)(x),col="blue",lty=3) 

legend("bottomright",
       c(expression(paste(hat('Y'))),
         expression(paste(hat('Y'),"|G=1")),
         expression(paste(hat('Y'),"|G=0")),
         expression(paste(hat(Y)[G \%<-\% 0],"|G=1 (true)")),
         expression(paste(hat(Y)[G \%<-\% 0],"|G=1 (est.)"))),
       col=c("black","red","blue","blue","blue"),
       lty=c(1,1,1,2,3),
       cex=0.5)



plot(0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="Value",
     ylab=expression(paste("Prop. ",hat('Y')," < x")))
lines(x,ecdf(dat$Yhat[which(dat$X==1)])(x),col="black") # CDF of Yhat|X=1
lines(x,ecdf(dat$Yhat[which(dat$G==1 & dat$X==1)])(x),col="red") # Yhat|G=1,X=1
lines(x,ecdf(dat$Yhat[which(dat$G==0 & dat$X==1)])(x),col="blue") # Yhat|G=0,X=1

# True counterfactual Yhat_{G <- 0}|G=1,X=1
lines(x,ecdf(YhatG1X1)(x),col="blue",lty=2) 

# Estimated counterfactual Yhat_{G <- 0}|G=1,X=1
lines(x,ecdf(YhatG1X1_resample)(x),col="blue",lty=3) 

legend("bottomright",
       c(expression(paste(hat('Y|X=1'))),
         expression(paste(hat('Y'),"|G=1,X=1")),
         expression(paste(hat('Y'),"|G=0,X=1")),
         expression(paste(hat(Y)[G \%<-\% 0],"|G=1,X=1 (true)")),
         expression(paste(hat(Y)[G \%<-\% 0],"|G=1,X=1 (est.)"))),
       col=c("black","red","blue","blue","blue"),
       lty=c(1,1,1,2,3),
       cex=0.5)

# In both plots, the estimated counterfactual CDF closely matches the CDF of the
#  true counterfactual.


}
