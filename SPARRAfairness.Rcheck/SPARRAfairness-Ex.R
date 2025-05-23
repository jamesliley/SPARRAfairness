pkgname <- "SPARRAfairness"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('SPARRAfairness')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("adjusted_fdr")
### * adjusted_fdr

flush(stderr()); flush(stdout())

### Name: adjusted_fdr
### Title: adjusted_fdr
### Aliases: adjusted_fdr

### ** Examples

# See vignette



cleanEx()
nameEx("adjusted_for")
### * adjusted_for

flush(stderr()); flush(stdout())

### Name: adjusted_for
### Title: adjusted_for
### Aliases: adjusted_for

### ** Examples

# See vignette



cleanEx()
nameEx("build_diff")
### * build_diff

flush(stderr()); flush(stdout())

### Name: build_diff
### Title: build_diff Prepares a data frame for a ggplot object to compare
###   differences using linear interpolation.
### Aliases: build_diff

### ** Examples

# Only used internally



cleanEx()
nameEx("cal_2panel")
### * cal_2panel

flush(stderr()); flush(stdout())

### Name: cal_2panel
### Title: cal_2panel Draws calibration curves (with legend) with a second
###   panel underneath showing predicted differences.
### Aliases: cal_2panel

### ** Examples

# See vignette



cleanEx()
nameEx("counterfactual_yhat")
### * counterfactual_yhat

flush(stderr()); flush(stdout())

### Name: counterfactual_yhat
### Title: counterfactual_yhat
### Aliases: counterfactual_yhat

### ** Examples


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
#  counterfactual Yhat{G<-1}|G=1 coincides with the conditional Yhat|G=0, since
#  the counterfactual G<-1 is equivalent to just conditioning on G=1.
#
# By contrast, if we change struct_G back to its original definition, but 
#  change the definition of struct_Yhat to
#  
#  struct_Yhat=function(u,x,n) (runif(n,0,1) + runif(n,0,u))/2 # Yhat|X,N ~ (U(0,1) + U(0,U))/2
#
# so Yhat depends on G only through the change in posterior of U from changing g, 
#  the counterfactual Yhat{G<01}|G=1 coincides with the conditional Yhat|G=1.


# Sample from complete causal model
U=background_U(N)
G=struct_G(U,N)
X=struct_X(U,G,N)
Yhat=struct_Yhat(U,X,N)
dat=data.frame(U,G,X,Yhat)


# True counterfactual Yhat{G <- 0}|G=1
w1=which(dat$G==1)
n1=length(w1)
UG1=dat$U[w1] # This is U|G=1
XG1=struct_X(UG1,rep(0,n1),n1)
YhatG1=struct_Yhat(UG1,XG1,n1)

# Estimated counterfactual Yhat{G <- 0}|G=1
ind_G1=counterfactual_yhat(dat,X="X",G="G",g = 1, gdash = 0)
YhatG1_resample=dat$Yhat[ind_G1]



# True counterfactual Yhat{G <- 0}|G=1,X=1
w11=which(dat$G==1 & dat$X==1)
n11=length(w11)
UG1X1=dat$U[w11] # This is U|G=1,X=1
XG1X1=struct_X(UG1X1,rep(0,n11),n11)
YhatG1X1=struct_Yhat(UG1X1,XG1X1,n11)

# Estimated counterfactual Yhat{G <- 0}|G=1
ind_G1X1=counterfactual_yhat(dat,X="X",G="G",g = 1, gdash = 0,x=1)
YhatG1X1_resample=dat$Yhat[ind_G1X1]




# Compare CDFs
x=seq(0,1,length=1000)
oldpar = par(mfrow=c(1,2))

plot(0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="Value",
     ylab=expression(paste("Prop. ",hat('Y')," < x")))
lines(x,ecdf(dat$Yhat)(x),col="black") # Unconditional CDF of Yhat
lines(x,ecdf(dat$Yhat[which(dat$G==1)])(x),col="red") # Yhat|G=1
lines(x,ecdf(dat$Yhat[which(dat$G==0)])(x),col="blue") # Yhat|G=0

# True counterfactual Yhat{G <- 0}|G=1
lines(x,ecdf(YhatG1)(x),col="blue",lty=2) 

# Estimated counterfactual Yhat{G <- 0}|G=1
lines(x,ecdf(YhatG1_resample)(x),col="blue",lty=3) 

legend("bottomright",
       c(expression(paste(hat('Y'))),
         expression(paste(hat('Y'),"|G=1")),
         expression(paste(hat('Y'),"|G=0")),
         expression(paste(hat(Y)[G %<-% 0],"|G=1 (true)")),
         expression(paste(hat(Y)[G %<-% 0],"|G=1 (est.)"))),
       col=c("black","red","blue","blue","blue"),
       lty=c(1,1,1,2,3),
       cex=0.5)



plot(0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="Value",
     ylab=expression(paste("Prop. ",hat('Y')," < x")))
lines(x,ecdf(dat$Yhat[which(dat$X==1)])(x),col="black") # CDF of Yhat|X=1
lines(x,ecdf(dat$Yhat[which(dat$G==1 & dat$X==1)])(x),col="red") # Yhat|G=1,X=1
lines(x,ecdf(dat$Yhat[which(dat$G==0 & dat$X==1)])(x),col="blue") # Yhat|G=0,X=1

# True counterfactual Yhat{G <- 0}|G=1,X=1
lines(x,ecdf(YhatG1X1)(x),col="blue",lty=2) 

# Estimated counterfactual Yhat{G <- 0}|G=1,X=1
lines(x,ecdf(YhatG1X1_resample)(x),col="blue",lty=3) 

legend("bottomright",
       c(expression(paste(hat('Y|X=1'))),
         expression(paste(hat('Y'),"|G=1,X=1")),
         expression(paste(hat('Y'),"|G=0,X=1")),
         expression(paste(hat(Y)[G %<-% 0],"|G=1,X=1 (true)")),
         expression(paste(hat(Y)[G %<-% 0],"|G=1,X=1 (est.)"))),
       col=c("black","red","blue","blue","blue"),
       lty=c(1,1,1,2,3),
       cex=0.5)

# In both plots, the estimated counterfactual CDF closely matches the CDF of the
#  true counterfactual.

# Restore parameters
par(oldpar)




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("dat2mat")
### * dat2mat

flush(stderr()); flush(stdout())

### Name: dat2mat
### Title: dat2mat
### Aliases: dat2mat

### ** Examples

# See vignette



cleanEx()
nameEx("demographic_parity")
### * demographic_parity

flush(stderr()); flush(stdout())

### Name: demographic_parity
### Title: demographic_parity
### Aliases: demographic_parity

### ** Examples

# See vignette



cleanEx()
nameEx("drawperson")
### * drawperson

flush(stderr()); flush(stdout())

### Name: drawperson
### Title: drawperson
### Aliases: drawperson

### ** Examples

plot(0,xlim=c(-1,1),ylim=c(-1,1),type="n")
drawperson(0,0,1,col="yellow",border="red",lwd=3,lty=2)



cleanEx()
nameEx("drawprop")
### * drawprop

flush(stderr()); flush(stdout())

### Name: drawprop
### Title: drawprop
### Aliases: drawprop

### ** Examples

# See vignette



cleanEx()
nameEx("for_breakdown")
### * for_breakdown

flush(stderr()); flush(stdout())

### Name: for_breakdown
### Title: for_breakdown
### Aliases: for_breakdown

### ** Examples


# See vignette



cleanEx()
nameEx("getcal")
### * getcal

flush(stderr()); flush(stdout())

### Name: getcal
### Title: getcal()
### Aliases: getcal

### ** Examples

# See vignette



cleanEx()
nameEx("getprc")
### * getprc

flush(stderr()); flush(stdout())

### Name: getprc
### Title: getprc()
### Aliases: getprc

### ** Examples

# See vignette



cleanEx()
nameEx("getroc")
### * getroc

flush(stderr()); flush(stdout())

### Name: getroc
### Title: getroc() Comprehensive plotting function for receiver-operator
###   characteristic curve. Also calculates AUROC and standard error.
### Aliases: getroc

### ** Examples

# See vignette



cleanEx()
nameEx("group_fairness")
### * group_fairness

flush(stderr()); flush(stdout())

### Name: group_fairness
### Title: group_fairness
### Aliases: group_fairness

### ** Examples

# See vignette



cleanEx()
nameEx("groupmetric_2panel")
### * groupmetric_2panel

flush(stderr()); flush(stdout())

### Name: groupmetric_2panel
### Title: groupmetric_2panel Draws plots of a group fairness metric with a
###   second panel underneath
### Aliases: groupmetric_2panel

### ** Examples

# See vignette



cleanEx()
nameEx("logistic")
### * logistic

flush(stderr()); flush(stdout())

### Name: logistic
### Title: Logistic
### Aliases: logistic

### ** Examples


# Plot
x=seq(-5,5,length=1000)
plot(x,logistic(x),type="l")



cleanEx()
nameEx("logit")
### * logit

flush(stderr()); flush(stdout())

### Name: logit
### Title: Logit
### Aliases: logit

### ** Examples


# Plot
x=seq(0,1,length=100)
plot(x,logit(x),type="l")

# Logit and logistic are inverses
x=seq(-5,5,length=1000)
plot(x,logit(logistic(x)),type="l")



cleanEx()
nameEx("plot.sparraCAL")
### * plot.sparraCAL

flush(stderr()); flush(stdout())

### Name: plot.sparraCAL
### Title: Plot function for class sparraCAL
### Aliases: plot.sparraCAL

### ** Examples

# See vignette



cleanEx()
nameEx("plot.sparraPRC")
### * plot.sparraPRC

flush(stderr()); flush(stdout())

### Name: plot.sparraPRC
### Title: Plot function for class above
### Aliases: plot.sparraPRC

### ** Examples

# See vignette



cleanEx()
nameEx("plot.sparraROC")
### * plot.sparraROC

flush(stderr()); flush(stdout())

### Name: plot.sparraROC
### Title: Plot function for class sparraROC
### Aliases: plot.sparraROC

### ** Examples

# See vignette



cleanEx()
nameEx("plot_decomp")
### * plot_decomp

flush(stderr()); flush(stdout())

### Name: plot_decomp
### Title: plot_decomp
### Aliases: plot_decomp

### ** Examples


# See vignette



cleanEx()
nameEx("prc_2panel")
### * prc_2panel

flush(stderr()); flush(stdout())

### Name: prc_2panel
### Title: prc_2panel Draws a PRC curve (with legend) with a second panel
###   underneath showing precision difference.
### Aliases: prc_2panel

### ** Examples

# See vignette



cleanEx()
nameEx("roc_2panel")
### * roc_2panel

flush(stderr()); flush(stdout())

### Name: roc_2panel
### Title: roc_2panel Draws a ROC curve (with legend) with a second panel
###   underneath showing sensitivity difference.
### Aliases: roc_2panel

### ** Examples

# See vignette



cleanEx()
nameEx("sim_pop_data")
### * sim_pop_data

flush(stderr()); flush(stdout())

### Name: sim_pop_data
### Title: sim_pop_data
### Aliases: sim_pop_data

### ** Examples


# Simulate data
dat=sim_pop_data(10000)
cor(dat[,1:7])

# See vignette



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
