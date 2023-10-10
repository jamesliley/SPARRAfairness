##**************************************************#
## R script for functions in SPARRAfairness      ####
##**************************************************#
##
## James Liley, Ioanna Thoma
## March 2023
##

##**************************************************#
## Colour schemes                                ####
##**************************************************#

##' phs_colours
##' 
##' Copied from github, "Public-Health-Scotland/phsstyles". Public Health Scotland
##'  colour scheme. Internal function.
##' 
##' @param colourname name of colour; usually something like phs-blue. If NULL returns all colours.
##' @param keep_names keep names of colours in return list. Defaults to false.
##' @export
##' @return vector of colours, optionally with names.
phs_colours=function (colourname = NULL, keep_names = FALSE) {
  phs_colour_values <- c(`phs-purple` = "#3F3685", `phs-magenta` = "#9B4393", `phs-blue` = "#0078D4", 
                         `phs-green` = "#83BB26", `phs-graphite` = "#948DA3", `phs-teal` = "#1E7F84", 
                         `phs-liberty` = "#6B5C85", `phs-rust` = "#C73918", `phs-purple-80` = "#655E9D", 
                         `phs-purple-50` = "#9F9BC2", `phs-purple-30` = "#C5C3DA", `phs-purple-10` = "#ECEBF3", 
                         `phs-magenta-80` = "#AF69A9", `phs-magenta-50` = "#CDA1C9", `phs-magenta-30` = "#E1C7DF", 
                         `phs-magenta-10` = "#F5ECF4", `phs-blue-80` = "#3393DD", `phs-blue-50` = "#80BCEA", 
                         `phs-blue-30` = "#B3D7F2", `phs-blue-10` = "#E6F2FB", `phs-green-80` = "#9CC951", 
                         `phs-green-50` = "#C1DD93", `phs-green-30` = "#DAEBBE", `phs-green-10` = "#F3F8E9", 
                         `phs-graphite-80` = "#A9A4B5", `phs-graphite-50` = "#CAC6D1", 
                         `phs-graphite-30` = "#DFDDE3", `phs-graphite-10` = "#F4F4F6", 
                         `phs-teal-80` = "#4B999D", `phs-teal-50` = "#8FBFC2", `phs-teal-30` = "#BCD9DA", 
                         `phs-teal-10` = "#E9F2F3", `phs-liberty-80` = "#897D9D", `phs-liberty-50` = "#B5AEC2", 
                         `phs-liberty-30` = "#D3CEDA", `phs-liberty-10` = "#F0EFF3", `phs-rust-80` = "#D26146", 
                         `phs-rust-50` = "#E39C8C", `phs-rust-30` = "#EEC4BA", `phs-rust-10` = "#F9EBE8"
  )
  if (is.null(colourname)) {
    phs_colour_values
  }
  else if (any(!colourname %in% names(phs_colour_values))) {
    col_not_list <- colourname[!colourname %in% names(phs_colour_values)]
    msg <- paste("These colours are not available:", paste(col_not_list, 
                                                           collapse = ","), "\nPlease run phs_colours() to see all the available colours")
    stop(msg)
  }
  else if (!keep_names) {
    unname(phs_colour_values[colourname])
  }
  else {
    phs_colour_values[colourname]
  }
}


##**************************************************#
## Small auxiliary functions                     ####
##**************************************************#

##' Logit
##'
##' @export
##' @name logit
##' @description Logit function: -log((1/x)-1)
##' @param x argument
##' @return value of logit(x); na if x is outside (0,1)
##' @examples
##'
##' # Plot
##' x=seq(0,1,length=100)
##' plot(x,logit(x),type="l")
##'
##' # Logit and logistic are inverses
##' x=seq(-5,5,length=1000)
##' plot(x,logit(logistic(x)),type="l")
logit=function(x) -log((1/x)-1)


##' Logistic
##'
##' @export
##' @name logistic
##' @description Logistic function: 1/(1+exp(-x))
##' @param x argument
##' @return value of logistic(x)
##' @examples
##'
##' # Plot
##' x=seq(-5,5,length=1000)
##' plot(x,logistic(x),type="l")
logistic=function(x) 1/(1+exp(-x))




##' integral() 
##' Quick form for trapezoidal integration over range of x
##'
##' @param x x co-ordinates, or nx2 matrix of points 
##' @param y y co-ordinates
##' @export 
##' @return trapezoidal estimate of integral of the xth value of y over range of x.
integral=function(x,y=NULL) {
  if (is.null(y)) {
    y=x[,2]; x=x[,1]
  }
  ox=order(x); xs=x[ox]; ys=y[ox]
  sum((xs[-1]-xs[-length(xs)])*(ys[-1]+ys[-length(ys)]))/2
}

##' ab() 
##' Shorthand to draw a red x-y line
##' @param ... passed to abline()
##' @return No return value, draws a figure
##' @export
ab=function(...) abline(0,1,col=phs_colours("phs-magenta"),...)




##**************************************************#
## Counterfactuals                               ####
##**************************************************#

##' counterfactual_yhat
##' 
##' Estimation of counterfactual quantities by resampling.
##' 
##' Counterfactual fairness is with respect to the causal graph:
##'
##'
##'   G <----------U
##'   |         /   |
##'   |     /       |
##'   | /           |
##'   VV            V
##'   X -------->  Yhat
##'   
##' where
##'
##' G=group (usually sensitive attribute);  Yhat=outcome; X=set of variables 
##'  through which G can act on Yhat, U=set of background variables;
##'
##' We want the counterfactual Yhat_{g' <- G} |X=x,G=g (or alternatively 
##'  Yhat_{g' <- G} |G=g).
##' 
##' This can be interpreted as:
##'  `the distribution of values of Yhat amongst indivdiduals whose values of U 
##'   are distributed as though they were in group G=g (and, optionally, had 
##'   values X=x, but whose value of G is g'`
##' 
##' Essentially, comparison of the counterfactual quantity above to the 
##'  conditional Yhat|G=g isolates the difference in Yhat due to the effect of 
##'  G on Yhat through X, removing any effect due to different distributions of 
##'  U due to different values of G.
##' 
##' To estimate Y'=Yhat_{g' <- G} | G=g, we need to
##'  1. Compute U'~(U|G=g)
##'  2. Compute the distribution X' as X'~(X|U~U',G=g')
##'  3. Sample Y'~(Yhat|X~X',U~U')
##' 
##' To estimate Y'=Yhat_{g' <- G} |X=x, G=g, we need to
##'  1. Compute U'~(U|G=g,X=x)
##'  2. Compute the distribution X' as X'~(X|U~U',G=g')
##'  3. Sample Y'~(Yhat|X~X',U~U')
##'  
##' This function approximates this samplying procedure as follows
##'  1. Look at individuals with G=g (and optionally X=x)
##'  2. Find the values of U for these individuals
##'  3. Find a second set of individuals with the same values of U but for whom G=g'
##'  4. Return the indices of these individuals
##'  
##' The values of Yhat for these individuals constitute a sample from the 
##'  desired counterfactual. 
##'  
##' @param dat data frame containing variables in X, Yhat, G, excl. Variables in U are assumed to be colnames(dat)\(X U Yhat U G U excl)
##' @param X set of variables and values on which to 'condition'; essentially we assume that any causal pathway from G to Yhat is through X
##' @param x values of variables X on which to condition; can be null in which case we use marginal distribution of X
##' @param G grouping variable, usually a sensitive attribute
##' @param g conditioned value of g
##' @param gdash counterfactual value of g
##' @param excl variable names to exclude from U
##' @param n number of samples; if NULL return all
##' @export 
##' @return indices representing sample(s) from counterfactual  Yhat_{g' <- G} |X=x,G=g
##' @examples
##' 
##' set.seed(23173)
##' N=10000
##' 
##' # Background variables sampler
##' background_U=function(n) runif(n) # U~U(0,1)
##' 
##' # Structural equations
##' struct_G=function(u,n) rbinom(n,1,prob=u) # G|U=u ~ Bern(u)
##' struct_X=function(u,g,n) rbinom(n,1,prob=u*(0.5 + 0.5*g)) # X|U=u,G=g ~ Bern(u(1+g)/2)
##' struct_Yhat=function(u,x,n) (runif(n,0,x) + runif(n,0,u))/2 # Yhat|X,N ~ (U(0,X) + U(0,U))/2
##' 
##' # To see that the counterfactual 'isolates' the difference in Yhat due to the 
##' #  causal pathway from G to Yhat through X, change the definition of struct_G to
##' #  
##' #  struct_G=function(u,n) rbinom(n,1,prob=1/2) # G|U=u ~ Bern(1/2)
##' #
##' # so the posterior of U|G=g does not depend on g. Note that, with this definition, the 
##' #  counterfactual Yhat_{G<01}|G=1 coincides with the conditional Yhat|G=0, since
##' #  the counterfactual G<-1 is equivalent to just conditioning on G=1.
##' #
##' # By contrast, if we change struct_G back to its original definition, but 
##' #  change the definition of struct_Yhat to
##' #  
##' #  struct_Yhat=function(u,x,n) (runif(n,0,1) + runif(n,0,u))/2 # Yhat|X,N ~ (U(0,1) + U(0,U))/2
##' #
##' # so Yhat depends on G only through the change in posterior of U from changing g, 
##' #  the counterfactual Yhat_{G<01}|G=1 coincides with the conditional Yhat|G=1.
##' 
##' 
##' # Sample from complete causal model
##' U=background_U(N)
##' G=struct_G(U,N)
##' X=struct_X(U,G,N)
##' Yhat=struct_Yhat(U,X,N)
##' dat=data.frame(U,G,X,Yhat)
##' 
##' 
##' # True counterfactual Yhat_{G <- 0}|G=1
##' w1=which(dat$G==1)
##' n1=length(w1)
##' UG1=dat$U[w1] # This is U|G=1
##' XG1=struct_X(UG1,rep(0,n1),n1)
##' YhatG1=struct_Yhat(UG1,XG1,n1)
##' 
##' # Estimated counterfactual Yhat_{G <- 0}|G=1
##' ind_G1=counterfactual_yhat(dat,X="X",G="G",g = 1, gdash = 0)
##' YhatG1_resample=dat$Yhat[ind_G1]
##' 
##' 
##' 
##' # True counterfactual Yhat_{G <- 0}|G=1,X=1
##' w11=which(dat$G==1 & dat$X==1)
##' n11=length(w11)
##' UG1X1=dat$U[w11] # This is U|G=1,X=1
##' XG1X1=struct_X(UG1X1,rep(0,n11),n11)
##' YhatG1X1=struct_Yhat(UG1X1,XG1X1,n11)
##' 
##' # Estimated counterfactual Yhat_{G <- 0}|G=1
##' ind_G1X1=counterfactual_yhat(dat,X="X",G="G",g = 1, gdash = 0,x=1)
##' YhatG1X1_resample=dat$Yhat[ind_G1X1]
##' 
##' 
##' 
##' 
##' # Compare CDFs
##' x=seq(0,1,length=1000)
##' oldpar = par(mfrow=c(1,2))
##' 
##' plot(0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="Value",
##'      ylab=expression(paste("Prop. ",hat('Y')," < x")))
##' lines(x,ecdf(dat$Yhat)(x),col="black") # Unconditional CDF of Yhat
##' lines(x,ecdf(dat$Yhat[which(dat$G==1)])(x),col="red") # Yhat|G=1
##' lines(x,ecdf(dat$Yhat[which(dat$G==0)])(x),col="blue") # Yhat|G=0
##' 
##' # True counterfactual Yhat_{G <- 0}|G=1
##' lines(x,ecdf(YhatG1)(x),col="blue",lty=2) 
##' 
##' # Estimated counterfactual Yhat_{G <- 0}|G=1
##' lines(x,ecdf(YhatG1_resample)(x),col="blue",lty=3) 
##' 
##' legend("bottomright",
##'        c(expression(paste(hat('Y'))),
##'          expression(paste(hat('Y'),"|G=1")),
##'          expression(paste(hat('Y'),"|G=0")),
##'          expression(paste(hat(Y)[G %<-% 0],"|G=1 (true)")),
##'          expression(paste(hat(Y)[G %<-% 0],"|G=1 (est.)"))),
##'        col=c("black","red","blue","blue","blue"),
##'        lty=c(1,1,1,2,3),
##'        cex=0.5)
##' 
##' 
##' 
##' plot(0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="Value",
##'      ylab=expression(paste("Prop. ",hat('Y')," < x")))
##' lines(x,ecdf(dat$Yhat[which(dat$X==1)])(x),col="black") # CDF of Yhat|X=1
##' lines(x,ecdf(dat$Yhat[which(dat$G==1 & dat$X==1)])(x),col="red") # Yhat|G=1,X=1
##' lines(x,ecdf(dat$Yhat[which(dat$G==0 & dat$X==1)])(x),col="blue") # Yhat|G=0,X=1
##' 
##' # True counterfactual Yhat_{G <- 0}|G=1,X=1
##' lines(x,ecdf(YhatG1X1)(x),col="blue",lty=2) 
##' 
##' # Estimated counterfactual Yhat_{G <- 0}|G=1,X=1
##' lines(x,ecdf(YhatG1X1_resample)(x),col="blue",lty=3) 
##' 
##' legend("bottomright",
##'        c(expression(paste(hat('Y|X=1'))),
##'          expression(paste(hat('Y'),"|G=1,X=1")),
##'          expression(paste(hat('Y'),"|G=0,X=1")),
##'          expression(paste(hat(Y)[G %<-% 0],"|G=1,X=1 (true)")),
##'          expression(paste(hat(Y)[G %<-% 0],"|G=1,X=1 (est.)"))),
##'        col=c("black","red","blue","blue","blue"),
##'        lty=c(1,1,1,2,3),
##'        cex=0.5)
##' 
##' # In both plots, the estimated counterfactual CDF closely matches the CDF of the
##' #  true counterfactual.
##' 
##' # Restore parameters
##' par(oldpar)
##' 
counterfactual_yhat=function(dat,X,x=NULL,G,g,gdash,excl=NULL,n=NULL) {
  
  U=setdiff(colnames(dat),c(X,G,excl))
  dat$id=1:dim(dat)[1]
  odat=dat[do.call(order,dat[,U]),]
  
  datXG=odat[,c("id",X,G)]
  w=1:dim(dat)[1]
  if (!is.null(x)) for (i in 1:length(x)) w=intersect(w,which(datXG[,X]==x[i]))
  w=intersect(w,which(datXG[,dim(datXG)[2]]==g))
  
  w2=which(datXG[[G]]==gdash)
  
  ax=approx(w2,1:length(w2),w,rule=2)$y
  wz1=w2[floor(ax)]
  wz2=w2[ceiling(ax)]
  
  w12=which(abs(wz2-w)<abs(wz1-w))
  wz=wz1; wz[w12]=wz2[w12]
  
  if (!is.null(n)) out=sample(datXG$id[wz],n) else out=datXG$id[wz]
  
  return(out)
}





##**************************************************#
## Compute curves                                ####
##**************************************************#

##' demographic_parity
##' 
##' Estimates demographic parity for a risk score (essentially cumulative distribution function)
##' 
##' @param scores vector of risk scores
##' @param group1 indices of group 1
##' @param group2 indices of group 2
##' @param cutoffs score cutoffs at which to estimate DP (default 100 evenly-spaced)
##' @export
##' @return matrix of dimension length(cutoffs)x4, with (i,2g-1)th entry the proportion of scores in group g which are less than or equal to the ith cutoff value and (i,2g)th entry the approximate standard error of the (i,2g-1)th entry
##' @examples 
##' # See vignette
demographic_parity=function(scores,group1,group2, cutoffs=seq(min(scores,na.rm=TRUE),max(scores,na.rm=TRUE),length=100)) {
  
  s1=scores[group1]
  s2=scores[group2]
  
  n1=length(group1)
  n2=length(group2)
  
  dp1=ecdf(s1)(cutoffs)
  dp2=ecdf(s2)(cutoffs)
  
  se1=sqrt(dp1*(1-dp1)/n1)
  se2=sqrt(dp2*(1-dp2)/n2)
  
  out=rbind(dp1,se1,dp2,se2)
  return(out)
}

##' group_fairness
##' 
##' Estimates group fairness metric according to a specification vector of the form
##' 
##' c(A1,B1,C1,A2,B2,C2)
##' 
##' encoding a probability
##' 
##' P(A1,B1,C1|A2,B2,C2)
##' 
##' where 
##' 
##'  A1/A2 are events coded by 1:'score>= cutoff'; 0: 'score<cutoff' and NA: 1/TRUE
##'  B1/B2 are events coded by 1:'target=TRUE'; 0: 'target=FALSE' and NA: 1/TRUE
##'  C1/C2 are events coded by 1:'group=g'; and NA: 1/TRUE
##' 
##' For example, specs=c(NA,1,NA,0,NA,1) would encode false omission rate:
##' 
##' P(target=TRUE|score<cutoff,group=g)
##' 
##' @param specs specification vector; see description
##' @param scores vector of risk scores
##' @param target vector of values of target (which risk score aims to predict)
##' @param group1 indices of group 1
##' @param group2 indices of group 2
##' @param cutoffs score cutoffs at which to estimate metric (default 100 evenly-spaced)
##' @export
##' @return matrix of dimension length(cutoffs)x4, with (i,2g-1)th entry the relevant fairness metric for group g at the ith cutoff value and (i,2g)th entry the approximate standard error of the (i,2g-1)th entry
##' @examples 
##' # See vignette
group_fairness=function(specs,scores, target, group1, group2, 
                        cutoffs=seq(min(scores,na.rm=TRUE),max(scores,na.rm=TRUE),length=100)) {
  
  ncut=length(cutoffs)
  ndat=length(scores)
  
  p_AA=rep(NA,ncut); p_BB=rep(NA,ncut)
  n_AA=rep(NA,ncut); n_BB=rep(NA,ncut)
  for (i in 1:ncut) {
    cond_AA=1:ndat
    cond_BB=1:ndat
    if (is.finite(specs[4])) {
      if (specs[4]==1) yhatc=which(scores>=cutoffs[i]) else yhatc=which(scores < cutoffs[i]) 
      cond_AA=intersect(cond_AA,yhatc)
      cond_BB=intersect(cond_BB,yhatc)
    }
    if (is.finite(specs[5])) {
      yc=which(target==specs[5])
      cond_AA=intersect(cond_AA,yc)
      cond_BB=intersect(cond_BB,yc)
    }
    if (is.finite(specs[6])) {
      cond_AA=intersect(cond_AA,group1)
      cond_BB=intersect(cond_BB,group2)
    }
    
    prob_AA=cond_AA
    prob_BB=cond_BB
    if (is.finite(specs[1])) {
      if (specs[1]==1) yhatc=which(scores>=cutoffs[i]) else yhatc=which(scores < cutoffs[i]) 
      prob_AA=intersect(prob_AA,yhatc)
      prob_BB=intersect(prob_BB,yhatc)
    }
    if (is.finite(specs[2])) {
      yc=which(target==specs[2])
      prob_AA=intersect(prob_AA,yc)
      prob_BB=intersect(prob_BB,yc)
    }
    if (is.finite(specs[3])) {
      prob_AA=intersect(prob_AA,group1)
      prob_BB=intersect(prob_BB,group2)
    }
    
    p_AA[i]=length(prob_AA)/length(cond_AA)
    p_BB[i]=length(prob_BB)/length(cond_BB)
    n_AA[i]=length(cond_AA)
    n_BB[i]=length(cond_BB)
  }
  
  se_AA=sqrt(p_AA*(1-p_AA)/n_AA)  
  se_BB=sqrt(p_BB*(1-p_BB)/n_BB)  
  
  out=rbind(p_AA,se_AA,p_BB,se_BB)
  
  return(out)
}







##' adjusted_for
##' 
##' Estimates false omission rate P(target=TRUE|score<=cutoff,group=g) 'adjusted' for some category. 
##' 
##' Namely, calculates 
##' 
##' sum ( P(target=TRUE|score<=cutoff,category=c,group=g)P(category=c|score<cutoff) )
##' 
##' where the sum is over categories c.
##' 
##' @param scores vector of risk scores
##' @param target vector of values of target (which risk score aims to predict)
##' @param group1 indices of group 1
##' @param group2 indices of group 2
##' @param category vector of categories
##' @param cutoffs score cutoffs at which to estimate metric (default 100 evenly-spaced)
##' @param nboot number of bootstrap samples for standard error
##' @export
##' @return matrix of dimension length(cutoffs)x4, with (i,2g-1)th entry the relevant fairness metric for group g at the ith cutoff value and (i,2g)th entry the approximate standard error of the (i,2g-1)th entry
##' @examples 
##' # See vignette
adjusted_for=function(scores, target, category, 
                      group1, group2, 
                      cutoffs=seq(min(scores,na.rm=TRUE),max(scores,na.rm=TRUE),length=100),
                      nboot=100) {
  
  ncut=length(cutoffs)
  ndat=length(scores)
  
  dsub=data.frame(score=scores,target=target,cat=category)
  dsub=dsub[order(as.character(dsub$cat)),]
  
  dt=table(dsub$cat)[unique(dsub$cat)]
  cdt=c(0,cumsum(dt))
  dsub$group=rep(NA,dim(dsub)[1])
  dsub$group[group1]="A"
  dsub$group[group2]="B"
  
  # Matrix of counts of A=a|Yhat<c and A=a|Yhat>= c for cutoffs c
  ucat=unique(dsub$cat)
  dt_forp=matrix(0,length(ucat),length(cutoffs))
  for (i in 1:length(cutoffs)) {
    dt_forp[,i]=table(dsub$cat[which(dsub$score<cutoffs[i])])[ucat]
  }
  dt_forp[which(is.na(dt_forp))]=0
  
  
  # Arrays of counts of A=a|Yhat<c and A=a|Yhat>= c for cutoffs c for bootstrap samples
  b_dt_forp=array(0,dim=c(length(ucat),length(cutoffs),nboot))
  b_n=array(0,dim=c(length(ucat),nboot))
  b_A=array(0,dim=c(length(ucat),nboot))
  b_B=array(0,dim=c(length(ucat),nboot))
  for (j in 1:nboot) {
    dsub_b=dsub[sample(dim(dsub)[1],replace=TRUE),]
    b_n[,j]=table(dsub_b$cat)[ucat]
    b_A[,j]=table(dsub_b$cat[which(dsub_b$group=="A")])[ucat]
    b_B[,j]=table(dsub_b$cat[which(dsub_b$group=="B")])[ucat]
    for (i in 1:length(cutoffs)) {
      b_dt_forp[,i,j]=table(dsub_b$cat[which(dsub_b$score<cutoffs[i])])[ucat]
    }
  }
  b_n[which(is.na(b_n))]=0
  b_A[which(is.na(b_A))]=0
  b_B[which(is.na(b_B))]=0
  b_dt_forp[which(is.na(b_dt_forp))]=0
  
  
  # FORP/FDRP for each category
  forp_A=matrix(NA,length(dt),length(cutoffs))
  forp_B=matrix(NA,length(dt),length(cutoffs))
  
  # Denominators for each category
  nforp_A=matrix(0,length(dt),length(cutoffs))
  nforp_B=matrix(0,length(dt),length(cutoffs))
  
  b_forp_A=array(NA,dim=c(length(dt),length(cutoffs),nboot))
  b_forp_B=array(NA,dim=c(length(dt),length(cutoffs),nboot))
  
  
  for (i in 1:length(dt)) {
    sub=dsub[(1+cdt[i]):(cdt[i+1]),]
    
    subA=sub[which(sub$group=="A"),]
    if (length(unique(subA$score))>2) {
      nA1=round(ecdf(subA$score)(cutoffs)*length(subA$score))
      nA2=round((1-ecdf(subA$score)(cutoffs))*length(subA$score))
      w1=which(nA1>0); w2=which(nA2>0)
      forp_A[i,w1]=cumsum(subA$target)[nA1[w1]]/nA1[w1]
    }
    for (j in 1:nboot) {
      subAc=subA[sort(sample(dim(subA)[1],b_A[i,j],replace=TRUE)),]
      if (length(unique(subAc$score))>2) {
        nA1=round(ecdf(subAc$score)(cutoffs)*length(subAc$score))
        nA2=round((1-ecdf(subAc$score)(cutoffs))*length(subAc$score))
        w1=which(nA1>0); w2=which(nA2>0)
        b_forp_A[i,w1,j]=cumsum(subAc$target)[nA1[w1]]/nA1[w1]
      }
    }
    
    subB=sub[which(sub$group=="B"),]
    if (length(unique(subB$score))>2) {
      nB1=round(ecdf(subB$score)(cutoffs)*length(subB$score))
      nB2=round((1-ecdf(subB$score)(cutoffs))*length(subB$score))
      w1=which(nB1>0); w2=which(nB2>0)
      forp_B[i,w1]=cumsum(subB$target)[nB1[w1]]/nB1[w1]
    }
    for (j in 1:nboot) {
      subBc=subB[sort(sample(dim(subB)[1],b_B[i,j],replace=TRUE)),]
      if (length(unique(subBc$score))>2) {
        nB1=round(ecdf(subBc$score)(cutoffs)*length(subBc$score))
        nB2=round((1-ecdf(subBc$score)(cutoffs))*length(subBc$score))
        w1=which(nB1>0); w2=which(nB2>0)
        b_forp_B[i,w1,j]=cumsum(subBc$target)[nB1[w1]]/nB1[w1]
      }
    }
  }
  
  FORP_A1=rep(0,length(cutoffs))
  FORP_B1=rep(0,length(cutoffs))
  b_FORP_A1=matrix(0,length(cutoffs),nboot)
  b_FORP_B1=matrix(0,length(cutoffs),nboot)
  
  for (i in 1:length(cutoffs)) {
    FORP_A1[i]=sum(forp_A[,i]*dt_forp[,i],na.rm=TRUE) 
    FORP_B1[i]=sum(forp_B[,i]*dt_forp[,i],na.rm=TRUE) 
    for (j in 1:nboot) {
      b_FORP_A1[i,j]=sum(b_forp_A[,i,j]*b_dt_forp[,i,j],na.rm=TRUE) 
      b_FORP_B1[i,j]=sum(b_forp_B[,i,j]*b_dt_forp[,i,j],na.rm=TRUE) 
    }
  }
  
  # Normalise
  FORP_A=rep(0,length(cutoffs))
  FORP_B=rep(0,length(cutoffs))
  b_FORP_A=matrix(0,length(cutoffs),nboot)
  b_FORP_B=matrix(0,length(cutoffs),nboot)
  for (i in 1:length(cutoffs)) {
    FORP_A[i]=FORP_A1[i]/sum(dt_forp[which(is.finite(forp_A[,i])),i])
    FORP_B[i]=FORP_B1[i]/sum(dt_forp[which(is.finite(forp_B[,i])),i])
    for (j in 1:nboot) {
      b_FORP_A[i,j]=b_FORP_A1[i,j]/sum(b_dt_forp[which(is.finite(b_forp_A[,i,j])),i,j])
      b_FORP_B[i,j]=b_FORP_B1[i,j]/sum(b_dt_forp[which(is.finite(b_forp_B[,i,j])),i,j])
    }
  }
  
  
  # Standard errors and confidence intervals are approximate, and take weights as fixed
  qx=function(x,alpha) if (any(is.finite(x))) quantile(x,alpha,na.rm=TRUE) else NA
  se_FORP_A=apply(b_FORP_A,1,function(x) qx(x,0.841)-qx(x,0.159))/2 # one SE
  se_FORP_B=apply(b_FORP_B,1,function(x) qx(x,0.841)-qx(x,0.159))/2 
  
  out=rbind(FORP_A,as.vector(se_FORP_A),FORP_B,as.vector(se_FORP_B))
  
  return(out)
  
}






##' adjusted_fdr
##' 
##' Estimates false discovery rate P(target=FALSE|score>cutoff,group=g) 'adjusted' for some category. 
##' 
##' Namely, calculates 
##' 
##' sum ( P(target=FALSE|score>cutoff,category=c,group=g)P(category=c|score<cutoff) )
##' 
##' where the sum is over categories c.
##' 
##' @param scores vector of risk scores
##' @param target vector of values of target (which risk score aims to predict)
##' @param group1 indices of group 1
##' @param group2 indices of group 2
##' @param category vector of categories
##' @param cutoffs score cutoffs at which to estimate metric (default 100 evenly-spaced)
##' @param nboot number of bootstrap samples for standard error
##' @export
##' @return matrix of dimension length(cutoffs)x4, with (i,2g-1)th entry the relevant fairness metric for group g at the ith cutoff value and (i,2g)th entry the approximate standard error of the (i,2g-1)th entry
##' @examples 
##' # See vignette
adjusted_fdr=function(scores, target, category, 
                      group1, group2, 
                      cutoffs=seq(min(scores,na.rm=TRUE),max(scores,na.rm=TRUE),length=100),
                      nboot=100) {
  
  ncut=length(cutoffs)
  ndat=length(scores)
  
  dsub=data.frame(score=scores,target=target,cat=category)
  dsub=dsub[order(as.character(dsub$cat)),]
  
  dt=table(dsub$cat)[unique(dsub$cat)]
  cdt=c(0,cumsum(dt))
  dsub$group=rep(NA,dim(dsub)[1])
  dsub$group[group1]="A"
  dsub$group[group2]="B"
  
  # Matrix of counts of A=a|Yhat<c and A=a|Yhat>= c for cutoffs c
  ucat=unique(dsub$cat)
  dt_fdrp=matrix(0,length(ucat),length(cutoffs))
  for (i in 1:length(cutoffs)) {
    dt_fdrp[,i]=table(dsub$cat[which(dsub$score>=cutoffs[i])])[ucat]
  }
  dt_fdrp[which(is.na(dt_fdrp))]=0
  
  
  # Arrays of counts of A=a|Yhat<c and A=a|Yhat>= c for cutoffs c for bootstrap samples
  b_dt_forp=array(0,dim=c(length(ucat),length(cutoffs),nboot))
  b_dt_fdrp=array(0,dim=c(length(ucat),length(cutoffs),nboot))
  b_n=array(0,dim=c(length(ucat),nboot))
  b_A=array(0,dim=c(length(ucat),nboot))
  b_B=array(0,dim=c(length(ucat),nboot))
  for (j in 1:nboot) {
    dsub_b=dsub[sample(dim(dsub)[1],replace=TRUE),]
    b_n[,j]=table(dsub_b$cat)[ucat]
    b_A[,j]=table(dsub_b$cat[which(dsub_b$group=="A")])[ucat]
    b_B[,j]=table(dsub_b$cat[which(dsub_b$group=="B")])[ucat]
    for (i in 1:length(cutoffs)) {
      b_dt_forp[,i,j]=table(dsub_b$cat[which(dsub_b$score<cutoffs[i])])[ucat]
      b_dt_fdrp[,i,j]=b_n[,j]-b_dt_forp[,i,j]
    }
  }
  b_n[which(is.na(b_n))]=0
  b_A[which(is.na(b_A))]=0
  b_B[which(is.na(b_B))]=0
  b_dt_fdrp[which(is.na(b_dt_fdrp))]=0
  
  
  # FORP/FDRP for each category
  fdrp_A=matrix(NA,length(dt),length(cutoffs))
  fdrp_B=matrix(NA,length(dt),length(cutoffs))
  
  # Denominators for each category
  nfdrp_A=matrix(0,length(dt),length(cutoffs))
  nfdrp_B=matrix(0,length(dt),length(cutoffs))
  
  b_fdrp_A=array(NA,dim=c(length(dt),length(cutoffs),nboot))
  b_fdrp_B=array(NA,dim=c(length(dt),length(cutoffs),nboot))
  
  
  for (i in 1:length(dt)) {
    sub=dsub[(1+cdt[i]):(cdt[i+1]),]
    
    subA=sub[which(sub$group=="A"),]
    if (length(unique(subA$score))>2) {
      nA1=round(ecdf(subA$score)(cutoffs)*length(subA$score))
      nA2=round((1-ecdf(subA$score)(cutoffs))*length(subA$score))
      w1=which(nA1>0); w2=which(nA2>0)
      fdrp_A[i,w2]=cumsum(rev(1-subA$target))[nA2[w2]]/nA2[w2]
    }
    for (j in 1:nboot) {
      subAc=subA[sort(sample(dim(subA)[1],b_A[i,j],replace=TRUE)),]
      if (length(unique(subAc$score))>2) {
        nA1=round(ecdf(subAc$score)(cutoffs)*length(subAc$score))
        nA2=round((1-ecdf(subAc$score)(cutoffs))*length(subAc$score))
        w1=which(nA1>0); w2=which(nA2>0)
        b_fdrp_A[i,w2,j]=cumsum(rev(1-subAc$target))[nA2[w2]]/nA2[w2]
      }
    }
    
    subB=sub[which(sub$group=="B"),]
    if (length(unique(subB$score))>2) {
      nB1=round(ecdf(subB$score)(cutoffs)*length(subB$score))
      nB2=round((1-ecdf(subB$score)(cutoffs))*length(subB$score))
      w1=which(nB1>0); w2=which(nB2>0)
      fdrp_B[i,w2]=cumsum(rev(1-subB$target))[nB2[w2]]/nB2[w2]
    }
    for (j in 1:nboot) {
      subBc=subB[sort(sample(dim(subB)[1],b_B[i,j],replace=TRUE)),]
      if (length(unique(subBc$score))>2) {
        nB1=round(ecdf(subBc$score)(cutoffs)*length(subBc$score))
        nB2=round((1-ecdf(subBc$score)(cutoffs))*length(subBc$score))
        w1=which(nB1>0); w2=which(nB2>0)
        b_fdrp_B[i,w2,j]=cumsum(rev(1-subBc$target))[nB2[w2]]/nB2[w2]
      }
    }
  }
  
  FDRP_A1=rep(0,length(cutoffs))
  FDRP_B1=rep(0,length(cutoffs))
  b_FDRP_A1=matrix(0,length(cutoffs),nboot)
  b_FDRP_B1=matrix(0,length(cutoffs),nboot)
  
  for (i in 1:length(cutoffs)) {
    FDRP_A1[i]=sum(fdrp_A[,i]*dt_fdrp[,i],na.rm=TRUE) 
    FDRP_B1[i]=sum(fdrp_B[,i]*dt_fdrp[,i],na.rm=TRUE) 
    for (j in 1:nboot) {
      b_FDRP_A1[i,j]=sum(b_fdrp_A[,i,j]*b_dt_fdrp[,i,j],na.rm=TRUE) 
      b_FDRP_B1[i,j]=sum(b_fdrp_B[,i,j]*b_dt_fdrp[,i,j],na.rm=TRUE) 
    }
  }
  
  # Normalise
  FDRP_A=rep(0,length(cutoffs))
  FDRP_B=rep(0,length(cutoffs))
  b_FDRP_A=matrix(0,length(cutoffs),nboot)
  b_FDRP_B=matrix(0,length(cutoffs),nboot)
  for (i in 1:length(cutoffs)) {
    FDRP_A[i]=FDRP_A1[i]/sum(dt_fdrp[which(is.finite(fdrp_A[,i])),i])
    FDRP_B[i]=FDRP_B1[i]/sum(dt_fdrp[which(is.finite(fdrp_B[,i])),i])
    for (j in 1:nboot) {
      b_FDRP_A[i,j]=b_FDRP_A1[i,j]/sum(b_dt_fdrp[which(is.finite(b_fdrp_A[,i,j])),i,j])
      b_FDRP_B[i,j]=b_FDRP_B1[i,j]/sum(b_dt_fdrp[which(is.finite(b_fdrp_B[,i,j])),i,j])
    }
  }
  
  
  # Standard errors and confidence intervals are approximate, and take weights as fixed
  qx=function(x,alpha) if (any(is.finite(x))) quantile(x,alpha,na.rm=TRUE) else NA
  se_FDRP_A=apply(b_FDRP_A,1,function(x) qx(x,0.841)-qx(x,0.159))/2
  se_FDRP_B=apply(b_FDRP_B,1,function(x) qx(x,0.841)-qx(x,0.159))/2
  
  out=rbind(FDRP_A,as.vector(se_FDRP_A),FDRP_B,as.vector(se_FDRP_B))
  
  return(out)
  
}





##' getroc() 
##' Comprehensive plotting function for receiver-operator characteristic curve. Also calculates AUROC and standard error. 
##' 
##' Rather than returning points corresponding to every cutoff, only returns a representative sample of equally-spaced points along the curve.
##'
##' SE of AUROC with no CV structure is from Hanley and McNeil 1982. SE of AUROC with CV folds is from LeDell et al 2012
##'
##' Does not plot anything. Object can be plotted in a default way.
##'
##' @param y class labels, 0/1 or logical
##' @param ypred predictions Pr(Y=1), numeric vector
##' @param cv cross-validation fold assignments, if relevant. Changes estimate of standard error.
##' @param res resolution. Returns this many equally-spaced points along the curve. Set res to null to return all points.
##' @export 
##' @return list containing: spec, specificity for res points in every cv fold; sens, sensitivity for res points in every cv fold; auc, areas under the curve for each fold and average (note length is 1 greater than number of CV folds); se, standard error for AUC in each fold and standard error for average auc (note length is 1 greater than number of CV folds)
##' @examples 
##' # See vignette
getroc=function(y,ypred,cv=NULL,res=100) {
  if (is.null(cv)) cv=rep(1,length(y))
  if (!(length(y)==length(ypred))) stop("Parameters y and ypred should have the same length")
  
  sens=c(); spec=c(); auc=c(); se=c(); cutoffs=c(); ncut=c()
  for (i in 1:max(cv)) {
    y0=y[which(cv==i)]; 
    ypred0=ypred[which(cv==i)]
    
    yt=sum(y0); yl=length(y0)
    opred=order(ypred0)
    #ipred=order(opred) # can use ipred to reorder in the order of original ypred
    
    sy=y0[opred]; sp=ypred0[opred]
    
    # Cutoffs and number of samples
    cutoffs0=sp
    ncut0=1:length(sp)
    
    # Disclosure control
    csy=cumsum(sy); csy[which(csy<5)]=0
    csy1=cumsum(1-sy); csy1[which(csy1<5)]=0
    ncut0[which(ncut0<5)]=0
    
    sens0=1- (csy/yt)
    spec0= csy1/(yl-yt)
    
    auc0=integral(sens0,spec0)
    se0=aucse(as.numeric(yt),as.numeric(yl-yt),auc0)
    
    if (!is.null(res)) {
      ds=cumsum(sqrt((spec0[1:(yl-1)]-spec0[2:yl])^2 + (sens0[1:(yl-1)]-sens0[2:yl])^2))
      ds=ds/ds[yl-1]
      lsp=(1:(yl-1))/yl
      sub=round(yl*approx(ds,lsp,n=res)$y)
      sens0=sens0[sub]
      spec0=spec0[sub]
      cutoffs0=cutoffs0[sub]
      ncut0=ncut0[sub]
    }
    
    auc=c(auc,auc0)
    se=c(se,se0)
    spec=rbind(spec,spec0)
    sens=rbind(sens,sens0)
    cutoffs=rbind(cutoffs,cutoffs0)
    ncut=rbind(ncut,ncut0)
  }
  
  if (length(auc)>1) {
    auc=c(auc,mean(auc))
    se=c(se,ci.cvAUC(ypred,y,folds=cv)$se)
  }
  
  out=list(sens=sens,spec=spec,cutoffs=cutoffs,ncut=ncut,auc=auc,se=se)
  class(out)="sparraROC"
  return(out)
}

# Internal function to compute SE of AUC
aucse=function(n1,n2,auc) {
  q1=auc/(2-auc); q2=2*(auc^2)/(1+auc)
  num=auc*(1-auc) + (n1-1)*(q1- (auc^2)) + (n2-1)*(q2-(auc^2))
  return(sqrt(num/(n1*n2)))
}


##' getprc() 
##' 
##' Comprehensive plotting function for precision-recall curve. Also calculates AUPRC and standard error. 
##' 
##' Rather than returning points corresponding to every cutoff, only returns a representative sample of equally-spaced points along the curve.
##'
##' Does not plot anything. Object can be plotted in a default way.
##'
##' @param y class labels, 0/1 or logical
##' @param ypred predictions Pr(Y=1), numeric vector
##' @param cv cross-validation fold assignments, if relevant. Changes estimate of standard error.
##' @param res resolution. Returns this many equally-spaced points along the curve. Set res to null to return all points.
##' @export 
##' @return list containing: ppv, ppv for res points in every cv fold; sens, sensitivity for res points in every cv fold; auc, areas under the curve for each fold and average (note length is 1 greater than number of CV folds); se, standard error for AUC in each fold and standard error for average auc (note length is 1 greater than number of CV folds)
##' @examples 
##' # See vignette
getprc=function(y,ypred,cv=NULL,res=100) {
  if (is.null(cv)) cv=rep(1,length(y))
  if (!(length(y)==length(ypred))) stop("Parameters y and ypred should have the same length")
  
  sens=c(); ppv=c(); auc=c(); se=c(); cutoffs=c(); ncut=c()
  for (i in 1:max(cv)) {
    y0=y[which(cv==i)]; 
    ypred0=ypred[which(cv==i)]
    
    yt=sum(y0); yl=length(y0)
    opred=order(ypred0)
    #ipred=order(opred) # can use ipred to reorder in the order of original ypred
    sy=y0[opred]; sp=ypred0[opred]
    
    cutoffs0=sp
    ncut0=1:length(sp)
    
    # Disclosure control
    crsy=cumsum(rev(sy)); crsy[which(crsy<5)] = 0
    csy=cumsum(sy); csy[which(csy<5)] = 0
    ncut0[which(ncut0<5)] = 0
    
    ppv0=rev(crsy/(1:length(sy)))
    sens0=1- (csy/yt)
    
    auc0=integral(sens0,ppv0)
    
    
    se0=sqrt(auc0*(1-auc0)/sum(y0))
    
    if (!is.null(res)) {
      ds=cumsum(sqrt((ppv0[1:(yl-1)]-ppv0[2:yl])^2 + (sens0[1:(yl-1)]-sens0[2:yl])^2))
      ds=ds/ds[yl-1]
      lsp=(1:(yl-1))/yl
      sub=suppressWarnings(round(yl*approx(ds,lsp,n=res)$y))
      sens0=sens0[sub]
      ppv0=ppv0[sub]
      cutoffs0=cutoffs0[sub]
      ncut0=ncut0[sub]
    }
    
    auc=c(auc,auc0)
    se=c(se,se0)
    ppv=rbind(ppv,ppv0)
    sens=rbind(sens,sens0)
    cutoffs=rbind(cutoffs,cutoffs0)
    ncut=rbind(ncut,ncut0)
  }
  
  if (length(auc)>1) {
    auc=c(auc,mean(auc))
    se=c(se,mean(se)/sqrt(3))
  }
  
  
  out=list(sens=sens,ppv=ppv,cutoffs=cutoffs,ncut=ncut,auc=auc,se=se)
  class(out)="sparraPRC"
  return(out)
}

##' getcal()
##' 
##' Produces a set of points for a calibration plot.
##' 
##' Uses either a binning method or a kernel method to determine height of points. 
##' 
##' In both methods, considers n equally spaced subintervals of (0,1)
##'
##' @param y class labels, 0/1 or logical
##' @param ypred predictions Pr(Y=1), numeric vector
##' @param n number of subintervals/points
##' @param kernel set to TRUE to use kernel method
##' @param kernel_sd kernel width for kernel method; see above
##' @param alpha return a pointwise confidence envolope for conservative 1-alpha confidence interval
##' @param c0 for computing maximum bias; assume true covariance function is of the form a0+ a1x + a2x^2, with |a0|<c0, |a2|<c2 (c1 does not matter)
##' @param c2 for computing maximum bias; assume true covariance function is of the form a0+ a1x + a2x^2, with |a0|<c0, |a2|<c2 (c1 does not matter)
##' @export
##' @return a list with components x (expected calibration), y (observed calibration), n (number of samples in bins, if relevant), lower/upper (confidence interval on y)
##' @examples 
##' # See vignette
getcal=function(y,ypred,n=10,kernel=FALSE,kernel_sd=0.05,alpha=0.05,c0=0,c2=0.1) {
  if (!kernel) {
    ycal=rep(0,n); xcal=ycal; ncal=ycal
    xup=rep(0,n); xdown=rep(0,n)
    for (i in 1:n) {
      sub=which((ypred> (i-1)/n) & (ypred< i/n))
      ycal[i]=mean(y[sub])
      xcal[i]=mean(ypred[sub])
      ncal[i]=length(sub)
      xse=sqrt(ycal[i]*(1-ycal[i])/length(sub))
      xup[i]=ycal[i] -qnorm(alpha/2)*xse
      xdown[i]=ycal[i] + qnorm(alpha/2)*xse
    }
    obj=list(x=xcal,y=ycal,n=ncal,upper=xup,lower=xdown)
  } else {  # use kernel method with given sd
    ypredsub=seq(0,1,length=n)
    kern=function(x,y) dnorm(x-y,sd=kernel_sd)
    wt=outer(ypredsub,ypred,kern)
    x1=(wt %*% ypred)
    y1=(wt %*% y)
    csub=ypredsub*y1/x1
    csub=pmax(pmin(csub,1),0)
    nsub=rep(0,length(csub))
    for (i in 1:(length(nsub)-1)) {
      sub=which(ypred>=csub[i] & ypred < csub[i+1])
    }
    
    # Confidence intervals
    wts=ypredsub*wt/as.vector(x1) # now csub= wts %*% y
    yvar=(wts^2 %*% as.vector(ypred*(1-ypred)))
    
    # Max bias
    bias=rowSums(outer(ypredsub,ypred,function(x,y) kern(x,y)*(c0*(y-x) + c2*(x^2 * y - y^2 * x))))/x1
    
    upper=csub - qnorm(alpha/2)*sqrt(yvar) + bias
    lower=csub + qnorm(alpha/2)*sqrt(yvar) - bias
    
    return(list(x=ypredsub,y=csub,n=nsub,lower=lower,upper=upper))
  }
  class(obj)="sparraCAL"
  return(obj)
}

##**************************************************#
## Plotting                                      ####
##**************************************************#


##' Plot function for class sparraCAL
##' @param x output from getcal()
##' @param cols colour to draw lines
##' @param add set to FALSE to add to existing plot
##' @param add_xy_line set to TRUE to draw an X-Y reference line.
##' @param ... passed to lines()
##' @return No return value, draws a figure
##' @export
##' @examples 
##' # See vignette
plot.sparraCAL=function(x,cols=rep(phs_colours("phs-blue"),dim(x$sens)[1]),add=FALSE,add_xy_line=TRUE,...) {
  if (!add) plot(0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="Predicted",ylab="Observed",...)
  ncv=dim(x$spec)[1]
  lines(x$x,x$y,...)
  if (add_xy_line) abline(0,1,col=phs_colours("phs-magenta"),lty=2)
}



##' Plot function for class sparraROC
##' @param x output from getroc()
##' @param addauc set to TRUE to add text to the plot showing the (mean) AUC and SE.
##' @param cols colour to draw lines
##' @param ... passed to plot()
##' @return No return value, draws a figure
##' @export
##' @examples 
##' # See vignette
plot.sparraROC=function(x,addauc=FALSE,cols=rep(phs_colours("phs-blue"),dim(x$sens)[1]),...) {
  plot(0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="1-Spec.",ylab="Sens.",...)
  ncv=dim(x$spec)[1]
  for (i in 1:ncv) lines(1-x$spec[i,],x$sens[i,],col=cols[i])
  abline(0,1,col=phs_colours("phs-magenta"),lty=2)
  auc=x$auc[length(x$auc)]
  se=x$se[length(x$se)]
  txx=paste0(signif(auc,digits=2),"+/-",signif(se,digits=2))
  if (addauc) text(0.6,0.4,txx)
}


##' Plot function for class above
##' @param x output from getprc()
##' @param addauc set to TRUE to add text to the plot showing the (mean) AUC and SE.
##' @param cols colour to draw lines
##' @param ... passed to plot()
##' @return No return value, draws a figure
##' @export
##' @examples 
##' # See vignette
plot.sparraPRC=function(x,addauc=FALSE,cols=rep(phs_colours("phs-blue"),dim(x$sens)[1]),...) {
  plot(0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="Recall",ylab="Precision",...)
  ncv=dim(x$sens)[1]
  for (i in 1:ncv) lines(x$sens[i,],x$ppv[i,],col=cols[i])
  auc=mean(x$auc)
  se=mean(x$se)/sqrt(3)
  txx=paste0(signif(auc,digits=2),"+/-",signif(se,digits=2))
  if (addauc) text(0.6,0.4,txx)
}




##' roc_2panel
##' Draws a ROC curve (with legend) with a second panel underneath showing sensitivity difference.
##' 
##' @param rocs list of sparraROC objects (if one object, plots folds separately)
##' @param labels labels to use in legend
##' @param col line colours
##' @param lty line type, defaults to 1
##' @param xy_lty line type for x-y line, defaults to 2 (dashed)
##' @param xy_col line colour for x-y line, defaults to red
##' @param highlight if non-null, add a point at this cutoff
##' @param mar_scale scale bottom and left margins by this amount. Also scales legend.
##' @param yrange_lower y range for lower plot. If NULL, generates automatically
##' @param ... other parameters passed to legend()
##' @return No return value, draws a figure
##' @export
##' @examples 
##' # See vignette
roc_2panel=function(rocs,labels,col=1:length(rocs),
                    lty=rep(1,length(col)),xy_lty=2,xy_col=phs_colours("phs-magenta"),
                    highlight=NULL,mar_scale=1,yrange_lower=NULL,...) {
  
  if ("sens" %in% names(rocs)) {
    nfold=length(rocs$auc)-1; np=length(rocs$sens)/nfold
    r0=list()
    for (i in 1:nfold) 
      r0[[i]]=list(sens=rocs$sens[i,],
                   spec=rocs$spec[i,])
    rocs=r0
  }
  
  # Set up plot parameters
  oldpar = par(no.readonly = TRUE)
  on.exit(par(oldpar)) 
  par(mar=c(1,mar_scale*4,0.1,0.1))
  layout(mat=rbind(matrix(1,4,4),matrix(2,2,4)))
  
  # Initialise
  plot(0,xlim=c(0,1),ylim=c(0,1),xaxt="n",ylab="Sensitivity",type="n")
  abline(0,1,col=xy_col,lty=xy_lty)
  
  # x-values to draw sensitivity at
  xspec=seq(0,1,length=100)[2:99]; xr=c()
  
  # Draw ROC curves on top panel
  for (i in 1:length(rocs)) {
    xy=rocs[[i]]
    lines(1-xy$spec,xy$sens,col=col[i],lty=lty[i])
    xsens=suppressWarnings(approx(1-xy$spec,xy$sens,xspec)$y)
    xr=rbind(xr,xsens)
    if (!is.null(highlight)) {
      w0=which.min(abs(xy$cutoffs-highlight))
      points(1-xy$spec[w0],xy$sens[w0],pch=16,col=col[i])
    }
  }
  
  # Add legend
  par(ps=par()$ps/mar_scale)
  legend("bottomright",legend=labels,col=col,lty=lty,cex=mar_scale,...)
  par(ps=par()$ps*mar_scale)
  
  # Bottom panel setup
  par(mar=c(mar_scale*4,mar_scale*4,0.1,0.1))
  if (is.null(yrange_lower)) {
    yrr=range(t(xr)-xr[1,],na.rm=TRUE); 
    if (!is.finite(sum(yrr))) yrr=c(-1,1)
  } else yrr=yrange_lower
  
  plot(0,xlim=c(0,1),ylim=yrr,type="n",
       xlab="1-Specificity",ylab=expression(paste(Delta,"(sens.)")),yaxt="n")
  axis(2,at=pretty(yrr,n=2))
  
  # Draw lines on bottom panel
  for (i in 1:length(rocs)) {
    lines(xspec,xr[i,]-xr[1,],col=col[i],lty=lty[i])
  }
  #abline(h=0,col=xy_col,lty=xy_lty)
  
}  














##' prc_2panel
##' Draws a PRC curve (with legend) with a second panel underneath showing precision difference.
##' 
##' @param prcs list of sparraPRC objects. If of length 1, splits into folds
##' @param labels labels to use in legend
##' @param col line colours
##' @param lty line type, defaults to 1
##' @param highlight if non-null, draw a point at a particular cutoff
##' @param mar_scale scale bottom and left margins by this amount. Also scales legend.
##' @param yrange_lower y range for lower plot. If NULL, generates automatically
##' @param ... other parameters passed to legend()
##' @return No return value, draws a figure
##' @export
##' @examples 
##' # See vignette
prc_2panel=function(prcs,labels,col=1:length(prcs),
                    lty=rep(1,length(col)),highlight=NULL,mar_scale=1,yrange_lower=NULL,...) {
  
  if ("ppv" %in% names(prcs)) {
    nfold=length(prcs$auc)-1; np=length(prcs$sens)/nfold
    r0=list()
    for (i in 1:nfold) 
      r0[[i]]=list(sens=prcs$sens[i,],
                   ppv=prcs$ppv[i,])
    prcs=r0
  }
  
  
  
  # Set up plot 
  oldpar = par(no.readonly = TRUE)
  on.exit(par(oldpar)) 
  par(mar=c(1,mar_scale*4,0.1,0.1))
  layout(mat=rbind(matrix(1,4,4),matrix(2,2,4)))
  
  # Initialise plot
  plot(0,xlim=c(0,1),ylim=c(0,1),xaxt="n",ylab="Precision",type="n")
  xr=c(); 
  
  # X-values at which difference in precision will be plotted
  xsens=seq(0,1,length=100)[2:99]
  
  # Draw curves on top panel
  for (i in 1:length(prcs)) {
    px=prcs[[i]]
    lines(px$sens,px$ppv,col=col[i],lty=lty[i])
    xppv=suppressWarnings(approx(px$sens,px$ppv,xsens)$y)
    xr=rbind(xr,xppv)
    if (!is.null(highlight)) {
      w0=which.min(abs(px$cutoffs-highlight))
      points(px$sens[w0],px$ppv[w0],pch=16,col=col[i])
    }
    
  }
  
  # Add legend
  par(ps=par()$ps/mar_scale)
  legend("topright",legend=labels,col=col,lty=lty,cex=mar_scale,...)
  par(ps=par()$ps*mar_scale)
  
  # Initialise bottom panel
  par(mar=c(mar_scale*4,mar_scale*4,0.1,0.1))
  if (is.null(yrange_lower)) {
    yrr=range(t(xr)-xr[1,],na.rm=TRUE); if (!is.finite(sum(yrr))) yrr=c(-1,1)
  } else yrr=yrange_lower
  
  plot(0,xlim=c(0,1),ylim=yrr,type="n",
       xlab="Recall",ylab=expression(paste(Delta,"(Prec.)")),
       yaxt="n")
  axis(2,at=pretty(range(t(xr)-xr[1,]),n=2))
  
  # Draw lines on bottom panel
  for (i in 1:length(prcs)) lines(xsens,xr[i,]-xr[1,],col=col[i],lty=lty[i])
}



##' cal_2panel
##' Draws calibration curves (with legend) with a second panel underneath showing predicted differences.
##' 
##' @param cals list of calibration objects, output from getcal(). 
##' @param labels labels to use in legend
##' @param col line colours
##' @param lty line type, defaults to 1
##' @param xy_lty line type for x-y line, defaults to 2 (dashed)
##' @param xy_col line colour for x-y line, defaults to phs-magenta
##' @param ci_col colours to draw confidence intervals on lower panel; NA to not draw. 
##' @param highlight if non-null, highlight a particular value
##' @param mar_scale scale bottom and left margins by this amount. Also scales legend.
##' @param yrange_lower y range for lower plot. If NULL, generates automatically
##' @param ... other parameters passed to legend()
##' @return No return value, draws a figure
##' @export
##' @examples 
##' # See vignette
cal_2panel=function(cals,labels,col=1:length(cals),
                    lty=rep(1,length(col)),xy_lty=2,xy_col=phs_colours("phs-magenta"),
                    ci_col=rep(NA,length(col)),highlight=NULL,mar_scale=1,yrange_lower=NULL,...) {
  
  
  # Setup plot
  oldpar = par(no.readonly = TRUE)
  on.exit(par(oldpar)) 
  par(mar=c(1,mar_scale*4,0.1,0.1))
  layout(mat=rbind(matrix(1,4,4),matrix(2,2,4)))
  
  # Initialise plot
  plot(0,xlim=c(0,1),ylim=c(0,1),xaxt="n",ylab="Observed",type="n")
  
  # Draw lines on top panel
  xr=c()
  for (i in 1:length(cals)) {
    cx=cals[[i]]
    lines(cx$x,cx$y,col=col[i],lty=lty[i])
    xr=rbind(xr,cx$y-cx$x)
    if (!is.null(highlight)) {
      w0=which.min(abs(cx$x-highlight))
      points(cx$x[w0],cx$y[w0],pch=16,col=col[i])
    }
  }
  abline(0,1,col=xy_col,lty=xy_lty,lwd=1)
  
  # Draw legend
  par(ps=par()$ps/mar_scale)
  legend("topleft",legend=labels,col=col,lty=lty,cex=mar_scale,...)
  par(ps=par()$ps*mar_scale)
  
  # Initialise bottom panel
  par(mar=c(mar_scale*4,mar_scale*4,0.1,0.1))
  xpred=seq(0,1,length=100)[2:99]
  
  if (is.null(yrange_lower)) {
    yrr=range(xr,na.rm=TRUE); if (!is.finite(sum(yrr))) yrr=c(-1,1)
  } else yrr=yrange_lower
  
  plot(0,xlim=c(0,1),ylim=yrr,type="n",
       xlab="Predicted",ylab=expression(paste(Delta,"(cal.)")),
       yaxt="n")
  axis(2,at=pretty(range(xr),n=3))
  
  # Draw lines on bottom panel
  for (i in 1:length(cals)) {
    cx=cals[[i]]
    if (!is.na(ci_col[i])) {
      xx=c(cx$x,rev(cx$x))
      yy=c(cx$lower,rev(cx$upper))-c(cx$x,rev(cx$x))
      wn=which(is.finite(xx+yy)); xx=xx[wn]; yy=yy[wn];
      polygon(xx,yy,
              col=ci_col[i],border=NA)
    }
    lines(cx$x,cx$y-cx$x,col=col[i],lty=lty[i])
  }
  abline(h=0,col=xy_col,lty=xy_lty)
  
}  







##' groupmetric_2panel
##' Draws plots of a group fairness metric with a second panel underneath
##' 
##' @param objs list of fairness objects. Each should contain sub-objects 'x', 'y' and 'ci', which specify x and y values and half-widths of confidence intervals around y.
##' @param labels labels to use in legend
##' @param col line colours
##' @param lty line type, defaults to 1
##' @param yrange limit of y axis; defaults to 0,1
##' @param ci_col confidence envelope colours. These will be transparent.
##' @param highlight if non-null, draw a point at a particular cutoff
##' @param logscale if TRUE, draw with log-scale.
##' @param mar_scale scale bottom and left margins by this amount. Also scales legend.
##' @param lpos legend position
##' @param yrange_lower y range for lower plot. If NULL, generates automatically
##' @param ... other parameters passed to legend()
##' @return No return value, draws a figure
##' @export
##' @examples 
##' # See vignette
groupmetric_2panel=function(objs,labels,col=1:length(objs),lty=rep(1,length(col)),yrange=NULL,
                            ci_col=1:length(objs),highlight=NULL,logscale=FALSE,
                            mar_scale=1,lpos=NULL,yrange_lower=NULL,...) {
  
  # Set up plot  
  oldpar = par(no.readonly = TRUE)
  on.exit(par(oldpar)) 
  par(mar=c(1,mar_scale*4,0.1,0.1))
  layout(mat=rbind(matrix(1,4,4),matrix(2,2,4)))
  
  if (logscale) {
    xvals=log(objs[[1]]$x); xvals[which(!is.finite(xvals))]=NA
    xrange=range(log(c(0.01,0.99)),na.rm=TRUE) 
  }
  else {
    xvals=objs[[1]]$x
    xrange=c(0,1)
  }
  
  if (is.null(yrange)) {
    y_all=c(); for (i in 1:length(objs)) y_all=c(y_all,objs[[i]]$y)
    yrange=range(y_all,na.rm=TRUE)
    if (!is.finite(sum(yrange))) yrange=c(0,1)
  }
  
  if (is.null(lpos)) lpos="topright"
  
  plot(0,xlim=xrange,ylim=yrange,type="n",xlab="",xaxt="n",ylab="Prob.",las=1)
  
  # Lines on top panel
  # X-values at which difference in precision will be plotted
  xsens=seq(0,1,length=100)[2:99]
  
  # Draw curves and confidence intervals on top panel
  yr=c()
  for (i in 1:length(objs)) {
    lines(xvals,objs[[i]]$y,col=col[i],lty=lty[i])
    
    # CI
    yci=c(objs[[i]]$y+objs[[i]]$ci,rev(objs[[i]]$y-objs[[i]]$ci)); wxx=which(is.finite(yci))
    cic=col2rgb(ci_col[i]); cicx=rgb(cic[1]/256,cic[2]/256,cic[3]/256,alpha=0.5)
    polygon(c(xvals,rev(xvals))[wxx],yci[wxx],col=cicx,border=NA)
    yr=rbind(yr,objs[[i]]$y-objs[[1]]$y)
  }
  
  # Legend
  par(ps=par()$ps/mar_scale)
  legend(lpos,labels,col=col,lty=lty,bg="white",cex=mar_scale,...)
  par(ps=par()$ps*mar_scale)
  
  # Add highlight line
  if (!is.null(highlight)) {
    if (logscale) abline(v=log(highlight)) else abline(v=highlight)
  }
  
  
  
  
  # Lower panel setup
  par(mar=c(mar_scale*4,mar_scale*4,0.1,0.1))
  
  if (is.null(yrange_lower)) {
    yrr=range(yr,na.rm=TRUE); if (!is.finite(sum(yrr))) yrr=c(-1,1)
  } else yrr=yrange_lower
  plot(0,xlim=xrange,ylim=yrr,type="n",
       xlab="Cutoffs",ylab=expression(paste(Delta,"prob")),xaxt="n",las=1)
  
  # Draw lines and confidence intervals on lower panel
  for (i in 1:length(objs)) {
    lines(xvals,objs[[i]]$y-objs[[1]]$y,col=col[i],lty=lty[i])
    
    # CI
    ci_dif=sqrt(objs[[i]]$ci^2 + objs[[1]]$ci^2); if (i==1) ci_dif=rep(0,length(objs[[i]]$ci))
    yci=c(objs[[i]]$y - objs[[1]]$y + ci_dif,rev(objs[[i]]$y - objs[[1]]$y - ci_dif)); wci=which(is.finite(yci))
    cic=col2rgb(ci_col[i]); cicx=rgb(cic[1]/256,cic[2]/256,cic[3]/256,alpha=0.5)
    polygon(c(xvals,rev(xvals))[wci],yci[wci],col=cicx,border=NA)
  }
  
  # Add highlight line
  if (!is.null(highlight)) {
    if (logscale) abline(v=log(highlight)) else abline(v=highlight)
  }
  
  # Add x axis labels
  pz=pretty(range(objs[[1]]$x),n=10)
  if (logscale) {
    axis(1,at=log(pz),labels=round(100*pz))
  } else {
    axis(1,at=pz,labels=round(100*pz))
  }
}









##' plot_decomp
##' 
##' Plots a bar graph of decomposition of FORP by cause of admission
##' 
##' Takes two matrices as input with the following specifications:
##'  Each matrix corresponds to one group
##'  Columns are named with the admission types to be plotted. Any admission 
##'   types including the string 'Died' are counted as deaths
##'  If the matrix has N rows, these are interpreted as corresponding to N
##'   score quantiles in increasing order.
##'  The (i,j)th entry of the matrix is the number of people admitted for 
##'   reason i with a score greater than or equal to (j-1)/N and less than (j/N) 
##'   who are in that group
##' 
##' @param decomp1 matrix for first group; see specification in description
##' @param decomp2 matrix for second group; see specification in description
##' @param threshold score threshold to plot (between 0 and 1)
##' @param labels labels for group 1 and group 2
##' @param inc_died set to TRUE to include a second panel showing 'death' type admissions
##' @param mar_scale scale margins by this amount. Also scales legend.
##' @return No return value, draws a figure
##' @export
##' @examples 
##' 
##' # See vignette
plot_decomp=function(decomp1,decomp2,threshold,labels,
                     inc_died=TRUE,mar_scale=1) {
  
  nquant=dim(decomp1)[1] # Number of score quantiles
  
  tq=floor(nquant*threshold) # Largest quantile less than score threshold

  dat1=colSums(decomp1[1:tq,]) # Group 1 with score < threshold
  dat2=colSums(decomp2[1:tq,]) # Group 2 with score < threshold
  ref1=colSums(decomp1) # Group 1, all
  ref2=colSums(decomp2) # Group 2, all
  all=colSums(rbind(decomp1,decomp2)) # Both groups, all
  
  dat=rbind(all,dat1,dat2,ref1,ref2)
  tot=outer(rowSums(dat),rep(1,dim(dat)[2]))
  pdat=dat/pmax(tot,1)
  ci=2*sqrt(pdat*(1-pdat)/tot)
  ci[which(is.na(ci))]=0
    
  g2=grep("Died",colnames(dat),value=TRUE)
  g1=setdiff(colnames(dat),g2)
  
  # Remove nonsense
  g2=setdiff(g2,"Died.of.unrecorded")
  g1=setdiff(g1,c("Not.recorded","Other"))
  
  w=order(-as.numeric(dat[1,g1]))
  w2=match(paste0("Died.of.",g1[w]),g2)
  
  oldpar = par(no.readonly = TRUE)
  on.exit(par(oldpar)) 
  
  
  lsp=0.1
  if (inc_died) {
    par(mfrow=c(2,1),mar=c(1,4*mar_scale,1,1))
  } else par(mar=c(8*mar_scale,4*mar_scale,1,1))
  
  plot(0,type="n",xlim=c(1,length(g1)),ylim=c(0,max(pdat)),xlab="",ylab="Prop. (adm)",xaxt="n")
  for (i in 1:dim(dat)[1]) {
    x=1:length(g1) + i*lsp
    y=as.numeric(pdat[i,g1[w]])
    yci=as.numeric(ci[i,g1[w]])
    segments(x,
             rep(0,length(g1)),
             x,
             y,col=i)
    segments(x, # Confidence intervals as thicker lines
             y-yci,
             x,
             y+yci,col=i,lwd=3)
    
    # Legend
    par(ps=par()$ps/mar_scale)
    legend("topright",c(
      "All adm.", 
      paste0(labels[1],"; score<",signif(threshold,digits=2)),
      paste0(labels[2],"; score<",signif(threshold,digits=2)),
      paste0("All adm. ",labels[1]),paste0("All adm. ",labels[2])),
      lty=1,col=1:5,cex=mar_scale)
    par(ps=par()$ps*mar_scale)
    
  }
  if (!inc_died) {
    axis(1,at=(1:length(g1)) + dim(dat)[1]*lsp/2,labels = gsub("."," ",g1[w],fixed=TRUE),las=2,cex.axis=0.7)
  } 
  if (inc_died) {
    par(mar=c(8*mar_scale,4*mar_scale,1,1))
    
    plot(0,type="n",xlim=c(1,length(g2)),ylim=c(0,max(pdat)/10),xlab="",ylab="Prop. (died)",xaxt="n")
    for (i in 1:dim(dat)[1]) {
      x=1:length(g2) + i*lsp
      y=as.numeric(pdat[i,g2[w2]])
      yci=as.numeric(ci[i,g2[w2]])
      segments(x,
               rep(0,length(g2)),
               x,
               y,col=i)
      segments(x, # Confidence intervals as thicker lines
               y-yci,
               x,
               y+yci,col=i,lwd=3)
    }
    axis(1,at=(1:length(g1)) + dim(dat)[1]*lsp/2,labels = gsub("."," ",g1[w],fixed=TRUE),las=2,cex.axis=0.7)
    
  }

}



##' for_breakdown
##' 
##' For a given category (e.g., 'male', 'over 65') considers 
##'   1) all admissions for people in that category
##'   2) all admissions for people in that category for which the 
##'       SPARRA score was less than some threshold (e.g., false
##'       negatives
##' 
##' For each of these groups, we consider the breakdown of medical
##'  admission types. We then plot the frequency of admission types
##'  in group 1 against the difference in frequencies between group 
##'  1 and group 2 (group 2 minus group 1).
##' An admission type which is relatively more common in group (1)
##'  indicates that, in the relevant category, the admission type
##'  tends to be associated with higher SPARRA scores (and is in a
##'  sense easier to predict). Such admission types will correspond
##'  to points below the line y=0.
##'  Admission types which are relatively more common in group 2 
##'  correspond to those which are relatively harder to predict. 
##'  These correspond to points above the line y=0
##' Since points are close together, only those greater than a 
##'  certain distance from 0 are marked.
##' 
##' Takes as an argument a matrix in which
##'   The matrix shows only data for the group in question
##'  Columns are named with the admission types to be plotted. Any admission 
##'   types including the string 'Died' are counted as deaths
##'  If the matrix has N rows, these are interpreted as corresponding to N
##'   score quantiles in increasing order.
##'  The (i,j)th entry of the matrix is the number of people admitted for 
##'   reason i with a score greater than or equal to (j-1)/N and less than (j/N) 
##'   who are in that group
##' 
##' @param decomp_table matrix for group; see specification in description
##' @param group name of group
##' @param threshold cutoff, rounded to nearest 0.05
##' @param inc_died set to TRUE to include a second panel showing 'death' type admissions
##' @param ldiff specifically label points this far from xy line
##' @param ci set to a value <1 to draw confidence intervals at that value, or FALSE to not draw confidence intervals.
##' @param xlimit limits for x axis; default c(-0.05,0.35)
##' @param ylimit limits for y axis; default c(-0.04,0.04)
##' @return ggplot figure (invisible)
##' @export
##' @examples
##' 
##' # See vignette
for_breakdown = function(decomp_table, group, threshold,
                         inc_died=TRUE,ldiff = 0.005,ci=0.95,
                         xlimit=c(-0.05,0.35),ylimit=c(-0.04, 0.04)) {
  
  # Remove 'died of' type admissions if specified
  if (!inc_died) decomp_table=decomp_table[,-grep("Died",colnames(decomp_table))]
  
  # Remove junk
  junk=c("Died.of.unrecorded","Not.recorded","Other")
  decomp_table=decomp_table[,setdiff(colnames(decomp_table),junk)]
  
  
  # Specify colours according to most common admission types
  pcol = c(1:8, rep(1, dim(decomp_table)[2] - 8))
  pcol = pcol[(rank(-colSums(decomp_table)))]
  
  nquant=dim(decomp_table)[1] # Number of score quantiles
  tq=floor(nquant*threshold) # Largest quantile less than score threshold
  
  # Frequencies by admission type
  xt = colSums(decomp_table)
  yt = colSums(decomp_table[1:tq,])
  
  # Sums (nxt=total admissions, nyt=total admissions with SPARRA<0.1)
  nxt = sum(xt)
  nyt = sum(yt)
  
  # Proportions of each admission type
  pxt = xt/max(1,nxt)
  pyt = yt/max(1,nyt)
  
  # Standard errors for pxt and pyt, taking nxt and nyt as constant
  se_pxt = sqrt(pxt * (1 - pxt)/nxt)
  se_pyt = sqrt(pyt * (1 - pyt)/nyt)
  
  # Standard error for difference pyt - pxt
  se_diff = sqrt(se_pxt^2 + se_pyt^2)
  
  # (asymptotic) confidence interval for difference pyt - pxt
  if (ci) {
    ci_level = ci
    z_alpha = -qnorm((1 - ci_level)/2)
    ci_diff_lower = (pyt - pxt) - z_alpha * se_diff
    ci_diff_upper = (pyt - pxt) + z_alpha * se_diff
  }
  
  # Normalise
  xt=xt/sum(xt)
  yt=yt/sum(yt)
  
  names(xt) = gsub(".", " ",names(xt),fixed=TRUE)
  data_df = data.frame(xt, yt, ci_diff_lower, ci_diff_upper)
  
  # Set up plot
  p = ggplot(data_df, aes(x = xt, y = yt - xt,color = factor(pcol))) +
    geom_hline(yintercept=0,linewidth=1/2) + 
    geom_point(pch = 16, size = 2)
  
  # Add confidence intervals
  if (ci){
    p = p +  geom_errorbar(
      aes(ymin = ci_diff_lower, ymax = ci_diff_upper),
      color = "blue",
      linewidth = 0.5)
  }
  
  
  p=suppressWarnings( p + 
    geom_label_repel(size = 2, label = ifelse(abs(xt - yt) > ldiff, names(xt), ""), nudge_x = 0.05, nudge_y = 0, box.padding = 0.5) +
    xlim(xlimit) + 
    ylim(ylimit) +
    labs(x = "Freq. in all adm.",
         y = paste0("(Freq. with score < ", threshold, ") - (Freq. in all adm.)"),
         title = paste0("Group: ", group)) +
    theme_minimal() +
    guides(color = "none"))
  suppressWarnings(print(p))
  invisible(p)
}




##**************************************************#
## Diagrams                                      ####
##**************************************************#


##' drawprop
##' 
##' Illustrates a proportion as a set of people who are blue rather than red.
##' 
##' Why anyone would want to think about a proportion this way is beyond the understanding of the authors, but the people have spoken.
##' 
##' @param prop the proportion to illustrate
##' @param ci half the 95% CI width of the proportion. 
##' @param nxy illustrate on an n x n grid
##' @param col1 colour to put the 'in' proportion
##' @param col2 the other colour
##' @param ... passed to 'plot'
##' @return No return value, draws a figure
##' @export
##' @examples 
##' # See vignette
drawprop=function(prop,ci,nxy=10,col1="maroon",col2="lightblue",...) {
  nprop=round(prop*nxy^2)
  nci_low=max(1,round((prop-ci)*nxy^2))
  nci_high=min(nxy^2,round((prop+ci)*nxy^2))
  plot(0,type="n",xlim=c(0,nxy),ylim=c(0,nxy),ann=FALSE,bty="n",xaxt="n",yaxt="n",...)
  segments(0:nxy,rep(0,1+nxy),0:nxy,rep(nxy,1+nxy),col="gray")
  segments(rep(0,1+nxy),0:nxy,rep(nxy,1+nxy),0:nxy,col="gray")
  inc=1
  if (nci_low<nci_high) ci_col=colorRampPalette(c(col1,col2))(nci_high-nci_low+1) else ci_col=col1
  for (i in nxy:1) for (j in 1:nxy) {
    if (inc<nci_low) xcol=col1
    if (inc>nci_high) xcol=col2
    if (inc>= nci_low & inc<= nci_high) xcol=ci_col[inc-nci_low+1]
    drawperson(j-0.5,i-0.5,col=xcol)
    inc=inc+1
  }
}


##' drawperson
##' 
##' Draws a figure at a particular location. With defaults, has centre at origin and fits in 1x1 box.
##' 
##' Dimensions customisable
##' 
##' @name drawperson
##' @description Draws a simple stock image of a person.
##' @param xloc x-axis offset from origin
##' @param yloc y-axis offset from origin
##' @param scale scale upwards from 1x1 box
##' @param headsize head size
##' @param headangle half angle of neck in terms of head
##' @param headloc location of centre of head relative to origin with scale 1
##' @param necklength neck length
##' @param shoulderwidth shoulder width
##' @param shouldersize size radius of arc for shoulder
##' @param armlength arm length
##' @param armangle angle of arm from horizontal
##' @param armwidth width of arm
##' @param leglength leg length
##' @param legangle angle of leg from horizontal
##' @param legwidth width of leg
##' @param torsolength length of torso
##' @param ... other parameters passed to polygon()
##' @export
##' @return invisibly returns co-ordinates
##' @examples
##' plot(0,xlim=c(-1,1),ylim=c(-1,1),type="n")
##' drawperson(0,0,1,col="yellow",border="red",lwd=3,lty=2)
drawperson=function(xloc=0,yloc=0,scale=1,
                    
                    headsize=0.16,
                    headangle=pi/8,
                    headloc=0.5,
                    
                    necklength=0.1,
                    shoulderwidth=0.1,
                    shouldersize=0.05,
                    
                    armlength=0.4,
                    armangle=7*pi/8,
                    armwidth=0.08,
                    
                    leglength=0.5,
                    legangle=9*pi/10,
                    legwidth=0.15,
                    
                    torsolength=0.4,
                    
                    ...
) {
  
  # head
  t_head=seq(0,pi-headangle,length=50)
  xv=headsize*sin(t_head)
  yv=headsize*cos(t_head) + headloc
  
  
  # neck
  xv=c(xv,xv[length(xv)])
  yv=c(yv,yv[length(yv)]-necklength)
  
  # shoulder
  xv=c(xv,xv[length(xv)] + shoulderwidth)
  yv=c(yv,yv[length(yv)])
  t_shoulder=seq(0,2*pi/5,length=20)
  xv=c(xv,shouldersize*sin(t_shoulder) + xv[length(xv)])
  yv=c(yv,shouldersize*cos(t_shoulder) + yv[length(yv)] - shouldersize)
  
  # arm
  xv=c(xv,xv[length(xv)] + armlength*sin(armangle))
  xv=c(xv,xv[length(xv)] + armwidth*sin(armangle + pi/2))
  xv=c(xv,xv[length(xv)] - 0.8*armlength*sin(armangle))
  yv=c(yv,yv[length(yv)] + armlength*cos(armangle))
  yv=c(yv,yv[length(yv)] + armwidth*cos(armangle + pi/2))
  yv=c(yv,yv[length(yv)] - 0.8*armlength*cos(armangle))
  
  # torso
  legstart=c(legwidth*sin(legangle - pi/2),yv[length(yv)]-torsolength + legwidth*cos(legangle - pi/2))
  xv=c(xv,legstart[1])
  yv=c(yv,legstart[2])
  
  
  # leg
  xv=c(xv,xv[length(xv)] + leglength*sin(legangle))
  xv=c(xv,xv[length(xv)] + legwidth*sin(legangle + pi/2))
  xv=c(xv,xv[length(xv)] - leglength*sin(legangle))
  
  yv=c(yv,yv[length(yv)] + leglength*cos(legangle))
  yv=c(yv,yv[length(yv)] + legwidth*cos(legangle + pi/2))
  yv=c(yv,yv[length(yv)] - leglength*cos(legangle))
  
  # both sides
  xp=c(xv,rev(-xv))
  yp=c(yv,rev(yv))
  
  # scale to fit in box
  xp=xp*0.5
  yp=yp*0.5
  
  # scale
  xp=xloc + scale*xp
  yp=yloc + scale*yp
  
  polygon(xp,yp,...)
  invisible(cbind(xp,yp))
}


