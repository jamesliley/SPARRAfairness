## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo=TRUE----------------------------------------------------------------

# Load packages
library(SPARRAfairness)
library(ranger)

# Get data
data(all_data)
data(decomposition_matrix)

# Set random seed
seed=463825
set.seed(seed)

# Simulate data
pop_data=sim_pop_data(10000)

# First few rows
head(pop_data)

## ----echo=TRUE----------------------------------------------------------------

# Fit model
sim_model=ranger(
  target~age+sexM+raceNW+prevAdm+SIMD+urban_rural+mainland_island,
  data=pop_data)

# Model predictions
score=predict(sim_model,pop_data)$predictions

# Affix to pop_data
pop_data$score=score

## ----echo=TRUE----------------------------------------------------------------
group1=which(pop_data$urban_rural==FALSE)
group2=which(pop_data$urban_rural==TRUE)

## ----echo=TRUE----------------------------------------------------------------
score_cutoffs=seq(0,1,length=101)

## ----echo=TRUE----------------------------------------------------------------
dem_par=demographic_parity(pop_data$score,group1,group2,cutoffs=score_cutoffs)

## ----echo=TRUE----------------------------------------------------------------

# Counterfactual sample (samples are rural, but 'resemble' urban samples)
cf_rural_ids=counterfactual_yhat(dat=pop_data,X=c("age","sexM"),
                                             x=NULL,G="urban_rural",
                                             g=FALSE,gdash=TRUE,
                                             excl=c("id","score","target"))

# Put together into data frame
cf_pop_data=rbind(pop_data[group1,],pop_data[cf_rural_ids,])
cf_group1=which(cf_pop_data$urban_rural==FALSE); 
cf_group2=which(cf_pop_data$urban_rural==TRUE)

cf_dem_par=demographic_parity(cf_pop_data$score,cf_group1,
                                 cf_group2,cutoffs=score_cutoffs)

## ----echo=TRUE----------------------------------------------------------------
roc1=getroc(pop_data$target[group1],pop_data$score[group1])
roc2=getroc(pop_data$target[group2],pop_data$score[group2])

## ----echo=TRUE----------------------------------------------------------------
prc1=getprc(pop_data$target[group1],pop_data$score[group1])
prc2=getprc(pop_data$target[group2],pop_data$score[group2])

## ----echo=TRUE----------------------------------------------------------------
cal1=getcal(pop_data$target[group1],pop_data$score[group1])
cal2=getcal(pop_data$target[group2],pop_data$score[group2])

## ----echo=TRUE----------------------------------------------------------------
spec_for=c(NA,1,NA,0,NA,1)
for_12=group_fairness(spec_for,pop_data$score,pop_data$target,
                      group1,group2,cutoffs = score_cutoffs)

## ----echo=TRUE----------------------------------------------------------------
spec_fdr=c(NA,0,NA,1,NA,1)
fdr_12=group_fairness(spec_fdr,pop_data$score,pop_data$target,
                      group1,group2,cutoffs = score_cutoffs)

## ----echo=TRUE----------------------------------------------------------------
cat_12=paste0(pop_data$age,"_",pop_data$sexM)
for_adjusted_12=adjusted_for(pop_data$score,pop_data$target,cat_12,
                             group1,group2,cutoffs=score_cutoffs,nboot=10)
fdr_adjusted_12=adjusted_fdr(pop_data$score,pop_data$target,cat_12,
                             group1,group2,cutoffs=score_cutoffs,nboot=10)

## ----echo=TRUE----------------------------------------------------------------
alpha=0.05 # We will plot 1-alpha confidence intervals
q=-qnorm(alpha/2) # Multiply standard errors by this to get half CI width
highlight_value=0.5 # Highlight what happens at a cutoff of 0.5/50%
colour_scheme <- phs_colours(c("phs-blue", 
                               "phs-purple", 
                               "phs-magenta")) # set colour scheme

## ----echo=TRUE----------------------------------------------------------------
obj_list=list(
  list(x=score_cutoffs,y=dem_par[1,],ci=q*dem_par[2,]),
  list(x=score_cutoffs,y=dem_par[3,],ci=q*dem_par[4,])
)
groupmetric_2panel(obj_list,
                   labels=c("Urban","Rural"),
                   col=phs_colours(c("phs-blue","phs-magenta")),
                   ci_col=phs_colours(c("phs-blue","phs-magenta")),
                   yrange=c(0,1),
                   lpos="topleft",
                   yrange_lower=c(-0.2,0.2),
                   highlight=highlight_value,
                   logscale=TRUE)

## ----echo=TRUE----------------------------------------------------------------
obj=all_data$dp_v3_Urban_rural_all

obj_list=list(
  list(x=exp(obj$xx),y=obj$yA,ci=obj$cA),
  list(x=exp(obj$xx),y=obj$yB,ci=obj$cB)
)
groupmetric_2panel(obj_list,
                   labels=c("Urban","Rural"),
                   col=phs_colours(c("phs-blue","phs-magenta")),
                   ci_col=phs_colours(c("phs-blue","phs-magenta")),
                   yrange=c(0,1),
                   lpos="topleft",
                   yrange_lower=c(-0.2,0.2),
                   highlight=highlight_value,
                   logscale=TRUE)


## ----echo=TRUE----------------------------------------------------------------
obj_list=list(
  list(x=score_cutoffs,y=cf_dem_par[1,],ci=q*cf_dem_par[2,]),
  list(x=score_cutoffs,y=cf_dem_par[3,],ci=q*cf_dem_par[4,])
)
groupmetric_2panel(obj_list,
                   labels=c("Urban","Rural"),
                   col=phs_colours(c("phs-blue","phs-magenta")),
                   ci_col=phs_colours(c("phs-blue","phs-magenta")),
                   yrange=c(0,1),
                   lpos="topleft",
                   yrange_lower=c(-0.2,0.2),
                   highlight=highlight_value,
                   logscale=TRUE)



## ----echo=TRUE----------------------------------------------------------------
obj=all_data$counterfactual_dp_v3_Urban_rural

obj_list=list(
  list(x=exp(obj$xx),y=obj$yA,ci=obj$cA),
  list(x=exp(obj$xx),y=obj$yB,ci=obj$cB)
)
groupmetric_2panel(obj_list,
                   labels=c("Urban","Rural"),
                   col=phs_colours(c("phs-blue","phs-magenta")),
                   ci_col=phs_colours(c("phs-blue","phs-magenta")),
                   yrange=c(0,1),
                   lpos="topleft",
                   yrange_lower=c(-0.2,0.2),
                   highlight=highlight_value,
                   logscale=TRUE)


## ----echo=TRUE----------------------------------------------------------------
roc_2panel(list(roc1,roc2),
           labels = c("Urban","Rural"),
           col=colour_scheme,
           highlight=highlight_value,
           yrange_lower=c(-0.2,0.2))

## ----echo=TRUE----------------------------------------------------------------
roc_2panel(all_data$roc_v3_Urban_rural[1:3],
           labels = c("Overall","Urban","Rural"),
           col=colour_scheme,
           highlight=highlight_value,
           yrange_lower=c(-0.2,0.2))

## ----echo=TRUE----------------------------------------------------------------
prc_2panel(list(prc1,prc2),
           labels = c("Urban","Rural"),
           col=colour_scheme,
           highlight=highlight_value,
           yrange_lower=c(-0.2,0.2))

## ----echo=TRUE----------------------------------------------------------------
prc_2panel(all_data$prc_v3_Urban_rural[1:3],
           labels = c("Overall","Urban","Rural"),
           col=colour_scheme,
           highlight=highlight_value,
           yrange_lower=c(-0.2,0.2))

## ----echo=TRUE----------------------------------------------------------------
cal_2panel(list(cal1,cal2),
           labels = c("Urban","Rural"),
           col=colour_scheme,
           highlight=highlight_value,
           yrange_lower=c(-0.2,0.2))

## ----echo=TRUE----------------------------------------------------------------
cal_2panel(all_data$cal_v3_Urban_rural[1:3],
           labels = c("Overall","Urban","Rural"),
           col=colour_scheme,
           highlight=highlight_value,
           yrange_lower=c(-0.2,0.2))

## ----echo=TRUE----------------------------------------------------------------
obj_list=list(
  list(x=score_cutoffs,y=for_12[1,],ci=for_12[2,]),
  list(x=score_cutoffs,y=for_12[3,],ci=for_12[4,])
)
groupmetric_2panel(obj_list,
                   labels=c("Urban","Rural"),
                   col=phs_colours(c("phs-blue","phs-magenta")),
                   ci_col=phs_colours(c("phs-blue","phs-magenta")),
                   yrange=c(0,0.5),
                   lpos="topleft",
                   yrange_lower=c(-0.1,0.1),
                   highlight=highlight_value)

## ----echo=TRUE----------------------------------------------------------------
obj_list=list(
  list(x=score_cutoffs,y=fdr_12[1,],ci=fdr_12[2,]),
  list(x=score_cutoffs,y=fdr_12[3,],ci=fdr_12[4,])
)
groupmetric_2panel(obj_list,
                   labels=c("Urban","Rural"),
                   col=phs_colours(c("phs-blue","phs-magenta")),
                   ci_col=phs_colours(c("phs-blue","phs-magenta")),
                   yrange=c(0,1),
                   lpos="topleft",
                   yrange_lower=c(-0.1,0.1),
                   highlight=highlight_value)

## ----echo=TRUE----------------------------------------------------------------
obj=all_data$forp_v3_Urban_rural_all

obj_list=list(
  list(x=obj$cutoffs,y=obj$p_AA,ci=obj$ci_AA),
  list(x=obj$cutoffs,y=obj$p_BB,ci=obj$ci_BB)
)
groupmetric_2panel(obj_list,
                   labels=c("Urban","Rural"),
                   col=phs_colours(c("phs-blue","phs-magenta")),
                   ci_col=phs_colours(c("phs-blue","phs-magenta")),
                   yrange=c(0,0.1),
                   lpos="topleft",
                   yrange_lower=c(-0.02,0.02),
                   highlight=highlight_value)


## ----echo=TRUE----------------------------------------------------------------
obj=all_data$fdrp_v3_Urban_rural_all

obj_list=list(
  list(x=obj$cutoffs,y=obj$p_AA,ci=obj$ci_AA),
  list(x=obj$cutoffs,y=obj$p_BB,ci=obj$ci_BB)
)
groupmetric_2panel(obj_list,
                   labels=c("Urban","Rural"),
                   col=phs_colours(c("phs-blue","phs-magenta")),
                   ci_col=phs_colours(c("phs-blue","phs-magenta")),
                   yrange=c(0,1),
                   lpos="topright",
                   yrange_lower=c(-0.1,0.1),
                   highlight=highlight_value)


## ----echo=TRUE----------------------------------------------------------------
obj_list=list(
  list(x=score_cutoffs,y=for_adjusted_12[1,],ci=for_adjusted_12[2,]),
  list(x=score_cutoffs,y=for_adjusted_12[3,],ci=for_adjusted_12[4,])
)
groupmetric_2panel(obj_list,
                   labels=c("Urban","Rural"),
                   col=phs_colours(c("phs-blue","phs-magenta")),
                   ci_col=phs_colours(c("phs-blue","phs-magenta")),
                   yrange=c(0,0.5),
                   lpos="topleft",
                   yrange_lower=c(-0.1,0.1),
                   highlight=highlight_value)

## ----echo=TRUE----------------------------------------------------------------
obj_list=list(
  list(x=score_cutoffs,y=fdr_adjusted_12[1,],ci=fdr_adjusted_12[2,]),
  list(x=score_cutoffs,y=fdr_adjusted_12[3,],ci=fdr_adjusted_12[4,])
)
groupmetric_2panel(obj_list,
                   labels=c("Urban","Rural"),
                   col=phs_colours(c("phs-blue","phs-magenta")),
                   ci_col=phs_colours(c("phs-blue","phs-magenta")),
                   yrange=c(0,1),
                   lpos="topleft",
                   yrange_lower=c(-0.1,0.1),
                   highlight=highlight_value)

## ----echo=TRUE----------------------------------------------------------------
obj=all_data$forp_adjusted_v3_Urban_rural

obj_list=list(
  list(x=obj$cutoffs,y=obj$p_AA,ci=obj$ci_AA),
  list(x=obj$cutoffs,y=obj$p_BB,ci=obj$ci_BB)
)
groupmetric_2panel(obj_list,
                   labels=c("Urban","Rural"),
                   col=phs_colours(c("phs-blue","phs-magenta")),
                   ci_col=phs_colours(c("phs-blue","phs-magenta")),
                   yrange=c(0,0.1),
                   lpos="topleft",
                   yrange_lower=c(-0.02,0.02),
                   highlight=highlight_value)


## ----echo=TRUE----------------------------------------------------------------
obj=all_data$fdrp_adjusted_v3_Urban_rural

obj_list=list(
  list(x=obj$cutoffs,y=obj$p_AA,ci=obj$ci_AA),
  list(x=obj$cutoffs,y=obj$p_BB,ci=obj$ci_BB)
)
groupmetric_2panel(obj_list,
                   labels=c("Urban","Rural"),
                   col=phs_colours(c("phs-blue","phs-magenta")),
                   ci_col=phs_colours(c("phs-blue","phs-magenta")),
                   yrange=c(0,1),
                   lpos="topright",
                   yrange_lower=c(-0.1,0.1),
                   highlight=highlight_value)


## ----echo=TRUE,fig.width=8,fig.height=6---------------------------------------
cutoff=0.5 # 10% risk score threshold
decomp_matrices=dat2mat(pop_data,
                        score=pop_data$score,
                        group1=which(pop_data$urban_rural==0),
                        group2=which(pop_data$urban_rural==1),
                        nquant=20)
plot_decomp(decomp_matrices$matrix1,
            decomp_matrices$matrix2,
            threshold=cutoff,
            labels=c("Urban","Rural"))

## ----echo=TRUE,fig.width=8,fig.height=6---------------------------------------
cutoff=0.5 # 10% risk score threshold
names_group1=paste0("v3_Urban_q",1:20)
names_group2=paste0("v3_Rural_q",1:20)
decomp1=decomposition_matrix[names_group1,]
decomp2=decomposition_matrix[names_group2,]
plot_decomp(decomp1,
            decomp2,
            threshold=cutoff,
            labels=c("Urban","Rural"))

## ----echo=TRUE----------------------------------------------------------------
plot(0,type="n",bty="n",ann=F,xlim=c(-1,1)/2,ylim=c(-1,1)/2)
drawperson(0,0,col="blue",border="gray",lwd=3)

## ----echo=TRUE,fig.width=6,fig.height=6---------------------------------------
drawprop(0.3,ci=c(0.2,0.4),col1 = "blue",col2="gray")

