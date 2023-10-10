##' sim_pop_data
##' 
##' Simulates population data with a reasonably realistic joint distribution
##' 
##' Simulates data for a range of people for the variables
##'  * Age (`age`)
##'  * Sex (`sexM`; 1 if male)
##'  * Race/ethnicity (`raceNW`: 1 if non-white ethnicity)
##'  * Number of previous hospital admissions (`PrevAdm`)
##'  * Deprivation decile (`SIMD`: 1 most deprived, 10 least deprived. NOTE - opposite to English IMD)
##'  * Urban-rural residence status (`urban_rural`: 1 for rural)
##'  * Mainland-island residence status (`mainland_island`: 1 for island)
##'  * Hospital admission (`target`: 1/TRUE if admitted to hospital in year following prediction date)
##' 
##' Can optionally add an ID column. 
##' 
##' Optionally includes an admission reason for samples with `target=1`. These admission reasons 
##'  roughly correspond to the first letters of ICD10 categories, and can either correspond to an
##'  admission or death. Admission reasons are simulated with a non-constant multinomial distribution
##'  which varies across age/sex/ethnicity/urban-rural/mainland-island/PrevAdm values in a randomly-
##'  chosen way. The distributions of admission reasons are *not* however chosen to reflect real 
##'  distributions, nor are systematic changes in commonality of admission types across categories
##'  intended to appear realistic.
##' 
##' 
##' @param npop population size
##' @param coef_adjust inverse scale for all (true) coefficients (default 4): lower means that hospital admissions are more predictable from covariates.
##' @param offset offset for logistic model (default 1): higher means a lower overall prevalence of admission
##' @param vcor a valid 5x5 correlation matrix (default NULL), giving correlation between variables. If 'NULL', values roughly represents realistic data.
##' @param coefs coefficients of age, male sex, non-white ethnicity, number of previous admissions, and deprivation decile on hospital admissions, Default (2,1,0,5,3). Divided through by coef_adjust.
##' @param seed random seed (default 12345)
##' @param incl_id include an ID column (default TRUE)
##' @param incl_reason include a column indicating reason for admission.
##' @export
##' @return data frame with realistic values.
##' @examples
##' 
##' # Simulate data
##' dat=sim_pop_data(10000)
##' cor(dat[,1:7])
##' 
##' # See vignette
sim_pop_data= function(
  npop,
  coef_adjust=4,
  offset=1,
  vcor=NULL,
  coefs=c(2,1,0,5,3,0,0),
  seed=12345,
  incl_id=TRUE,
  incl_reason=TRUE
) {
  
  # Seed
  set.seed(seed)
  
  
  # Variable names
  varnames=c("Age","SexM","RaceNW","PrevAdm","SIMD","urban_rural","mainland_island")
  
  
  
  # Scripts
  #source("~/Research/SPARRA/Results/Current/SPARRAv4/auxiliary.R")
  
  # Correlations: -10 to 10
  # E.g. cor(Age,SexM)=-3: higher age is less likely male
  #  cor(SexM,RaceNW)=3: individuals of non-white ethnicity more likely to be male
  if (is.null(vcor)) {
    vcor=rbind(
      # Age     SexM   RaceNW  PA    SIMD   UR   MI
      c(0,      -3,    -3,     6,    3,     5,   4  ),  # Age
      c(0,      0,     3,      1,    1,     2,   1  ),  # SexM
      c(0,      0,     0,      -2,   -2,    -4,  -5 ),  # RaceNW
      c(0,      0,     0,      0,    -3,    -2,  -1 ),  # Prev Admissions
      c(0,      0,     0,      0,    0,     5,   4  ),  # SIMD
      c(0,      0,     0,      0,    0,     0,   8  ),  # Urban-rural
      c(0,      0,     0,      0,    0,     0,   0  )   # Mainland-island
    )  
    vcor=vcor+ 12*diag(dim(vcor)[1])
    vcor=vcor + t(vcor) - vcor*diag(dim(vcor)[1]) # Symmetrise
    vcor=vcor/max(vcor*diag(dim(vcor)[1])) # Now a correlation matrix
    #print(eigen(vcor)$values) # check positive definite
  }
  
  rownames(vcor)=varnames
  colnames(vcor)=varnames
  
  # Outcome specification
  names(coefs)=varnames
  coefs=coefs/coef_adjust
  
  # Now sample latent variables. All marginals are N(0,1)
  sp_l=rmvnorm(npop,mean = rep(0,dim(vcor)[1]),sigma=vcor)
  colnames(sp_l)=varnames
  
  # Outcome
  z_out=(sp_l %*% coefs)
  p_out=logistic(z_out-offset)
  target=rbinom(npop,1,prob=p_out)
  
  # Transform to realistic values
  
  # Age
  lat_age=pnorm(sp_l[,"Age"]) # Now uniform between 0 and 1
  age=round(90*tan(lat_age*0.8)) # Roughly realistic age distribution
  
  # SexM
  lat_sexM=sp_l[,"SexM"]
  sexM=(lat_sexM>0)
  
  # RaceNW: in Scotland about 8%
  lat_raceNW=sp_l[,"RaceNW"]
  raceNW=(lat_raceNW>qnorm(1-0.08))
  
  # Prev admissions
  lat_prevAdm=sp_l[,"PrevAdm"]
  prevAdm=round(tan(pnorm(lat_prevAdm)*1.5))
  
  # SIMD
  lat_simd=sp_l[,"SIMD"]
  SIMD=round(0.5 + 10*(pnorm(lat_simd)))
  
  # Urban_rural
  lat_ur=sp_l[,"urban_rural"]
  urban_rural=(lat_ur>1)
  
  # Mainland-island
  lat_mi=sp_l[,"mainland_island"]
  mainland_island=(lat_ur>2)
  
  # Design matrix
  dmat=data.frame(age,sexM,raceNW,prevAdm,SIMD,urban_rural,mainland_island,target)
  
  # ID
  if (incl_id) dmat$id=1:dim(dmat)[1]
  
  # Reason
  if (incl_reason) {
    
    dprob=0.2 # This proportion of admissions were deaths
    
    cats0=c("Abnormality.NEC", "Blood", "Circulatory", "Congenital", "Digestive", 
            "Ear", "Endocrine.metabolic", "External", "Eye", "Genitourinary", 
            "Infectious.disease", "Mental.behavioural", "Musculoskeletal", 
            "Neoplasm", "Nervous.system", "Not.recorded", "Obstetric.puerperium", 
            "Other", "Perinatal", "Respiratory", "Skin")
    cats=c(cats0,paste0("Died.of.",cats0))
    
    # Order dmat by variables in a randomly-chosen way, then use a different 
    #  multinomial distribution for each successive third of the dataset. This means
    #  that age/sex/etc are associated with slightly different multinomial 
    #  distributions of admission causes.
    o_coef=rnorm(7); 
    dmat2=dmat[,1:7]; for (i in 1:dim(dmat2)[2]) dmat2[,i]=(dmat2[,i]-mean(dmat2[,i]))/sd(dmat2[,i])
    xord=(as.matrix(dmat2) %*% o_coef) + rnorm(dim(dmat2)[1])
    ox=order(xord)
    
    b1=round(npop/3); b2=round(2*npop/3)
    q1=ox[1:b1]; q2=ox[(b1+1):b2]; q3=ox[(b2+1):npop]
    
    adm_reason=rep(NA,dim(dmat)[1])
    for (i in 1:3) {
      q=get(paste0("q",i))
      vx0=runif(length(cats0)) # (Relative) multinomial probabilities for admision types
      vx=c(vx0*(1-dprob),vx0*dprob) # Relative multinomial probabilities for admission types and deaths
      vx=vx/sum(vx) # Technically unnecessary normalisation
      reason_i=sample(cats,length(q),prob=vx,replace=TRUE)
      adm_reason[q]=reason_i
    }
    adm_reason[which(dmat$target==0)]=NA
    dmat$reason=adm_reason
  }
  
  return(dmat)
  
}




##' dat2mat
##' 
##' Generates matrices for decomposition of admission type which can be used in `plot_decomp`
##' 
##' Generates two matrices with the following specifications:
##'  Each matrix corresponds to one group
##'  Columns are named with the admission types to be plotted. Any admission 
##'   types including the string 'Died' are counted as deaths
##'  If the matrix has N rows, these are interpreted as corresponding to N
##'   score quantiles
##'  The (i,j)th entry of the matrix is the number of people admitted for 
##'   reason i with a score greater than or equal to (j-1)/N and less than (j/N) 
##'   who are in that group
##' 
##' @param dat data frame with population data, such as output from sim_pop_data. Must include a column `reason`
##' @param score risk scores corresponding to `dat`
##' @param group1 indices for group 1
##' @param group2 indices for group 2
##' @param nquant number of quantiles of code to use; default 20
##' @param cats vector of strings giving names of admission categories; default the unique values in dat$reason. Can include NAs.
##' @export
##' @return list with two objects `matrix1` and `matrix2` giving output matrices
##' @examples
##' # See vignette
dat2mat=function(dat,score,group1,group2,nquant=20,cats=unique(dat$reason)) {

  cats=setdiff(cats,NA)
  
  out1=matrix(0,nquant,length(cats))
  out2=matrix(0,nquant,length(cats))
  colnames(out1)=cats
  colnames(out2)=cats
  
  for (i in 1:nquant) {
    sub=which((score >= (i-1)/nquant)  & (score < i/nquant) & (dat$target==1))
    sub1=intersect(sub,group1)
    sub2=intersect(sub,group2)
    if (length(sub1)>0) out1[i,]=table(dat$reason[sub1])[cats]
    if (length(sub2)>0) out2[i,]=table(dat$reason[sub2])[cats]
  }
  out1[which(is.na(out1))]=0
  out2[which(is.na(out2))]=0
  
  return(list(matrix1=out1,matrix2=out2))
}