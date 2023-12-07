# Code for testing Type I error rate under the null
library(tidyverse)
library(locfit)
library(locpol)
library(pscl)
library(gridExtra)
library(rdd)

logistic = function(x) { exp(x)/(1+exp(x))}
loglinear = function(x) {exp(x)}

### local likelihood with IK bandwidth (for binomial) from Xu 2017

rd_mnl_final<-function(DAT,c,H0_t,H0_R,H0_r,level)
{
  n=length(DAT[,1])
  J=length(DAT[1,])-2
  
  dat=DAT[,1:2] #dat is tailored for the command "multinom"
  
  for (j in 1:J) {
    for (i in 1:n) {if (DAT[i,j+2]==1) dat[i,2]=j+1}
  }
  
  for (i in 1:n) {if (dat[i,1]>=c) break}
  
  L=dat[1:(i-1),]  # the left sample
  R=dat[i:n,]      # the right sample
  
  n_f=length(L[,1])
  n_r=length(R[,1])
  
  g=matrix(0,J,1)
  g2=matrix(0,J,1)
  g3=matrix(0,J,1)
  mu=matrix(0,J,1)
  Gamma=matrix(0,J,J)
  
  library(nnet)
  ss=relevel(factor(L[,2]),ref=J+1)
  #ss=relevel(factor(L[,2]),ref=3)
  res=multinom(ss ~ L[,1]+I(L[,1]^2)+I(L[,1]^3))
  res=summary(res)$coefficients
  res=matrix(res,nrow=J,ncol=4)
  
  sum=1
  
  sum
  for (j in 1:J) {
    g[j]=res[j,1]+res[j,2]*c+res[j,3]*c^2+res[j,4]*c^3
    g2[j]=3*2*res[j,4]*c+2*res[j,3] # 2nd derivative
    g3[j]=3*2*res[j,4] # 3rd derivative
    sum=sum+exp(g[j])
  }
  
  sum
  
  for (j in 1:J) {
    mu[j]=exp(g[j])/sum
  }
  
  for (j in 1:J){
    for (k in 1:J){
      if (j==k) Gamma[j,k]=mu[j]*(1-mu[j])
      else Gamma[j,k]=-mu[j]*mu[k]
    }
    
  }
  
  den=density(DAT[,1],from=c,to=c)
  den_c=den$y[1]
  
  K0=2.702^5 # the kernel constant (flat kernel is used) (for bandwidth to estimate mu)
  K2=3.557^7 # the kernel constant (flat kernel is used) (for bandwidth to estimate g2, the 2nd derivative of g)
  h_mu_f=K0*sum(diag(Gamma))/(n*den_c*(norm(Gamma%*%g2,type="F"))^2)
  h_mu_f=h_mu_f^(1/5)
  h_g2_f=K2*sum(diag(solve(Gamma)))/(n*den_c*(norm(g3,type="F"))^2)
  h_g2_f=h_g2_f^(1/7)
  
  ss=relevel(factor(R[,2]),ref=J+1)
  res=multinom(ss ~ R[,1]+I(R[,1]^2)+I(R[,1]^3))
  res=summary(res)$coefficients
  res=matrix(res,nrow=J,ncol=4)
  
  sum=1
  for (j in 1:J) {
    g[j]=res[j,1]+res[j,2]*c+res[j,3]*c^2+res[j,4]*c^3
    g2[j]=3*2*res[j,4]*c+2*res[j,3] # 2nd derivative
    g3[j]=3*2*res[j,4] # 3rd derivative
    sum=sum+exp(g[j])
  }
  
  for (j in 1:J) mu[j]=exp(g[j])/sum
  
  
  Gamma=matrix(0,J,J)
  for (j in 1:J) {
    for (k in 1:J) {
      if (j==k) Gamma[j,k]=mu[j]*(1-mu[j])
      else Gamma[j,k]=-mu[j]*mu[k]
    }
  }
  
  h_mu_r=K0*sum(diag(Gamma))/(n*den_c*(norm(Gamma%*%g2,type="F"))^2)
  h_mu_r=h_mu_r^(1/5)
  h_g2_r=K2*sum(diag(Gamma^(-1)))/(n*den_c*(norm(g3,type="F"))^2)
  h_g2_r=h_g2_r^(1/7)
  
  for (i in 1:n_f) {if (L[i,1]>=c-h_mu_f) break}
  L1=L[i:n_f,] # the left sample used to estimate mu
  for (i in 1:n_f) {if (L[i,1]>=c-h_g2_f) break}
  L2=L[i:n_f,] # the left sample used to estimate the 2nd derivative of g
  
  for (i in 1:n_r) {if (R[i,1]>=c+h_mu_r) break}
  R1=R[1:i,] # the right sample used to estimate mu
  for (i in 1:n_f) {if (R[i,1]>=c+h_g2_r) break}
  R2=R[1:i,] # the right sample used to estimate the 2nd derivative of g
  
  # now we use these 4 subsamples to estimate Gamma_f (using L1), Gamma_r (using R1), g2_f (using L2) and g2_r (using R2).
  
  ss=relevel(factor(L1[,2]),ref=J+1)
  res=multinom(ss ~ L1[,1])
  res=summary(res)$coefficients
  res=matrix(res,nrow=J,ncol=2)
  
  sum=1
  for (j in (1:J)) {
    g[j]=res[j,1]+res[j,2]*c
    sum=sum+exp(g[j])
  }
  
  for (j in 1:J) mu[j]=exp(g[j])/sum
  
  Gamma_f=matrix(0,J,J)
  for (j in 1:J) {
    for (k in 1:J) {
      if (j==k) Gamma_f[j,k]=mu[j]*(1-mu[j])
      else Gamma_f[j,k]=-mu[j]*mu[k]
    }
  }
  
  ss=relevel(factor(R1[,2]),ref=J+1)
  res=multinom(ss ~ R1[,1])
  res=summary(res)$coefficients
  res=matrix(res,nrow=J,ncol=2)
  
  sum=1
  for (j in 1:J)
  {
    g[j]=res[j,1]+res[j,2]*c;
    sum=sum+exp(g[j]);
    mu[j]=exp(g[j]);
  }
  
  for (j in 1:J) {mu[j]=exp(g[j])/sum}
  
  Gamma_r=matrix(0,J,J)
  for (j in 1:J) {
    for (k in 1:J) {
      if (j==k) Gamma_r[j,k]=mu[j]*(1-mu[j])
      else Gamma_r[j,k]=-mu[j]*mu[k]
    }
  }
  
  g2_f=matrix(0,J,1)
  g2_r=matrix(0,J,1)
  ss=relevel(factor(L2[,2]),ref=J+1)
  res=multinom(ss ~ L2[,1]+I(L2[,1]^2))
  res=summary(res)$coefficients
  res=matrix(res,nrow=J,ncol=3)
  
  for (j in 1:J) {g2_f[j]=2*res[j,3]}
  
  ss=relevel(factor(R2[,2]),ref=J+1)
  res=multinom(ss ~ R2[,1]+I(R2[,1]^2))
  res=summary(res)$coefficients
  res=matrix(res,nrow=J,ncol=3)
  
  for (j in 1:J) {g2_r[j]=2*res[j,3]}
  
  # now to calculate one final product: optimal bandwidth
  
  h_opt=2.702^5*sum(diag(Gamma_f+Gamma_r))/(n*den_c*(norm(Gamma_r%*%g2_r-Gamma_f%*%g2_f,type="F"))^2)
  h_opt=h_opt^(0.2)
  
  
  #Reporting all results
  list(h_opt=h_opt)
}

#For empirical researchers

#Step 1: Input the data, naming it as "DAT". (e.g. DAT has dimension 40750 by 4, as in the CA infant birth weight example)
#Step 2: Define the cutoff c (e.g. c=1.5, as in the infant birth weight example)
#Step 3: Run all the lines above (i.e. define a function "rd_mnl_final" in R).
#Step 4: Run the following:





### Local linear regression case

one_linear_ft_coef_and_var=function(data, value, bw, dire='left'){
  
  if (dire=='left'){
    result=data %>% filter(x<value & x>value-bw)
  }
  else if (dire=='right'){
    result=data %>% filter(x>=value & x<value+bw)
  }
  else {
    result='error : dire must be left or right'
  }
  
  if(dim(result)[1] < 2){
    return(list(NA, NA))
  }else{
    K= dnorm(x=value, mean=result$x, sd=bw)
    data_revision = result %>% mutate(x=x-value)
    est=lm(y~x, data=data_revision)
    
    
    return(list(est$coefficients, vcov(est)))
  }
  
}


# compute tau_SRD by local linear regression



SRD_lm_coef_var=function(data, value, bw){
  
  
  result_left = one_linear_ft_coef_and_var(data=data, value=value, bw=bw, dire='left')
  result_right = one_linear_ft_coef_and_var(data=data, value=value, bw=bw, dire='right')
  
  
  tauhat=result_right[[1]][1]-result_left[[1]][1]
  
  
  sigma_left=result_left[[2]]
  sigma_right=result_right[[2]]
  
  
  
  nl=sum(data$x<value & data$x>value-bw)
  nr=sum(data$x>value & data$x<value+bw)
  
  
  asym_var=sigma_left[1, 1]+sigma_right[1, 1]
  
  
  return(c(tauhat, asym_var))
  
  
}


### GLM case

# local likelihood with glm

one_glm_ft_coef_and_var = function(data, value, bw, dire='left', type){
  
  
  if (dire=='left'){
    result=data %>% filter(x<value & x>value-bw)
  }
  else if (dire=='right'){
    result=data %>% filter(x>=value & x<value+bw)
  }
  else {
    result='error : dire must be left or right'
  }
  
  K= dnorm(x=value, mean=result$x, sd=bw)
  data_revision = result %>% mutate(x=x-value)
  est=glm(y~x, data=data_revision, family = type)
  
  
  return(list(est$coefficients, vcov(est)))
  
  
}



# compute tau_SRD by local likelihood



SRD_glm_coef_var = function(data, value, bw, type){
  
  
  result_left = one_glm_ft_coef_and_var(data=data, value=value, bw=bw, dire='left', type)
  result_right = one_glm_ft_coef_and_var(data=data, value=value, bw=bw, dire='right', type)
  
  if (type == 'poisson'){
    
    mu_left = exp(result_left[[1]][1])
    mu_right = exp(result_right[[1]][1])
    
    tau_hat = mu_right - mu_left
    
    sigma_left = result_left[[2]]
    sigma_right = result_right[[2]]
    
    asym_var = mu_left^2*sigma_left[1, 1] + mu_right^2*sigma_right[1, 1]
    
  }
  
  else if (type == 'binomial'){
    
    mu_left = exp(result_left[[1]][1])/(1+exp(result_left[[1]][1]))
    mu_right = exp(result_right[[1]][1])/(1+exp(result_right[[1]][1]))
    
    tau_hat = mu_right - mu_left
    
    sigma_left =result_left[[2]]
    sigma_right = result_right[[2]]
    
    asym_var = (exp(result_left[[1]][1])/(1+exp(result_left[[1]][1]))^2)^2 * sigma_left[1, 1] + (exp(result_right[[1]][1])/(1+exp(result_right[[1]][1]))^2)^2 * sigma_right[1, 1]
    
  }
  
  else if (type == 'quasipoisson'){
    
    mu_left = exp(result_left[[1]][1])
    mu_right = exp(result_right[[1]][1])
    
    tau_hat = mu_right - mu_left
    
    sigma_left = result_left[[2]]
    sigma_right = result_right[[2]]
    
    asym_var = mu_left^2*sigma_left[1, 1] + mu_right^2*sigma_right[1, 1]
    
  }
  
  
  return(c(tau_hat, asym_var))
  
}


### Zero infaltion case

#logistic, loglinear function

logistic = function(x) { exp(x)/(1+exp(x))}

loglinear = function(x) {exp(x)}



### zeroinfl-function

one_zeroinfl_ft_coef_and_var=function(data, value, bw, dire='left'){
  
  if (dire=='left'){
    result=data %>% filter(x<value & x>value-bw)
  }
  else if (dire=='right'){
    result=data %>% filter(x>=value & x<value+bw)
  }
  else {
    result='error : dire must be left or right'
  }
  
  K= dnorm(x=value, mean=result$x, sd=bw)
  data_revision = result %>% mutate(x=x-value)
  est=zeroinfl(y~x, data=data_revision)
  
  
  return(list(est$coefficients, est$vcov))
  
  
}


#logit : alpha, log : beta

# gradient of g(theta) - for delta method

gradient_ft = function(alpha0, beta0){
  
  result = c( exp(beta0)/(1+exp(alpha0)), 0, -exp(alpha0+beta0)/(1+exp(alpha0))^2, 0)
  
  return(as.matrix(result, nrow=4))
  
}



# compute tau_SRD and asymptotic variance by local likelihood


SRD_zeroinfl_coef_var=function(data, value, bw){
  
  
  result_left = one_zeroinfl_ft_coef_and_var(data=data, value=value, bw=bw, dire='left')
  result_right = one_zeroinfl_ft_coef_and_var(data=data, value=value, bw=bw, dire='right')
  
  pi_left = exp(result_left[[1]]$zero[1])/(1+exp(result_left[[1]]$zero[1]))
  lambda_left=exp(result_left[[1]]$count[1])
  
  pi_right = exp(result_right[[1]]$zero[1])/(1+exp(result_right[[1]]$zero[1]))
  lambda_right=exp(result_right[[1]]$count[1])
  
  
  sigma_left=result_left[[2]]
  sigma_right=result_right[[2]]
  grad_left=gradient_ft(alpha0=result_left[[1]]$zero[1], beta0=result_left[[1]]$count[1])
  grad_right=gradient_ft(alpha0=result_right[[1]]$zero[1], beta0=result_right[[1]]$count[1])
  
  
  
  nl=sum(data$x<value & data$x>value-bw)
  nr=sum(data$x>value & data$x<value+bw)
  
  tauhat = (1-pi_right)*lambda_right - (1-pi_left)*lambda_left
  
  asym_var=t(grad_left)%*%sigma_left%*%grad_left +t(grad_right)%*%sigma_right%*%grad_right
  
  
  return(c(tauhat, asym_var))
  
  
}


### Cross validation

### Cross-validation


### Finding quantile

quntilerange=function(data, q,cutoff){
  
  return(c(quantile(data$x[data$x<cutoff], q), quantile(data$x[data$x>=cutoff], 1-q)))
}




### Cross-validation

cv_rd_srd = function(data, cutoff, q, bw, type='ZIP'){
  
  
  cvrange=quntilerange(data=data, q=q, cutoff=cutoff)
  
  
  check_df=data[data$x>=cvrange[1] & data$x<=cvrange[2], ]
  check_x=check_df$x
  
  result=NULL
  
  
  if (type=='ZIP'){
    
    for (i in 1:length(check_x)){
      if (check_x[i]<cutoff){
        re = one_zeroinfl_ft_coef_and_var(data=data, value=check_x[i], bw=bw, dire='left')[[1]]
      }
      else if (check_x[i]>=cutoff){
        
        re = one_zeroinfl_ft_coef_and_var(data=data, value=check_x[i], bw=bw, dire='right')[[1]]
      }
      
      pi=exp(re$zero[1])/(1+exp(re$zero[1]))
      lambda=exp(re$count[1])
      
      check_y = check_df[i, 'y']
      
      like = ifelse(check_y==0, log(pi +(1-pi)*exp(-lambda)), log(1-pi) - lambda + check_y*log(lambda)-log(factorial(check_y)))
      
      result=c(result, like)
      
    }
  }
  
  
  else if (type=='linear'){
    
    for (i in 1:length(check_x)){
      
      if (check_x[i]<cutoff){
        
        re=one_linear_ft_coef_and_var(data=data, value=check_x[i], bw=bw, dire='left')[[1]]
      }
      else if (check_x[i]>=cutoff){
        
        
        re=one_linear_ft_coef_and_var(data=data, value=check_x[i], bw=bw, dire='right')[[1]]
      }
      
      mu=re[1]
      result=c(result, mu)
    }
    
    return(list(result, mean((check_df$y-result)^2)))
  }
  
  
  else {
    
    for (i in 1:length(check_x)){
      
      if (check_x[i]<cutoff){
        re = one_glm_ft_coef_and_var(data=data, value=check_x[i], bw=bw, dire='left', type=type)[[1]]
        
      }
      else{
        
        re = one_glm_ft_coef_and_var(data=data, value=check_x[i], bw=bw, dire='right', type=type)[[1]]
        
      }
      
      check_y = check_df[i, 'y']
      
      if (type =='binomial'){
        
        pi = exp(re[1])/(1+exp(re[1]))
        like = check_y * log(pi/(1-pi)) + log(1-pi)
        
      }
      else{
        lambda = exp(re[1])
        like = -lambda + check_y*log(lambda) - log(factorial(check_y))
      }
      
      result = c(result, like)
    }
    
    
  }
  
  return(list(result, sum(result)))
  
  
  
  
  
}


### Simulation function

simulation_ft = function(lambda, pi, cutoff, q, range, n, trial, type='ZIP', truediff){
  
  if (type=='ZIP'){
    
    
    
    # cross validation
    
    
    
    #tau estimation with replacement
    
    srd_bd_zip=NULL
    srd_bd_lm=NULL
    srd_est=NULL
    srd_var=NULL
    srd_lm_est=NULL
    srd_lm_var=NULL  
    srd_lm_est_IK_t=NULL
    srd_lm_var_IK_t=NULL  
    srd_lm_est_IK_r=NULL
    srd_lm_var_IK_r=NULL  
    
    for (ii in 1:trial){
      
      #data generateion
      
      
      x=runif(n, range[1], range[2])
      lamb=lambda(x)
      p=pi(x)
      
      z=rbernoulli(n, p=p)
      zpoi=rpois(n, lambda=lamb)
      
      y=(1-z)*zpoi
      data=data.frame(y=y, x=x)
      
      bandwidth=c(seq(0.1, 3, by=0.1))
      cv_result_lm=NULL
      cv_result_zip=NULL
      IKbw_r = IKbandwidth(data$x, data$y, cutpoint = 0, kernel="rectangular")
      IKbw_t = IKbandwidth(data$x, data$y, cutpoint = 0)
      
      # cross validation
      
      
      for (j in 1:length(bandwidth)){
        
        result_bd_zip=tryCatch(cv_rd_srd(data=data, cutoff=cutoff, q=q, bw=bandwidth[j], type='ZIP')[[2]], error=function(x) return(-Inf))
        result_bd_lm=tryCatch(cv_rd_srd(data=data, cutoff=cutoff, q=q, bw=bandwidth[j], type='linear')[[2]], error=function(x) return(Inf))
        if (is.nan(result_bd_zip)){
          result_bd_zip = -Inf
        }
        cv_result_zip=c(cv_result_zip, result_bd_zip)
        cv_result_lm=c(cv_result_lm, result_bd_lm)
        print('working')
        print(result_bd_zip)
        print(result_bd_lm)
        print('working')
        
        
        
      }
      
      print(paste0('Trial ', ii, ' is finished.'))
      
      bd_zip=bandwidth[which.max(cv_result_zip)]
      bd_lm=bandwidth[which.min(cv_result_lm)]
      
      
      
      
      result_zip=tryCatch(SRD_zeroinfl_coef_var(data=data, value=cutoff, bw=bd_zip), warning=function(x) return(c(Inf, Inf)), error=function(x) return(c(Inf, Inf)))
      
      result_lm=tryCatch(SRD_lm_coef_var(data, value=cutoff, bw=bd_lm), warning=function(x) return(c(Inf, Inf)), error=function(x) return(c(Inf, Inf)))
      result_lm_IK_t=tryCatch(SRD_lm_coef_var(data, value=cutoff, bw=IKbw_t), warning=function(x) return(c(Inf, Inf)), error=function(x) return(c(Inf, Inf)))
      result_lm_IK_r=tryCatch(SRD_lm_coef_var(data, value=cutoff, bw=IKbw_r), warning=function(x) return(c(Inf, Inf)), error=function(x) return(c(Inf, Inf)))
      
      srd_bd_zip[ii]=bd_zip
      srd_bd_lm[ii]=bd_lm
      srd_est[ii]=result_zip[1]
      srd_var[ii]=result_zip[2]
      srd_lm_est[ii]=result_lm[1]
      srd_lm_var[ii]=result_lm[2]
      srd_lm_est_IK_t[ii]=result_lm_IK_t[1]
      srd_lm_var_IK_t[ii]=result_lm_IK_t[2]
      srd_lm_est_IK_r[ii]=result_lm_IK_r[1]
      srd_lm_var_IK_r[ii]=result_lm_IK_r[2]
    }
    
    
    final_result=data.frame('bd_zip'=srd_bd_zip, 'bd_lm'=srd_bd_lm, 'zipest'=srd_est, 'zipvar'=srd_var, 'linest'=srd_lm_est, 'linvar'=srd_lm_var, 'linest_IK_t'=srd_lm_est_IK_t, 'linvar_IK_t'=srd_lm_var_IK_t, 'linest_IK_r'=srd_lm_est_IK_r, 'linvar_IK_r'=srd_lm_var_IK_r)# %>% filter(zipest!=Inf, linest!=Inf)
    final_result=final_result %>% mutate(zip_ci_left=zipest-qnorm(0.975, 0, 1)*sqrt(zipvar),
                                         zip_ci_right=zipest+qnorm(0.975, 0, 1)*sqrt(zipvar),
                                         lin_ci_left=linest-qnorm(0.975, 0, 1)*sqrt(linvar),
                                         lin_ci_right=linest+qnorm(0.975, 0, 1)*sqrt(linvar),
                                         lin_ci_left_IK_t=linest_IK_t-qnorm(0.975, 0, 1)*sqrt(linvar_IK_t),
                                         lin_ci_right_IK_t=linest_IK_t+qnorm(0.975, 0, 1)*sqrt(linvar_IK_t),
                                         lin_ci_left_IK_r=linest_IK_r-qnorm(0.975, 0, 1)*sqrt(linvar_IK_r),
                                         lin_ci_right_IK_r=linest_IK_r+qnorm(0.975, 0, 1)*sqrt(linvar_IK_r))
    
    final_result=final_result %>% mutate(zip_ci_in=(zip_ci_left < truediff & truediff<zip_ci_right),
                                         lin_ci_in=(lin_ci_left < truediff & truediff<lin_ci_right),
                                         lin_ci_in_IK_t=(lin_ci_left_IK_t < truediff & truediff<lin_ci_right_IK_t),
                                         lin_ci_in_IK_r=(lin_ci_left_IK_r < truediff & truediff<lin_ci_right_IK_r))
    
    
    
    
    
    
    
  }
  
  else if (type =='poisson' | type== 'quasipoisson'){
    
    
    # data generation
    
    x=runif(n, range[1], range[2])
    
    lamb=lambda(x)
    
    
    zpoi=rpois(n, lambda=lamb)
    
    y=zpoi
    
    data = data.frame(y=y, x=x)
    
    
    # cross validation
    
    
    #tau estimation with replacement
    
    srd_bd_zip=NULL
    srd_bd_lm=NULL
    srd_est=NULL
    srd_var=NULL
    srd_lm_est=NULL
    srd_lm_var=NULL  
    srd_lm_est_IK_t=NULL
    srd_lm_var_IK_t=NULL  
    srd_lm_est_IK_r=NULL
    srd_lm_var_IK_r=NULL  
    
    for (ii in 1:trial){
      
      
      x=runif(n, range[1], range[2])
      lamb=lambda(x)
      
      zpoi=rpois(n, lambda=lamb)
      
      y=zpoi
      data=data.frame(y=y, x=x)
      
      
      
      bandwidth=c(seq(0.1, 3, by=0.1))
      cv_result_lm=NULL
      cv_result_poi=NULL
      IKbw_r = IKbandwidth(data$x, data$y, cutpoint = 0, kernel="rectangular")
      IKbw_t = IKbandwidth(data$x, data$y, cutpoint = 0)
      
      for (j in 1:length(bandwidth)){
        
        result_bd_poi=tryCatch(cv_rd_srd(data=data, cutoff=cutoff, q=q, bw=bandwidth[j], type=type)[[2]], error=function(x) return(-Inf))
        result_bd_lm=tryCatch(cv_rd_srd(data=data, cutoff=cutoff, q=q, bw=bandwidth[j], type='linear')[[2]], error=function(x) return(Inf))
        if (is.nan(result_bd_poi)){
          result_bd_poi = -Inf
        }
        cv_result_poi=c(cv_result_poi, result_bd_poi)
        cv_result_lm=c(cv_result_lm, result_bd_lm)
        print('working')
        print(result_bd_poi)
        print(result_bd_lm)
        print('working')
        
      }
      
      print(paste0('Trial ', ii, ' is finished.'))
      
      bd_zip=bandwidth[which.max(cv_result_poi)]
      bd_lm=bandwidth[which.min(cv_result_lm)]
      
      
      
      result_zip=tryCatch(SRD_glm_coef_var(data=data, value=cutoff, bw=bd_zip, type='poisson'), warning=function(x) return(c(Inf, Inf)), error=function(x) return(c(Inf, Inf)))
      
      result_lm=tryCatch(SRD_lm_coef_var(data, value=cutoff, bw=bd_lm), warning=function(x) return(c(Inf, Inf)), error=function(x) return(c(Inf, Inf)))
      result_lm_IK_t=tryCatch(SRD_lm_coef_var(data, value=cutoff, bw=IKbw_t), warning=function(x) return(c(Inf, Inf)), error=function(x) return(c(Inf, Inf)))
      result_lm_IK_r=tryCatch(SRD_lm_coef_var(data, value=cutoff, bw=IKbw_r), warning=function(x) return(c(Inf, Inf)), error=function(x) return(c(Inf, Inf)))
      
      srd_bd_zip[ii]=bd_zip
      srd_bd_lm[ii]=bd_lm
      srd_est[ii]=result_zip[1]
      srd_var[ii]=result_zip[2]
      srd_lm_est[ii]=result_lm[1]
      srd_lm_var[ii]=result_lm[2]
      srd_lm_est_IK_t[ii]=result_lm_IK_t[1]
      srd_lm_var_IK_t[ii]=result_lm_IK_t[2]
      srd_lm_est_IK_r[ii]=result_lm_IK_r[1]
      srd_lm_var_IK_r[ii]=result_lm_IK_r[2]
      
      
      
    }
    
    
    final_result=data.frame('bd_zip'=srd_bd_zip, 'bd_lm'=srd_bd_lm, 'zipest'=srd_est, 'zipvar'=srd_var, 'linest'=srd_lm_est, 'linvar'=srd_lm_var, 'linest_IK_t'=srd_lm_est_IK_t, 'linvar_IK_t'=srd_lm_var_IK_t, 'linest_IK_r'=srd_lm_est_IK_r, 'linvar_IK_r'=srd_lm_var_IK_r)# %>% filter(zipest!=Inf, linest!=Inf)
    final_result=final_result %>% mutate(zip_ci_left=zipest-qnorm(0.975, 0, 1)*sqrt(zipvar),
                                         zip_ci_right=zipest+qnorm(0.975, 0, 1)*sqrt(zipvar),
                                         lin_ci_left=linest-qnorm(0.975, 0, 1)*sqrt(linvar),
                                         lin_ci_right=linest+qnorm(0.975, 0, 1)*sqrt(linvar),
                                         lin_ci_left_IK_t=linest_IK_t-qnorm(0.975, 0, 1)*sqrt(linvar_IK_t),
                                         lin_ci_right_IK_t=linest_IK_t+qnorm(0.975, 0, 1)*sqrt(linvar_IK_t),
                                         lin_ci_left_IK_r=linest_IK_r-qnorm(0.975, 0, 1)*sqrt(linvar_IK_r),
                                         lin_ci_right_IK_r=linest_IK_r+qnorm(0.975, 0, 1)*sqrt(linvar_IK_r))
    
    final_result=final_result %>% mutate(zip_ci_in=(zip_ci_left < truediff & truediff<zip_ci_right),
                                         lin_ci_in=(lin_ci_left < truediff & truediff<lin_ci_right),
                                         lin_ci_in_IK_t=(lin_ci_left_IK_t < truediff & truediff<lin_ci_right_IK_t),
                                         lin_ci_in_IK_r=(lin_ci_left_IK_r < truediff & truediff<lin_ci_right_IK_r))
    
    
    
    
    
    
  }
  
  
  else if(type=='binomial'){
    
    x=runif(n, range[1], range[2])
    
    prob = pi(x)
    
    z = rbernoulli(n, p=prob)
    data = data.frame(x=x, y=z)
    
    srd_bd_bin=NULL
    srd_bd_lm=NULL
    srd_bd_IK=NULL
    srd_est=NULL
    srd_var=NULL
    srd_est_IK=NULL
    srd_var_IK=NULL
    srd_lm_est=NULL
    srd_lm_var=NULL  
    srd_lm_est_IK_t=NULL
    srd_lm_var_IK_t=NULL  
    srd_lm_est_IK_r=NULL
    srd_lm_var_IK_r=NULL  
    
    for (ii in 1:trial){
      
      
      
      x=runif(n, range[1], range[2])
      prob = pi(x)
      
      z = rbernoulli(n, p=prob)
      data = data.frame(x=x, y=z)
      
      
      bandwidth=c(seq(0.1, 3, by=0.1))
      cv_result_lm=NULL
      cv_result_bin=NULL
      IKbw_r = IKbandwidth(data$x, data$y, cutpoint = 0, kernel="rectangular")
      IKbw_t = IKbandwidth(data$x, data$y, cutpoint = 0)
      
      DAT  =  data %>% arrange(x)
      DAT$y  =  ifelse(DAT$y==T, 0, 1)
      DAT$z  =  1- DAT$y
      J = length(DAT[1,])-2
      IKbw = tryCatch(rd_mnl_final(DAT,0,matrix(0,J,1),diag(J),matrix(0,J,1),0.90)$h_opt, error=function(x) return(Inf))
      
      for (j in 1:length(bandwidth)){
        
        result_bd_bin=tryCatch(cv_rd_srd(data=data, cutoff=cutoff, q=q, bw=bandwidth[j], type=type)[[2]], error=function(x) return(-Inf))
        result_bd_lm=tryCatch(cv_rd_srd(data=data, cutoff=cutoff, q=q, bw=bandwidth[j], type='linear')[[2]], error=function(x) return(Inf))
        if (is.nan(result_bd_bin)){
          result_bd_bin = -Inf
        }
        cv_result_bin=c(cv_result_bin, result_bd_bin)
        cv_result_lm=c(cv_result_lm, result_bd_lm)
        print('working')
        print(result_bd_bin)
        print(result_bd_lm)
        print('working')
        
        
      }
      print(paste0('Trial ', ii, ' is finished.'))
      bd_bin=bandwidth[which.max(cv_result_bin)]
      print(bd_bin)
      bd_lm=bandwidth[which.min(cv_result_lm)]
      print(bd_lm)
      print(IKbw)
      
      
      result_bin=tryCatch(SRD_glm_coef_var(data=data, value=cutoff, bw=bd_bin, type='binomial'), warning=function(x) return(c(Inf, Inf)), error=function(x) return(c(Inf, Inf)))
      result_bin_IK=tryCatch(SRD_glm_coef_var(data=data, value=cutoff, bw=IKbw, type='binomial'), warning=function(x) return(c(Inf, Inf)), error=function(x) return(c(Inf, Inf)))
      result_lm=tryCatch(SRD_lm_coef_var(data, value=cutoff, bw=bd_lm), warning=function(x) return(c(Inf, Inf)), error=function(x) return(c(Inf, Inf)))
      result_lm_IK_t=tryCatch(SRD_lm_coef_var(data, value=cutoff, bw=IKbw_t), warning=function(x) return(c(Inf, Inf)), error=function(x) return(c(Inf, Inf)))
      result_lm_IK_r=tryCatch(SRD_lm_coef_var(data, value=cutoff, bw=IKbw_r), warning=function(x) return(c(Inf, Inf)), error=function(x) return(c(Inf, Inf)))
      
      srd_bd_bin[ii]=bd_bin
      srd_bd_lm[ii]=bd_lm
      srd_bd_IK[ii]=IKbw
      srd_est[ii]=result_bin[1]
      srd_var[ii]=result_bin[2]
      srd_est_IK[ii]=result_bin_IK[1]
      srd_var_IK[ii]=result_bin_IK[2]
      srd_lm_est[ii]=result_lm[1]
      srd_lm_var[ii]=result_lm[2]
      srd_lm_est_IK_t[ii]=result_lm_IK_t[1]
      srd_lm_var_IK_t[ii]=result_lm_IK_t[2]
      srd_lm_est_IK_r[ii]=result_lm_IK_r[1]
      srd_lm_var_IK_r[ii]=result_lm_IK_r[2]
      
      
      
      
      
      
    }
    
    
    final_result=data.frame('bd_bin'=srd_bd_bin, 'bd_lm'=srd_bd_lm, 'bd_IK'=srd_bd_IK, 'zipest'=srd_est, 'zipvar'=srd_var, 'zipest_IK'=srd_est_IK, 'zipvar_IK'=srd_var_IK, 'linest'=srd_lm_est, 'linvar'=srd_lm_var, 'linest_IK_t'=srd_lm_est_IK_t, 'linvar_IK_t'=srd_lm_var_IK_t, 'linest_IK_r'=srd_lm_est_IK_r, 'linvar_IK_r'=srd_lm_var_IK_r)# %>% filter(zipest!=Inf, linest!=Inf)
    final_result=final_result %>% mutate(zip_ci_left=zipest-qnorm(0.975, 0, 1)*sqrt(zipvar),
                                         zip_ci_right=zipest+qnorm(0.975, 0, 1)*sqrt(zipvar),
                                         zip_ci_left_IK=zipest-qnorm(0.975, 0, 1)*sqrt(zipvar_IK),
                                         zip_ci_right_IK=zipest+qnorm(0.975, 0, 1)*sqrt(zipvar_IK),
                                         lin_ci_left=linest-qnorm(0.975, 0, 1)*sqrt(linvar),
                                         lin_ci_right=linest+qnorm(0.975, 0, 1)*sqrt(linvar),
                                         lin_ci_left_IK_t=linest_IK_t-qnorm(0.975, 0, 1)*sqrt(linvar_IK_t),
                                         lin_ci_right_IK_t=linest_IK_t+qnorm(0.975, 0, 1)*sqrt(linvar_IK_t),
                                         lin_ci_left_IK_r=linest_IK_r-qnorm(0.975, 0, 1)*sqrt(linvar_IK_r),
                                         lin_ci_right_IK_r=linest_IK_r+qnorm(0.975, 0, 1)*sqrt(linvar_IK_r))
    
    final_result=final_result %>% mutate(zip_ci_in=(zip_ci_left < truediff & truediff<zip_ci_right),
                                         zip_ci_in_IK=(zip_ci_left_IK < truediff & truediff<zip_ci_right_IK),
                                         lin_ci_in=(lin_ci_left < truediff & truediff<lin_ci_right),
                                         lin_ci_in_IK_t=(lin_ci_left_IK_t < truediff & truediff<lin_ci_right_IK_t),
                                         lin_ci_in_IK_r=(lin_ci_left_IK_r < truediff & truediff<lin_ci_right_IK_r))
    
    
    
    
    
    
    
    
  }
  
  
  
  return(final_result)
  
}


make_csv = function(lambda, pi, cutoff, q, range, n_case, trial, type='ZIP', truediff){
  
  
  name = paste0(n_case, '_', 'type','_', type, '_', 'trial','_',trial, '_NULL.csv')
  
  result = simulation_ft(lambda, pi, cutoff, q, range, n_case, trial, type, truediff)
  
  write.csv(result, name)
  
}


lambdaft=function(x) {return(loglinear(-0.4*x+1.1 + (x>=0)*0+0.2*x*(x>=0)))}
pift = function(x) {return(logistic(0.3*x+1.5-(x>=0)*0))}


truediff_zip = exp(1.1)/(1+exp(1.5)) - exp(1.1)/(1+exp(1.5))
truediff_binomial =  exp(1.5)/(1+exp(1.5)) - exp(1.5)/(1+exp(1.5))
truediff_poi = exp(1.1) - exp(1.1)

# Vary n_case and type: (100,300,500,1000,3000); ("ZIP", "binomial", "poisson")
make_csv(lambda = lambdaft, pi = pift, cutoff =0, q=0.5, range=c(-3, 3), n_case=1000, trial = 1, type='ZIP', truediff = truediff_zip)

