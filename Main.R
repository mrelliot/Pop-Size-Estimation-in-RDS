library(Rcpp)
library(foreach)
library(partitions)#restrictedparts
library(RcppAlgos)#permuteGeneral
library(doParallel)
library(purrr)
library(dplyr)
library(truncnorm)
library(gtools)
library(iterpc)
sourceCpp("calProbDup.cpp")

########################### calProb_ss #############################
calProb_ss<-function(sample,pk_dist,N_list){
  # K<-nrow(pk_dist)
  d_ave<-sum(pk_dist$d_act*pk_dist$pk)
  uni_loc<-which(!duplicated(sample$IDrecruited))
  logp_list<-sapply(N_list,function(N)
    sum(log(sample$Drecruited[uni_loc]))-sum(log(N*d_ave-c(0,cumsum(sample$Drecruited[uni_loc][1:(length(uni_loc)-1)])))))
  return(logp_list)
}



calProb_obs<-function(sample,pk_dist,N_list){
  
  ##### approximated denominator ############
  num_seed<-sum(sample$wave==1)
  tb_wave2<-sample%>%filter(wave==2) %>%group_by(IDrecruiter)%>%summarise(n=n())
  tb_afterwave2<-sample%>%filter(wave>2) %>%group_by(IDrecruiter)%>%summarise(n=n())
  num_subtract<-data.frame(subtract=0:max(num_seed,tb_wave2$n,tb_afterwave2$n+1),count=0)
  ## wave 1 (seed) 
  num_subtract[num_subtract$subtract%in%(0:(num_seed-1)),"count"]<-1 
  ## wave 2
  tb<-as.data.frame(table(Reduce("c", sapply(tb_wave2$n,function(x) rep(1:x)))))
  num_subtract[num_subtract$subtract%in%tb$Var1,"count"]<-num_subtract[num_subtract$subtract%in%tb$Var1,"count"]+tb$Freq
  
  ## after wave 2
  tb<-as.data.frame(table(Reduce("c", sapply(tb_afterwave2$n,function(x) rep(2:(x+1))))))
  num_subtract[num_subtract$subtract%in%tb$Var1,"count"]<-num_subtract[num_subtract$subtract%in%tb$Var1,"count"]+tb$Freq
  ########################### pass parameters #################
  ######## aggregate pk #######################
  pk_dist_new<-pk_dist[,c("d_act","pk")]
  i=1
  while(i<nrow(pk_dist_new)){
    accumulate<-0
    while(pk_dist_new$pk[i]<.25-accumulate|abs(pk_dist_new$pk[i]-(.25-accumulate))<1e-10){
      accumulate<-accumulate+pk_dist_new$pk[i]
      i<-i+1
      if(i>nrow(pk_dist_new))
        break
    }
    if(accumulate<0.25&abs(accumulate-.25)>1e-10){
      if(i<nrow(pk_dist_new)){
        pk_dist_new[(i+2):(nrow(pk_dist_new)+1),]<-pk_dist_new[(i+1):nrow(pk_dist_new),]
        pk_dist_new[i+1,]<-c(pk_dist_new[i,"d_act"], pk_dist_new[i,"pk"]-(.25-accumulate))
        pk_dist_new[i,"pk"]<-.25-accumulate
        i<-i+1
      }else{
        pk_dist_new[(nrow(pk_dist_new)+1),]<-c(pk_dist_new[i,"d_act"], pk_dist_new[i,"pk"]-(.25-accumulate))
        pk_dist_new[i,"pk"]<-.25-accumulate
        i<-i+1
      }
      
      
    }
  }
  newPk<-c()
  t<-lapply(cumsum(pk_dist_new$pk),function(x)which(abs(x-seq(0,1,0.25))<1e-10))
  loc<-c(0,which(Reduce("c",lapply(t,length))==1))
  for(i in 2:length(loc)){
    newPk<-c(newPk, sum(pk_dist_new$d_act[(loc[i-1]+1):loc[i]]*pk_dist_new$pk[(loc[i-1]+1):loc[i]]/sum(pk_dist_new$pk[(loc[i-1]+1):loc[i]])))
  }
  pk_dist_agg<-data.frame(d_act=newPk,pk=1/length(newPk))
  pk_dist_agg_mtx<-as.matrix(pk_dist_agg[,c("d_act","pk")])
  colnames(pk_dist_agg_mtx)<-NULL
  rownames(pk_dist_agg_mtx)<-NULL
  
  num_uni<-length(unique(sample$IDrecruited))
  K<-nrow(pk_dist_agg)
  decompose_num<-restrictedparts(num_uni,m=K)
  decompose_mtx<-as.matrix(decompose_num)
  colnames(decompose_mtx)<-NULL
  rownames(decompose_mtx)<-NULL
  
  num_dup<-nrow(sample)-num_uni
  if(num_dup>0){
    decompose_num2<-restrictedparts(num_dup,m=K)
    decompose_mtx2<-as.matrix(decompose_num2)
    colnames(decompose_mtx2)<-NULL
    rownames(decompose_mtx2)<-NULL
  }
  else
    decompose_mtx2<-matrix(NA,1,1)
  
  num_subtract_mtx<-as.matrix(num_subtract)
  colnames(num_subtract_mtx)<-NULL
  rownames(num_subtract_mtx)<-NULL
  
  dup<-calProb_dup(pk_dist_agg_mtx,N_list,decompose_mtx,decompose_mtx2,recruits_num=nrow(sample),num_subtract_mtx)
  log_ss<-calProb_ss(sample,pk_dist,N_list)
  d_obs_unique<-sample$Drecruited[!duplicated(sample$IDrecruited)]
  nk_list<-as.vector(table(d_obs_unique))
  
  v1<-t(outer(N_list,pk_dist$pk,"*"))
  v2<-v1-nk_list+1
  v1[v1<v2]<-0
  v2[v1<v2]<-0
  log_coeff_list<-sapply(1:ncol(v1),function(col)mapply(function(x,y) sum(log(x:y)), v1[,col],v2[,col]))
  log_coeff_list[is.infinite(log_coeff_list)]<-0
  log_coeff<-colSums(log_coeff_list)
  
  out<-dup*exp(log_ss+log_coeff)/N_list
  return(log(out))
} 



calProb_pk<-function(N,pk_dist,pi_dist,nk_list,d_obs_unique){
  #### p(pk_list|d_obs,dupVec,N)
  Nkleft_list<-N*pk_dist$pk-nk_list
  D<-N*sum(pk_dist$pk*pk_dist$d_act)
  if(any(D-sum(d_obs_unique)-pk_dist$d_act<0))
    return(NA)
  log_p<-lgamma(N-length(d_obs_unique)+1)-sum(lgamma(Nkleft_list+1))+sum(Nkleft_list*log(pi_dist$pi))+#log(multichoose(Nkleft_list, bigz = T))+
    sum(Nkleft_list*log(D-sum(d_obs_unique)-pk_dist$d_act))-sum(Nkleft_list*log(D-sum(d_obs_unique)))+
    sum(log(d_obs_unique))-sum(log(D-c(0,cumsum(d_obs_unique)[-length(d_obs_unique)])))
  
  
  return(log_p)
}

MCMC <-function(niter,Nmin,Nmax,sample){
  degrees<-unique(sample$Drecruited)
  degrees<-degrees[order(degrees)]
  K<-length(degrees)
  
  d_obs_unique<-sample$Drecruited[!duplicated(sample$IDrecruited)]
  nk_list<-as.vector(table(d_obs_unique))
  results<- data.frame(matrix(ncol=2*K+8, nrow = 0))
  colnames(results)<-c("iter","N","log_P(N|pk)","naccept_N","df_N",
                       paste0("Pi",1:K),paste0("P",1:K),"log_P(pk|N)","naccept_pk","df_pk")
  ### initialize N and pi, pk
  N<-Nmax/2
  pi_list<-nk_list/degrees/(sum(nk_list/degrees))
  pi_dist<-data.frame(d_act=degrees,pi=pi_list)
  pk_dist<-data.frame(d_act=degrees,pk=c(rmultinom(1,N,pi_dist$pi))/N)
  while(any(N*pk_dist$pk<nk_list)){
    # print("Wrong initialization")
    pk_dist<-data.frame(d_act=degrees,pk=c(rmultinom(1,N,pi_dist$pi))/N)
  }
  
  prob_N<-calProb_obs(sample,pk_dist,N)
  prob_pk<-calProb_pk(N,pk_dist,pi_dist,nk_list,d_obs_unique)
  
  df_N<-100
  df_pk<-0.01
  every200<-1
  results[1,]<-results[2,]<- results[3,]<-rep(NA,ncol(results))
  results[4,]<-c(1,N,prob_N,1,df_N,pi_dist$pi,pk_dist$pk,prob_pk,1,df_pk)
  
  iter<-2
  while(iter <=niter){
    N_cur <- as.numeric(results[4*(iter-1),"N"])
    pk_cur <- as.numeric(results[4*(iter-1),paste0("P",1:K)])
    pk_dist_cur<-data.frame(d_act=degrees,pk=pk_cur)
    pi_cur <- as.numeric(results[4*(iter-1),paste0("Pi",1:K)])
    pi_dist_cur <- data.frame(d_act=degrees,pi=pi_cur)
    
    ## generate proposed new pi,pk
    pi_proposed<-AcceptPi<- c(rdirichlet(1,N_cur*pk_cur+1))
    pi_dist_proposed<-data.frame(d_act=degrees,pi=pi_proposed)
    results[4*(iter-1)+1,]<-c(iter,N_cur,NA,NA,NA,
                              pi_proposed,pk_cur,calProb_pk(N_cur,pk_dist_cur,pi_dist_proposed,nk_list,d_obs_unique),NA,NA)

    pk_proposed<-c(rmultinom(1,N_cur,pi_proposed))/N_cur#+epsilon
    pk_dist_proposed<-data.frame(d_act=degrees,pk=pk_proposed)
    ## generate proposed new N
    N_proposed<-round(rtruncnorm(1, a=Nmin, b=Nmax, mean =N_cur, sd = df_N))
    
    ### draw pk
    ## if proposed pk out of range, reject new pk
    if(any(pi_proposed<0)|any(N_cur*pk_proposed<nk_list)){
      AcceptP<-pk_cur
      naccep_pk<-as.numeric(results[4*(iter-1),"naccept_pk"])
      prob_pk<-0
      prob_N<-results$`log_P(N|pk)`[4*(iter-1)]
      prob_Nproposed<-NA#calProb_obs(sample,pk_dist_cur,N_proposed)
      prob_pk_final_temp<-results[4*(iter-1)+1,"log_P(pk|N)"]
    }else{
      prob_pk<-calProb_pk(N_cur,pk_dist_proposed,pi_dist_proposed,nk_list,d_obs_unique)
      A<-min(1,exp(prob_pk-as.numeric(results[4*(iter-1)+1,"log_P(pk|N)"])+
                     dmultinom(round(pk_cur*N_cur),sum(round(pk_cur*N_cur)),pi_proposed,log=T)-#+epsilon
                     dmultinom(pk_proposed*N_cur,N_cur,pi_proposed,log=T)))#+epsilon
      if(runif(1)<A){# accept move with probabily min(1,A)
        AcceptP<-pk_proposed
        naccep_pk<-as.numeric(results[4*(iter-1),"naccept_pk"])+1
        prob_Ns<-calProb_obs(sample,pk_dist_proposed,c(N_cur,N_proposed))
        prob_N<-prob_Ns[1]
        prob_Nproposed<-prob_Ns[2]
        prob_pk_final_temp<-prob_pk
      }else { # otherwise "reject" move, and stay where we are
        
        AcceptP<-pk_cur
        naccep_pk<-as.numeric(results[4*(iter-1),"naccept_pk"])
        prob_N<-results$`log_P(N|pk)`[4*(iter-1)]
        prob_Nproposed<-NA#calProb_obs(sample,pk_dist_cur,N_proposed)
        prob_pk_final_temp<-results[4*(iter-1)+1,"log_P(pk|N)"]
        
      }
    }
    pk_dist_accept<-data.frame(d_act=degrees,pk=AcceptP)
    pi_dist_accept<-data.frame(d_act=degrees,pi=AcceptPi)
    results[4*(iter-1)+2,]<-c(iter,N_cur,prob_N,NA,NA,
                              pi_proposed,pk_proposed, prob_pk,naccep_pk,df_pk)
    
    
    ### draw N
    
    if(any(N_proposed*AcceptP<nk_list)){
      AcceptN<-N_cur
      naccept_N<-as.numeric(results[4*(iter-1),"naccept_N"])
      prob_Nproposed<-0
      prob_N_final<-results[4*(iter-1)+2,"log_P(N|pk)"]
      prob_pk_final<-prob_pk_final_temp
    }else{
      if(is.na(prob_Nproposed))
        prob_Nproposed<-calProb_obs(sample,pk_dist_accept,N_proposed)
      A<-min(1,exp(prob_Nproposed-as.numeric(results[4*(iter-1)+2,"log_P(N|pk)"]))*
               dtruncnorm(N_cur,a=Nmin,b=Nmax,mean=N_proposed,sd=df_N)/dtruncnorm(N_proposed,a=Nmin,b=Nmax,mean=N_cur,sd=df_N))
      if(runif(1)<A){# accept move with probabily min(1,A)
        AcceptN<-N_proposed
        naccept_N<-as.numeric(results[4*(iter-1),"naccept_N"])+1
        prob_N_final<-prob_Nproposed
        prob_pk_final<-calProb_pk(AcceptN,pk_dist_accept,pi_dist_accept,nk_list,d_obs_unique)
      }else { # otherwise "reject" move, and stay where we are
        AcceptN<-N_cur
        naccept_N<-as.numeric(results[4*(iter-1),"naccept_N"])
        prob_N_final<-results[4*(iter-1)+2,"log_P(N|pk)"]
        prob_pk_final<-prob_pk_final_temp
      }
    }
    results[4*(iter-1)+3,]<-c(iter,N_proposed,prob_Nproposed,naccept_N,df_N,
                              rep(NA,2*K),NA,NA,NA)
    
    
    results[4*(iter-1)+4,]<-c(iter,AcceptN,prob_N_final,naccept_N,df_N,
                              AcceptPi,AcceptP,prob_pk_final,naccep_pk,df_pk)
    
    iter<-iter+1
    if(ceiling(iter/200)>every200){
      every200<-every200+1
      Naccept_rate<-(as.numeric(results$naccept_N[4*(iter-1)])-as.numeric(results$naccept_N[4*(iter-200)]))/200
      if(Naccept_rate<0.5&Naccept_rate>.4)
        df_N<-df_N
      if(Naccept_rate>=0.5)
        df_N<-df_N*2
      if(Naccept_rate<=0.4)
        df_N<-df_N/2  
      Paccept_rate<-(as.numeric(results$naccept_pk[4*(iter-1)])-as.numeric(results$naccept_pk[4*(iter-200)]))/200#as.numeric(results$naccept_pk[3*(iter-1)])/(iter-1)
      if(Paccept_rate<0.3&Paccept_rate>.15)
        df_pk<-df_pk
      if(Paccept_rate>=0.3)
        df_pk<-min(1/K,df_pk+0.005)
      if(Paccept_rate<=0.15)
        df_pk<-max(0.001,df_pk-0.005)
      
    }
  }
  return(results)
}


out<-MCMC(niter,Nmin,Nmax,sample)