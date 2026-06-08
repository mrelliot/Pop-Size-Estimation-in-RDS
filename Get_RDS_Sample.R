library(dplyr)
############### network-construction ######################
data_prep<-function(trueN,maxK,minK,interval){
  try<-TRUE
  while(try){
    deg<-data.frame(id=1:trueN,d_exp=0,d_act=0)
    # assume degree ~ uniform distribution
    set.seed(123)
    deg$d_exp<-sample(seq(minK,maxK,interval),size=trueN,replace=T)
    net<-matrix(0,nrow=trueN,ncol=trueN)
    for(i in 1:(trueN-1)){
      d_left<-deg$d_exp[i]-sum(net[1:i,i])
      if(d_left<1)
        next
      friend<-(i+1):trueN
      friend<-friend[deg$d_act[friend]<deg$d_exp[friend]]
      link<-sample(friend,min(d_left,trueN-i,length(friend)),replace = F)
      net[i,link]<-1
      deg$d_act[link]<-deg$d_act[link]+1
    }
    net[lower.tri(net)] <-t(net)[lower.tri(net)]
    deg$d_act<-apply(net,1,sum)
    uni<-unique(deg$d_act)
    uni<-uni[order(uni)]
    if(all(diff(uni)==interval)){
      try<-FALSE
    }
      
  }
   pk_dist<-deg%>%
    group_by(d_act)%>%
    summarise(n=n())%>%
    mutate(pk=n/sum(n))%>%
     select(-n)
  assign("net",net,envir=.GlobalEnv)
  #write.csv(net,"net.csv",row.names = F)
  assign("deg",deg,envir=.GlobalEnv)
  saveRDS(deg,"deg.RDS")
  write.table(deg,"deg.txt",sep='\t',row.names = FALSE,col.names = TRUE)
  assign("pk_dist",pk_dist,envir=.GlobalEnv)
  saveRDS(pk_dist,"pk_dist.RDS")
  write.table(pk_dist,"pk_dist.txt",sep='\t',row.names = FALSE,col.names = TRUE)
}

########### recruitment process ##########################
recruiting<-function(net,deg,num_seed,num_coupon,require_size){

  out<-data.frame(chain=integer(),wave=integer(),dup=logical(),IDrecruited=integer(),IDrecruiter=integer())
  ## select seeds
  for(i in 1:num_seed){
    seed<-sample(setdiff(deg$id,out$IDrecruited),size=1,prob=deg$d_act[setdiff(deg$id,out$IDrecruited)])
    out[nrow(out)+1,]<-c(i,1,FALSE,seed,-1)
  }

  alive_chain<-1:num_seed
  current_wave<-1
  ## recruit subsequent units
  while(length(unique(out$IDrecruited))<require_size&length(alive_chain)>0){
    for(i in alive_chain){
      if(length(unique(out$IDrecruited))>=require_size)
        break
      TBrecruiter<-out[out$chain==i&out$wave==current_wave&out$dup==FALSE,]
      for(j in 1:nrow(TBrecruiter)){
        recruiter<-TBrecruiter$IDrecruited[j]
        link<-deg$id[which(net[recruiter,]==1)]
        link<-setdiff(link,TBrecruiter$IDrecruiter[j])
        if(length(link)<1){
          next
        }
        else{
          recruits<-sample(link,size=min(length(link),sample(1:num_coupon,1)),prob = deg$d_act[link])
          if(length(recruits)>=(require_size-length(unique(out$IDrecruited)))){
            recruits<-recruits[1:(require_size-length(unique(out$IDrecruited)))]
            checkDup<-recruits%in%out$IDrecruited
            out<-rbind(out,data.frame(chain=i,wave=current_wave+1,dup=checkDup,
                                      IDrecruited=recruits,IDrecruiter=recruiter))
            break
          }
          else{
            checkDup<-recruits%in%out$IDrecruited
            out<-rbind(out,data.frame(chain=i,wave=current_wave+1,dup=checkDup,
                                      IDrecruited=recruits,IDrecruiter=recruiter))
          }
         
          
        }
      }
      alive_recruits<-sum(out$dup[out$chain==i&out$wave==current_wave+1]==FALSE)
      if(alive_recruits==0)
        alive_chain<-setdiff(alive_chain,i)
      
      
      
     
      
    }
    current_wave<-current_wave+1
  }
  out$Drecruited<-deg$d_act[out$IDrecruited]
  out$Drecruiter<-apply(out,1, function(x) ifelse(x["IDrecruiter"]!=-1,deg$d_act[x["IDrecruiter"]],NA))
  return(out)
  
}

#########################################################
get_samples<-function(sim_size,net,deg,num_seed,num_coupon,require_size){
  sim<-1
  SMRY<-as.data.frame(matrix(rep(NA,5),nrow=1))
  names(SMRY)<-c("sampleID","num_dup","uni_size","uni_deg","max_wave")
  while(sim<=sim_size){
    set.seed(sim)
    out<-recruiting(net,deg,num_seed,num_coupon,require_size)
    while(length(unique(out$IDrecruited))<require_size){
      out<-recruiting(net,deg,num_seed,num_coupon,require_size)
    }
    saveRDS(out,paste0("sample_",sim,".RDS"))
    num_dup<-sum(duplicated(out$IDrecruited))
    uni_size<-nrow(out)-num_dup
    uni_deg<-length(unique(out$Drecruited))
    max_wave<-max(out$wave)
    SMRY[nrow(SMRY)+1,]<-c(sim,num_dup,uni_size,uni_deg,max_wave)
    sim<-sim+1
  }
  SMRY<-SMRY[-1,]
  rownames(SMRY)<-1:sim_size
  saveRDS(SMRY,"SMRY.RDS")
  # write.table(SMRY,"SMRY.txt",sep='\t',row.names = FALSE,col.names = TRUE)
  
}



toge<-function(trueN,pct,num_seed,num_coupon,maxK,minK,interval,sim_size){
  ### generate population
  data_prep(trueN,maxK,minK,interval)
  ### generate samples
  get_samples(sim_size,net,deg,num_seed,num_coupon,pct*trueN)
  
}

toge(trueN,pct,num_seed,num_coupon,maxK,minK,interval,sim_size)

