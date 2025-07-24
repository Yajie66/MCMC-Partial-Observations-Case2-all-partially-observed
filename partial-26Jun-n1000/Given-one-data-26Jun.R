
ssa.f.total<-function(theta,S0){
  
  Sini<-S0
  
  names(Sini)<-c('S1','S2')
  
  S.can<-S0
  
  t<-0
  
  Reaction<-matrix(c(-1,0,
                     0,+1,
                     0,-1),byrow=TRUE,ncol=2)  
  #Jump direction of (mRNA, protein) 
  
  
  # while(t[length(t)] <= Mt){
  while( !all(S.can == c(0, 0)) ){
    
    
    theta.can<-c(theta[1]*S.can[1], theta[2]*S.can[1], theta[3]*S.can[2]) 
    #calculate the current intensity function \lambda(x)=theta*x
    
    # Check if theta.can is c(0,0,0) and break the loop
    if (all(theta.can == c(0, 0, 0))) {
      break  # Exit the while loop
    }
    
    lambda<-sum(theta.can, na.rm = T) 
    
    wait<-rexp(1,rate=lambda) #waiting is exp(sum_intensities)
    
    
    t[length(t)+1] <-t[length(t)]+wait
    
    
    prob<-theta.can/lambda   
    
    R.can<-sample(size =1,x=c(1,2,3), prob = prob) #next jump direction
    
    
    
    if(R.can==1){
      S.can<-S.can+Reaction[1,] 
    }else{
      
      if(R.can==3){
        S.can<-S.can+Reaction[3,]
      }else{
        
        S.can<-S.can+Reaction[2,]
      }
    }
    
    
    Sini<-rbind.data.frame(Sini,S.can)
    
  }
  
  
  names(Sini)<-paste0('Species',c('mRNA','Protein'))
  
  jump.time.species<-cbind.data.frame(jump.times=t,Sini)
  return(jump.time.species)
}


true.theta=  c(0.5,80,1)
N=1000
theta1=true.theta[1]     #true parameter
 theta2.o= true.theta[2]
 theta3= true.theta[3]
 
 S0 <- c(N, 0)
 theta.ori<-c(theta1,theta2.o/N,theta3)
 
 one.data.out<- ssa.f.total(theta=theta.ori, S0=S0)
 
 
 
 file_name <- paste0("one_given_data_N=",N, "_",
                     format(Sys.Date(), "%Y%m%d"), ".csv")
 file_name
 write.csv(one.data.out, file = file_name, row.names = FALSE)
 cat("CSV saved to", normalizePath( file_name ), "\n") 
 
 
 
 
 