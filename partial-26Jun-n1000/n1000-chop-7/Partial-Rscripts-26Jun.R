

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



step.fn<- function(all.data,chop.number.times){
  time.r<- all.data$jump.times
  
  cond <- all( all.data[(nrow(all.data)-1),2:3] == c(0,0) )   #if the species all die out before T, the reaction final is die out time
  
  if( cond ){
    reaction.final.time <-  time.r[length(time.r)-1]
  }else{  reaction.final.time<- time.r[length(time.r)]  }  #if the species not die out at time T, 
  
  
  chop.numebr= 2^{chop.number.times}
  
  each.inter.len <- reaction.final.time/(chop.numebr+1)   #interval number
  
  partial.equal.obse.time<- rep(NA,l=(chop.numebr+2))
  
  for(i in  1:chop.numebr){
    partial.equal.obse.time[i+1]<- each.inter.len * i
  }
  partial.equal.obse.time[1]<-time.r[1]
  partial.equal.obse.time[length(partial.equal.obse.time)]<-reaction.final.time
  
  return(list(partial.equal.obse.time=partial.equal.obse.time,chop.numebr=chop.numebr ))
}


chop.time.equally.fn<- function( all.data, chop.number.times,  know.reaction.final.index, T.protein.unknown ){
  #if we known the whole reaction die out time, the [0,T] is the given by the last line of data
  if(know.reaction.final.index == TRUE){
    
    T.protein = all.data$jump.times[nrow(all.data)]
    with.T.all.data<- all.data
    
  }else{
    T.protein = T.protein.unknown
    #else we chop the time respect to some random observation time
    
    time.r<- all.data$jump.times
    if(time.r[length(time.r)]<= T.protein){
      index<- length(time.r)+1 #since the index is +1
      #if all false, the which.max will return 1
    }else{ index<- which.max( time.r >= T.protein)}
    
    
    all.data  <- all.data[1:(index-1),]
    last.interval<- as.data.frame( c(T.protein,  all.data [nrow( all.data ),-1] )  )
    colnames(last.interval)<- colnames(all.data)
    
    with.T.all.data <- rbind(all.data,last.interval)
    
  }
  ###################################################################################################
  #give them jump events index
  with.T.all.data$SpeciesProtein<- as.numeric(with.T.all.data$SpeciesProtein )
  with.T.all.data$SpeciesmRNA <- as.numeric(with.T.all.data$SpeciesmRNA)
  with.T.all.data$jump.type<- rep(NA, l= nrow(with.T.all.data))
  for(j in 2:nrow(with.T.all.data) ){
    
    if ( with.T.all.data$SpeciesProtein[j] - with.T.all.data$SpeciesProtein[j-1] == +1 ){
      with.T.all.data$jump.type[j] <- c("birth.event")
      
    }
    if ( with.T.all.data$SpeciesmRNA[j] - with.T.all.data$SpeciesmRNA[j-1] == -1 ){
      with.T.all.data$jump.type[j]<- c("mrna.death.event")}
    if ( with.T.all.data$SpeciesProtein[j] - with.T.all.data$SpeciesProtein[j-1] == -1 ){
      with.T.all.data$jump.type[j]<- c("death.event")}
  }
  
  
  with.T.all.data[,1:3]  <- data.frame(
    lapply(with.T.all.data[,1:3], function(col) as.numeric(col)),
    check.names = FALSE
  )
  ###################################################################################################    
  #partial data
  all.out<- step.fn(all.data = with.T.all.data, chop.number.times = chop.number.times)
  
  partial.obse.times<- all.out$partial.equal.obse.time
  
  chop.number<- all.out$chop.numebr
  
  partial.all.species.data<- as.data.frame(matrix(NA, nrow= length(partial.obse.times), ncol=4))
  colnames( partial.all.species.data)<- c("partial.observed.time", "obs.SpeciesmRNA" ,   "obs.SpeciesProtein","jump.type")
  
  
  for(j in 1:length(partial.obse.times)){
    
    index <- findInterval(partial.obse.times[j], with.T.all.data$jump.times)
    one.row <-  cbind.data.frame("partial.observed.time"=partial.obse.times[j], 
                                 "obs.SpeciesmRNA"= as.numeric( with.T.all.data[index,2]), 
                                 "obs.SpeciesProtein"= as.numeric( with.T.all.data[index,3]),
                                 "jump.type"="Partial.obse.times", deparse.level =0)
    partial.all.species.data[j,]<-  one.row 
    
    
  }
  
  
  return(list(partial.all.species.data=partial.all.species.data, 
              with.T.all.data= with.T.all.data, 
              T.protein= T.protein,  
              chop.number=  chop.number) )
}

cond.para.based.full.data.fn<- function(all.data.T){
  
  #extract the components in ith interval 
  time.r <-   as.numeric( all.data.T$jump.times  )  #include the head and tail of the partial observed data
  protein.full <-  as.numeric(  all.data.T$SpeciesProtein  )
  mrna.full <- as.numeric( all.data.T$SpeciesmRNA  )
  #extract the na,nb,nd and integral of X1,X2  
  l.t<- length(time.r) #total is with the head and tail T
  int.X2<- 0
  int.X1<-0
  nd =0  #protein death events number
  na =  mrna.full[1] - mrna.full[length(mrna.full)]  #mRNA degradation number
  nb=0
  for(i in 2:l.t){
    
    int.X2<- int.X2+ (time.r[i]-time.r[i-1])*protein.full[i-1]
    
    int.X1<- int.X1 +(time.r[i]-time.r[i-1])*mrna.full[i-1]
    
    if(protein.full[i]-protein.full[i-1] == -1){  nd = nd +1}
    if( protein.full[i]-protein.full[i-1] == +1){nb=nb+1}
    
  }
  
  #nb= (l.t-2)- na - nd #protein birth events number
  
  return(list(int.X1=int.X1,int.X2=int.X2, na=na,nb=nb,nd=nd))     
}  

full.cond.theta.gamma.fn<- function(full.data.summary.out.T, N,alp1, alp2, bet1, bet2,alp3, bet3 ){
  
  na= full.data.summary.out.T$na
  nb= full.data.summary.out.T$nb
  nd= full.data.summary.out.T$nd
  int.X1= full.data.summary.out.T$int.X1
  int.X2= full.data.summary.out.T$int.X2
  
  #return scaled version theta2
  a1<- na + alp1
  b1<-int.X1  +bet1
  
  the1.can<-rgamma(1, shape=a1,rate=b1)
  
  a2<- nb +alp2
  b2<-  int.X1/N +   bet2 
  #here the theta2 is the the scaled version, like 10/n, our theta2.out is 10
  the2.can<-rgamma(1,shape=a2,rate=b2)
  
  a3=nd  +alp3
  b3=int.X2  +bet3 
  the3.out <- rgamma(n=1,shape=a3 , rate=b3 )
  
  theta.all.can<- c( the1.can, the2.can,the3.out )
  return(theta.all.can)  
}

ith.interval.mrna.samples.fixed.jumps.fn<- function(ith.all.data.T, theta,N){
  theta1<- theta[1]
  theta2.o<- theta[2]
  a<- theta1+ theta2.o/N  #this is the parameters needed to be updated afterwards
  
  protein.born.index<- which( ith.all.data.T$jump.type == "birth.event" )
  protein.birth.times<- as.numeric( ith.all.data.T$jump.times[protein.born.index] )
  
  mrna.ith.jumps.index<- which( ith.all.data.T$jump.type == "mrna.death.event")
  
  mrna.ith.all.states  <- rbind( ith.all.data.T[1,], 
                                 ith.all.data.T[mrna.ith.jumps.index,], 
                                 ith.all.data.T[nrow(ith.all.data.T),] ) #with head and tail
  
  
  # yi.update.len = length(mrna.ith.jumps.index)+1
  yi.update.len = length(mrna.ith.jumps.index)+1
  if(yi.update.len==1){yi.update.len=2}
  
  for(i in 2:yi.update.len) {
    
    yi.minus.1 <- mrna.ith.all.states$jump.times[i-1]
    yi.plus.1 <-  mrna.ith.all.states$jump.times[i+1]
    
    
    
    which.protein.index <- protein.birth.times > yi.minus.1 & protein.birth.times <  yi.plus.1
    
    Region.i.protein.born <- protein.birth.times[which.protein.index]
    #now the region becomes [ yi.minus.1, #Region.i.protein.bron, yi.plus.1 ]  
    
    Region.i.protein.bron.number<- length(Region.i.protein.born)
    #length of subintervals is Region.i.protein.bron.number+1
    
    sub.interval.yi<-c(  yi.minus.1,Region.i.protein.born, yi.plus.1)
    
    
    #before the yi the jumps, with the X_1(y_{i-1}) weight
    weight.i.minus.mRNA <-     mrna.ith.all.states$SpeciesmRNA[i-1]
    
    #after the yi the jumps, with the X_1(y_{i-1}) weight
    weight.i.mRNA<-  mrna.ith.all.states$SpeciesmRNA[i]
    
    
    
    if(weight.i.mRNA == 0){
      #the last death event of mrna molecule need to later than protein bith events
      ratio.sub.interval<- rep(0,l=(Region.i.protein.bron.number+1) )
      ratio.sub.interval[length(ratio.sub.interval)]<- 1
      
    }else{
      ################################################################################   
      #Each subinterval have each weight
      log.total.weight.vec<-rep(NA,l=(Region.i.protein.bron.number+1) )
      
      log.subinterval.size.vec<- rep(NA,l=(Region.i.protein.bron.number+1) )
      
      #calculate the sub-interval size 
      #here we normalized the weight
      
      #####################################################################################################      
      
      for(k in 1:(Region.i.protein.bron.number+1)){
        
        # total.weight.vec[k]<- (weight.i.minus.mRNA^(k-1))*( weight.i.mRNA^(Region.i.protein.bron.number-(k-1)) )  
        
        log.total.weight.vec[k]<- (k-1)*log(weight.i.minus.mRNA)+ (Region.i.protein.bron.number-(k-1))*log( weight.i.mRNA )  
        
        log.subinterval.size.vec[k]<- log.total.weight.vec[k]+ log((exp(-a*sub.interval.yi[k])-exp(-a*sub.interval.yi[k+1]))*a^{-1})
        
        # if(is.na(log.subinterval.size.vec[k]) ){ print(log.total.weight.vec[k]) }
        
      }
      
      max.log.size<- max(log.subinterval.size.vec)
      
      log_sum_all_interval_size <- max.log.size + 
        log( sum( exp(log.subinterval.size.vec -max.log.size )) )
      
      # ratio.sub.interval<- subinterval.size.vec/sum(subinterval.size.vec)
      
      #if(is.na(log_sum_all_interval_size) ){ print(max.log.size) }
      
      
      ratio.sub.interval<- exp( log.subinterval.size.vec - log_sum_all_interval_size)
      #####################################################################################################   
    }     
    
    
    
    #update into this region
    index.yi.updated.region<- sample(x=(Region.i.protein.bron.number+1),size=1, prob=ratio.sub.interval)
    
    pos.updated.yi.region<- sub.interval.yi[index.yi.updated.region:(index.yi.updated.region+1)]
    
    #inverse sampling  
    u=runif(n=1)
    #this is new yi
    yi.can<- -a^{-1}*log( exp(-a*pos.updated.yi.region[1])-u*(exp(-a*pos.updated.yi.region[1])-exp(-a*pos.updated.yi.region[2]) ) )
    
    #put it back to the original ith interval data
    mrna.ith.all.states$jump.times[i]<- yi.can #where the mrna and protein states unchange
    
  }
  
  #update back to the ith interval data
  
  ith.all.data.T.del.mrna <- ith.all.data.T[-c(mrna.ith.jumps.index),]
  
  
  for(j in 2:(nrow(mrna.ith.all.states)-1) ){
    
    new.mrna.jumps <- mrna.ith.all.states$jump.times[j]    #without head and tail,just mrna jumps
    
    index.mrna<- findInterval( new.mrna.jumps,  as.numeric(ith.all.data.T.del.mrna$jump.times))
    
    ith.all.data.T.del.mrna <-  rbind(ith.all.data.T.del.mrna[1:index.mrna,],
                                      c( new.mrna.jumps, NA, NA, "mrna.death.event"),
                                      ith.all.data.T.del.mrna[(index.mrna+1):nrow(ith.all.data.T.del.mrna),]
    )
    
  }
  re.ith.all.data.T<- all.species.reorder.fn(ith.all.data.T= ith.all.data.T.del.mrna )
  
  return(re.ith.all.data.T)
  
}


extract.ith.data.T.fn <-function(partial.data.out, all.data.T, k){
  
  y.k1<- partial.data.out$partial.observed.time[k]
  y.k<-  partial.data.out$partial.observed.time[k+1]
  
  error<- 1e-12
  
  row.head <- cbind(partial.data.out[k,])
  row.tail<- cbind(partial.data.out[(k+1),])
  
  index.chop.all.data<-   which( all.data.T$jump.times > y.k1+error  & all.data.T$jump.times < y.k-error   )
  
  if( length(index.chop.all.data) == 0){
    ith.all.data.T <- rbind(  row.head,
                              row.tail, deparse.level = 0 )
  }else{
    
    chop.ith.all.data.T <- all.data.T[index.chop.all.data[1]:index.chop.all.data[length(index.chop.all.data)], ]
    #this the ith interval data
    colnames( row.head)<- colnames(all.data.T)
    colnames( row.tail)<- colnames(all.data.T)
    
    ith.all.data.T <- rbind(  row.head,
                              chop.ith.all.data.T,
                              row.tail, deparse.level = 1 )
    
  }
  colnames(ith.all.data.T)<- colnames(all.data.T)
  
  ith.all.data.T[,1:3]  <- data.frame(
    lapply(ith.all.data.T[,1:3], function(col) as.numeric(col)),
    check.names = FALSE
  )
  
  return(ith.all.data.T)
}

ful.cond.mrna.jumps.updated.fixed.na.fn<- function( theta, partial.data.out, all.data.T,N){
  
  all.data.T[,1:3]  <- data.frame(
    lapply(all.data.T[,1:3], function(col) as.numeric(col)),
    check.names = FALSE
  )
  theta1<- theta[1]
  theta2.o<- theta[2]
  theta3<- theta[3]
  chop.number<- length(partial.data.out[,1])-2
  updated.all.data.T<- all.data.T[1,]
  
  #for each interval protein update
  ##################################################################################################
  for(k in 1:(chop.number+1) ){
    
    #extract the all data within the ith partial time observation time 
    ith.all.data.T<- extract.ith.data.T.fn(partial.data.out = partial.data.out, all.data.T=  all.data.T, k=k )
    #give the ith data
    ##################################################################################################
    
    #no mrna jumps, no mran samples
    if( !any( ith.all.data.T$jump.type== c("mrna.death.event") ) ){
      mrna.updated.ith.all.data.T<- ith.all.data.T
    }else{
      #let us update all the ith mran jump times during each partial time observation
      mrna.updated.ith.all.data.T<-  ith.interval.mrna.samples.fixed.jumps.fn(ith.all.data.T =ith.all.data.T, theta = theta, N=N )
    }
    #print(k)
    
    colnames(mrna.updated.ith.all.data.T)<- colnames(ith.all.data.T)
    
    updated.all.data.T<- rbind(updated.all.data.T, mrna.updated.ith.all.data.T[-1, ])
    
  }
  
  return(updated.all.data.T)
}
all.species.reorder.fn<-function(ith.all.data.T){
  
  #reorder the species molecule number
  ith.all.data.T$SpeciesmRNA<- as.numeric(ith.all.data.T$SpeciesmRNA)
  ith.all.data.T$SpeciesProtein<- as.numeric(ith.all.data.T$SpeciesProtein)
  
  for(l in 2:(nrow(ith.all.data.T)-1)){
    
    if ( ith.all.data.T$jump.type[l] == "mrna.death.event" ) {
      
      ith.all.data.T$SpeciesProtein[l]<- ith.all.data.T$SpeciesProtein[l-1]
      ith.all.data.T$SpeciesmRNA[l]<- ith.all.data.T$SpeciesmRNA[l-1]-1
      
    }else{
      
      if( ith.all.data.T$jump.type[l] == "birth.event"){
        ith.all.data.T$SpeciesProtein[l]<-  ith.all.data.T$SpeciesProtein[l-1] +1 
        ith.all.data.T$SpeciesmRNA[l]<- ith.all.data.T$SpeciesmRNA[l-1]
      } 
      if(ith.all.data.T$jump.type[l] == "death.event" ){
        ith.all.data.T$SpeciesProtein[l]<-  ith.all.data.T$SpeciesProtein[l-1]  - 1 
        ith.all.data.T$SpeciesmRNA[l]<- ith.all.data.T$SpeciesmRNA[l-1]
      }
    }
  }
  ith.all.data.T[,1:3]  <- data.frame(
    lapply(ith.all.data.T[,1:3], function(col) as.numeric(col)),
    check.names = FALSE
  )
  
  return(ith.all.data.T)
  
}




protein.adding.pair.fn<- function(ith.all.data.T){
  
  #new birth and death data  
  yi.minus.1<- ith.all.data.T$jump.times[1]
  yi<-ith.all.data.T$jump.times[nrow(ith.all.data.T)] 
  
  add.sams <-  runif(n=2, min=yi.minus.1, max=yi)
  birth.time.add <- add.sams[1]
  death.time.add<- add.sams[2]
  
  
  birth.time.index <-   findInterval( birth.time.add, as.numeric( ith.all.data.T$jump.times)  )  # the ith interavl is index by [i-1,i]
  
  
  ith.all.data.T  <- rbind(ith.all.data.T[1:birth.time.index, ], 
                           c(birth.time.add, NA, NA, "birth.event"),   
                           ith.all.data.T[(birth.time.index+1):nrow(ith.all.data.T), ])
  
  #this need to be later, where the ith.all.data.T is reordered!
  
  death.time.index <-   findInterval( death.time.add, as.numeric( ith.all.data.T$jump.times)  )
  
  ith.all.data.T<- rbind(ith.all.data.T[1:death.time.index, ], 
                         c(death.time.add, NA, NA, "death.event"),   ith.all.data.T[(death.time.index+1):nrow(ith.all.data.T), ])
  
  re.ith.all.data.T<- all.species.reorder.fn(ith.all.data.T=ith.all.data.T)
  return( re.ith.all.data.T )
}

protein.deleting.pair.fn<- function(ith.all.data.T ){
  
  
  birth.index <-  which(ith.all.data.T$jump.type == "birth.event")
  death.index <-  which(ith.all.data.T$jump.type == "death.event") 
  
  
  if( length(birth.index)>1 ){
    birth.delete.sam <-  sample( birth.index, size=1)}else{  birth.delete.sam <-  birth.index }
  
  if( length(death.index) >1 ){
    death.delete.sam  <- sample(death.index, size=1)
  }else{death.delete.sam  <-  death.index  }
  
  
  ith.all.data.T<- ith.all.data.T[-c( birth.delete.sam, death.delete.sam), ]
  
  
  re.ith.all.data.T<- all.species.reorder.fn(ith.all.data.T=ith.all.data.T)
  
  return(re.ith.all.data.T)
  
}


protein.shift.pair.fn<- function( ith.all.data.T){
  
  
  index.all.events <- which( ith.all.data.T$jump.type %in% c("birth.event", "death.event"))
  
  if(length(index.all.events)==1){
    shift.index<- index.all.events
  }else{
    shift.index <- sample(c(index.all.events), size = 1)
  }
  
  ith.all.data.T$jump.times<- as.numeric(ith.all.data.T$jump.times)
  
  r.t<- runif(n=1,min=ith.all.data.T$jump.times[1],
              max= ith.all.data.T$jump.times[nrow(ith.all.data.T)])
  
  
  
  shifted.row <-  ith.all.data.T[shift.index, ] #keep this shifted row information
  shifted.row$jump.times = r.t
  shifted.row$SpeciesmRNA = NA
  
  ith.all.data.T<- ith.all.data.T[-shift.index, ]
  
  added.index<- findInterval(r.t, ith.all.data.T$jump.times    )
  
  ith.all.data.T  <- rbind(ith.all.data.T[1:added.index, ], 
                           shifted.row, #added the previous shifted one
                           ith.all.data.T[(added.index+1):nrow(ith.all.data.T), ])
  
  
  #reorder the protein value
  
  re.ith.all.data.T<- all.species.reorder.fn(ith.all.data.T=ith.all.data.T)
  
  return(re.ith.all.data.T)
  
}
log.density.protein.cond.mrna.jumps.fn<- function(ith.all.data.T,N,theta){
  
  if( !all(ith.all.data.T$SpeciesProtein >=0) ){
    data2.cond.data1= -Inf 
    return( data2.cond.data1 )
  }
  
  
  
  theta1= theta[1]
  theta2.o= theta[2]
  theta3= theta[3]
  #############################################################################################
  ith.protein.d <-  as.numeric( ith.all.data.T$SpeciesProtein )
  ith.mrna.d <- as.numeric( ith.all.data.T$SpeciesmRNA )
  ith.times<- as.numeric( ith.all.data.T$jump.times )
  int.X1.ith<- 0 
  int.X2.ith<- 0 
  
  log.prod.x1<- 0  #this may not needed, because mran does not change during the ith interval
  log.prod.x2<-0
  
  
  for(j in 2:length(ith.times)){
    int.X1.ith<- int.X1.ith +   ith.mrna.d[j-1]*( ith.times[j]- ith.times[j-1] )
    int.X2.ith<- int.X2.ith + ith.protein.d[j-1]*(ith.times[j]- ith.times[j-1] )
    
    if(ith.protein.d[j]- ith.protein.d[j-1] == +1 ){
      log.prod.x1 <-  log.prod.x1 + log( theta2.o*ith.mrna.d[j]/N )}
    
    if(ith.protein.d[j]- ith.protein.d[j-1] == -1){
      log.prod.x2<- log.prod.x2 + log(theta3*ith.protein.d[j-1] )
    }
  }
  
  if(theta1>0 &&  theta2.o>0 && theta3>0 ){
    
    term1<-  log.prod.x1 - ( theta2.o/N )*int.X1.ith 
    term2<- log.prod.x2 - theta3*int.X2.ith  
    
    data2.cond.data1<- term1+term2
    
  }else{  data2.cond.data1=-Inf }
  return( data2.cond.data1 )
}


case1.00.move.protein.fn<- function(ith.all.data.T, weight.prob ){
  
  new.ith.all.data.T<- protein.adding.pair.fn(ith.all.data.T = ith.all.data.T)
  
  #calculate the weight proposal ratio  
  ##########################################################  
  #when nb=bd=0, only adding prob.weig=1,others are 0, the reverse are not, three possibles, with weight p_2 of deleting
  
  yi.min.1<- ith.all.data.T$jump.times[1]
  yi<- ith.all.data.T$jump.times[nrow(ith.all.data.T)]
  
  proporsal.weight.ratio<- log( weight.prob[2]) + 2*log( yi - yi.min.1 ) 
  
  move.type<- "add_a_pair"
  
  return(list(new.ith.all.data.T=new.ith.all.data.T , proporsal.weight.ratio= proporsal.weight.ratio, move.type=move.type)) 
}




case2.one.jump.type.protein<-function(ith.all.data.T, weight.prob){
  
  norm.prob<-c(weight.prob[1], weight.prob[3])/(weight.prob[1]+weight.prob[3])
  
  u.2<- sample(c(1,2), size=1, prob =norm.prob) #equal of add or shift
  
  yi.min.1<- ith.all.data.T$jump.times[1]
  yi<- ith.all.data.T$jump.times[nrow(ith.all.data.T)]
  
  
  if(u.2 == 1){
    new.ith.all.data.T<-   protein.adding.pair.fn(ith.all.data.T = ith.all.data.T)
    
    #calculate the weight proposal ratio  
    ##########################################################    
    new.nb <- length(which( new.ith.all.data.T$jump.type == "birth.event"))
    new.nd <- length(which( new.ith.all.data.T$jump.type == "death.event"))
    
    proporsal.weight.ratio<- log( weight.prob[2]/norm.prob[1]) + 2*log( yi - yi.min.1 )- log(new.nb*new.nd )
    move.type<- "add_a_pair"
    
    ##########################################################       
  }else{
    
    new.ith.all.data.T<- protein.shift.pair.fn( ith.all.data.T = ith.all.data.T)
    #calculate the weight proposal ratio  
    ##########################################################     
    proporsal.weight.ratio<- log(1)
    ##########################################################
    move.type<- "shift_a_event"
  }
  
  
  return(list(new.ith.all.data.T= new.ith.all.data.T, proporsal.weight.ratio= proporsal.weight.ratio,  move.type=  move.type)) 
  
}

case3.ge1move.protein.fn<- function(ith.all.data.T, weight.prob ){
  
  three.rv <- sample(c(1,2,3), size=1, prob= weight.prob)
  norm.prob= c(weight.prob[1], weight.prob[3])/(weight.prob[1]+weight.prob[3])
  
  yi.min.1<- ith.all.data.T$jump.times[1]
  yi<- ith.all.data.T$jump.times[nrow(ith.all.data.T)]
  
  #adding 
  if( three.rv == 1 ){
    # print("add a pair.by.three.possibles")
    
    new.ith.all.data.T<- protein.adding.pair.fn(ith.all.data.T = ith.all.data.T)
    new.nb <- length(which( new.ith.all.data.T$jump.type == "birth.event"))
    new.nd <- length(which( new.ith.all.data.T$jump.type == "death.event"))
    
    #calculate the weight proposal ratio  
    ##########################################################    
    proporsal.weight.ratio<- log( weight.prob[2]/weight.prob[1]) + 2*log(yi -yi.min.1)- log(new.nb*new.nd )}
  ##########################################################     
  move.type<- "add_a_pair"
  
  #shifting      
  if( three.rv == 3 ){
    
    #print("shift a event.by.three.possibles")
    new.ith.all.data.T<- protein.shift.pair.fn( ith.all.data.T = ith.all.data.T)
    #calculate the weight proposal ratio  
    ##########################################################  
    proporsal.weight.ratio<- log(1)
    ########################################################## 
    move.type<- "shift_a_event"
  }
  
  #deleting   
  if( three.rv == 2 ){
    
    # print("delete a pair.by.three.possibles")
    new.ith.all.data.T<- protein.deleting.pair.fn( ith.all.data.T = ith.all.data.T) 
    
    #calculate the weight proposal ratio  
    ##########################################################
    old.nb <- length(which( ith.all.data.T$jump.type == "birth.event"))
    old.nd <- length(which( ith.all.data.T$jump.type == "death.event"))
    move.type<- "delete_a_pair"
    
    if( old.nb== 1 && old.nd ==1){
      norm.prob.ratio <-  log(1/weight.prob[2])
    }else{
      if( (old.nb >=2 && old.nd==1 ) || (old.nb == 1 && old.nd>= 2 )  ){
        norm.prob.ratio <- log(norm.prob[1]/weight.prob[2])
      }
      
      if( old.nb>= 2 && old.nd >=2 ){
        norm.prob.ratio<- log(weight.prob[1]/ weight.prob[2])
      }
      
    }
    
    proporsal.weight.ratio<-  norm.prob.ratio +  log(old.nb*old.nd)- 2*log(yi -yi.min.1)
    
  }
  
  
  return(list(new.ith.all.data.T= new.ith.all.data.T, proporsal.weight.ratio= proporsal.weight.ratio, move.type=move.type)) 
}


protein.all.data.update.fn <- function(all.data.T, partial.data.out, weight.prob, theta, N){
  
  updated.all.data.T <- all.data.T[1,] #store all the updated data
  # norm.prob= c(weight.prob[1], weight.prob[3])/(weight.prob[1]+weight.prob[3]) 
  chop.number<- length(partial.data.out[,1])-2
  #for each interval protein update
  count=0
  
  interval.protein.birth.number<- rep(NA,l=(nrow(partial.data.out)-1) )
  
  
  events.move <- as.data.frame( matrix(0, nrow= 1, ncol=3) )
  colnames(events.move)<- c("add","shift","delete")
  
  ##################################################################################################
  for(k in 1:(chop.number+1) ){
    
    #kth interval is from the index of k:k+1
    ith.all.data.T<- extract.ith.data.T.fn(partial.data.out = partial.data.out, k=k, all.data.T = all.data.T)
    
    #let us start the three protein move types: ####################################################################################################################################################################################################
    #case 1:(0,0)
    cond1<- all(ith.all.data.T$jump.type %in% c("Partial.obse.times", "mrna.death.event") )
    #case 3: (>=1,>=1), there are three possible moves
    cond3<- all( c("birth.event", "death.event") %in% ith.all.data.T$jump.type ) 
    #case2: their complements
    cond2<- !(cond1 |cond3)
    
    
    #this the all data within the kth interval of partially observed data
    #if there is no events of protein, only add a pair is possible
    if(cond1  ){
      case1.out<- case1.00.move.protein.fn(ith.all.data.T= ith.all.data.T, weight.prob=weight.prob)
      new.ith.all.data.T<-  case1.out$new.ith.all.data.T
      new.proporsal.weight.ratio<- case1.out$proporsal.weight.ratio
      propose.type<- case1.out$move.type
      ##########################################################
    } 
    
    #case 2: (0,>=1) or (>=1, 0)  
    #if there is no death/birth events of protein during this interval, only adding or shift weight equally
    ##############################################################################################################
    
    if(cond2){
      case2.out<-  case2.one.jump.type.protein(ith.all.data.T= ith.all.data.T, weight.prob=weight.prob)
      new.ith.all.data.T<-  case2.out$new.ith.all.data.T
      new.proporsal.weight.ratio<- case2.out$proporsal.weight.ratio
      propose.type<- case2.out$move.type
    } 
    
    
    #########################################################################################################
    if( cond3 ){
      case3.out<- case3.ge1move.protein.fn(ith.all.data.T= ith.all.data.T, weight.prob=weight.prob)
      new.ith.all.data.T<-  case3.out$new.ith.all.data.T
      new.proporsal.weight.ratio<- case3.out$proporsal.weight.ratio
      propose.type<- case3.out$move.type
    }
    
    
    #finished cases
    
    #finished the proposal move, let us decide to accept it or not
    ##################################################################################################
    ##################################################################################################
    #MH of protein update
    #old ith density
    #input mrna data into the density of protein birth events
    old.d2.cond.d1.density <- log.density.protein.cond.mrna.jumps.fn( ith.all.data.T = ith.all.data.T, theta= theta, N=N)
    
    #new ith density
    new.d2.cond.d1.density <- log.density.protein.cond.mrna.jumps.fn( ith.all.data.T = new.ith.all.data.T, theta= theta, N=N)
    
    density.log.ratio <-  new.d2.cond.d1.density-  old.d2.cond.d1.density
    
    ##################################################################################################
    
    u.rv<- log( runif(n=1) )
    
    alpha.i<-  density.log.ratio +  new.proporsal.weight.ratio  
    
    if(u.rv < alpha.i ){
      current.ith<- new.ith.all.data.T
      #print(paste("accept the new", k,"th interval data"))
      count= count+1
      
      if( propose.type ==  "add_a_pair"){ events.move$add<- events.move$add +1 }
      if(propose.type == "shift_a_event" ){events.move$shift <-events.move$shift +1 }
      if(propose.type == "delete_a_pair"  ){ events.move$delete <- events.move$delete +1  }
      
      
    }else{
      # print(paste("Not accept the new", k,"th interval data"))
      current.ith<-ith.all.data.T
    }
    
    
    ##################################################################################################
    ##################################################################################################
    interval.protein.birth.number[k]<- length(which(current.ith$jump.type == "birth.event"))
    ##################################################################################################
    
    updated.all.data.T<- rbind(updated.all.data.T, current.ith[2:nrow(current.ith), ] )
    
  }
  
  #finshed ith interval updates
  acceptance.protein<- count/(chop.number+1)
  # print(paste("accpetance rate of protein for one realization is", acceptance.protein))
  
  #this include the partial observed time
  return(list(updated.all.data.T=updated.all.data.T, 
              acceptance.protein= acceptance.protein,
              interval.protein.birth.number=interval.protein.birth.number,
              events.move= events.move ))
}


one.mcmc.run.fn<- function(all.data.T, partial.data.out, weight.prob, para.iter, initial.para.states, N, alp1,alp2, alp3, bet1, bet2,bet3 ){
  
  #input the scaled theta, where theta2 is not theta2/N here
  
  one.para.out<- matrix(NA, nrow= (para.iter+1), ncol=3)
  one.para.out[1,]<- initial.para.states
  
  total.events.move.out<- as.data.frame( matrix(NA, nrow= para.iter, ncol=3) )
  colnames(total.events.move.out)<- c("add","shift","delete")
  
  #start each data, one row of parameter corresponding to one of data matrix
  store.all.data.T<- vector("list", (para.iter+1))
  store.all.data.T[[1]]<- all.data.T
  
  #acceptance.rate.protein<- rep(NA,l=para.iter)
  protein.birth.number<- rep(NA,l=para.iter)
  #start from true protein birth number
  protein.birth.number[1]<- length(which(all.data.T$jump.type == "birth.event"))
  
  total.interval.protein.birth.number<- matrix(NA, nrow=para.iter, ncol= (nrow(partial.data.out)-1) )
  
  for(i in 2:(para.iter+1)){
    #step 1: Gibbs of paramters 
    #conditional on the [0,T] of full-observed data
    #gibbs update: which is consider the full trajectory not in the ith interval!
    ##################################################################################################
    #this is the info for the summary of data to input the  the full-conditional of parameter 
    full.data.summary.out.T.value<- cond.para.based.full.data.fn(all.data.T = all.data.T)
    
    
    #this is the function of full conditional of parameter out given the full conditional data
    theta.all.can<- full.cond.theta.gamma.fn(full.data.summary.out.T= full.data.summary.out.T.value,
                                             alp1=alp1, alp2=alp2, bet1 = bet1, bet2 = bet2,alp3=alp3, bet3=bet3, N=N)
    
    
    one.para.out[i,]<- theta.all.can
    ##################################################################################################
    
    
    #step 2: Gibbs of mRNA
    #conditional on the [0,T] of protein birth data & N jumps of the total reaction times of mRNA data 
    #because the mRNA jumps is conditional on the protein birth events
    #2.1 input of protein birth data(X2.data)
    ##################################################################################################
    
    mrna.upated.all.data.T <- ful.cond.mrna.jumps.updated.fixed.na.fn(theta= theta.all.can, partial.data.out = partial.data.out, all.data.T = all.data.T,  N=N)
    ##################################################################################################
    #step3: MH of protein  
    #this is the update of all intervals 
    ##################################################################################################
    
    
    protein.all.data.T <- protein.all.data.update.fn(all.data.T= mrna.upated.all.data.T, partial.data.out=partial.data.out, weight.prob=weight.prob, theta=theta.all.can, N=N )
    
    
    all.data.T<-  protein.all.data.T$updated.all.data.T  #update the current species data
    
    total.events.move.out[i-1, ] <- protein.all.data.T$events.move
    
    #store the acceptance of protein
    acceptance.rate.protein[i-1]<- protein.all.data.T$acceptance.protein
    #each one parameter, store one data set 
    
    #store the ith interval protein birth number
    total.interval.protein.birth.number[i-1,]<- protein.all.data.T$interval.protein.birth.number
    #store all the data
    store.all.data.T[[i]]<- all.data.T
    #store the total protein birth number
    protein.birth.number[i]<- length(which(all.data.T$jump.type == "birth.event"))
    
  }
  
  one.sample.mean<- colMeans(one.para.out)
  
  return(list(store.all.data.T=store.all.data.T, 
              one.para.out=one.para.out,
              one.sample.mean=one.sample.mean, 
              acceptance.rate.protein= acceptance.rate.protein,
              protein.birth.number= protein.birth.number,
              total.interval.protein.birth.number= total.interval.protein.birth.number ,
              total.events.move.out= total.events.move.out) )
}

#############################################################################################################################

protein.all.data.update.simple.para.acceptance.fn<- function(all.data.T, partial.data.out, weight.prob, theta, N){
  
  updated.all.data.T <- all.data.T[1,] #store all the updated data
  
  chop.number<- length(partial.data.out[,1])-2
  #for each interval protein update
  count=0
  
  # interval.protein.birth.number<- rep(NA,l=(nrow(partial.data.out)-1) )
  
  
  events.move <- as.data.frame( matrix(0, nrow= 1, ncol=3) )
  colnames(events.move)<- c("add","shift","delete")
  
  ##################################################################################################
  for(k in 1:(chop.number+1) ){
    
    #kth interval is from the index of k:k+1
    ith.all.data.T<- extract.ith.data.T.fn(partial.data.out = partial.data.out, k=k, all.data.T = all.data.T)
    
    #let us start the three protein move types: #####################################################################################################################################################
    #case 1:(0,0)
    cond1<- all(ith.all.data.T$jump.type %in% c("Partial.obse.times", "mrna.death.event") )
    #case 3: (>=1,>=1), there are three possible moves
    cond3<- all( c("birth.event", "death.event") %in% ith.all.data.T$jump.type ) 
    #case2: their complements
    cond2<- !(cond1 |cond3)
    
    
    #this the all data within the kth interval of partially observed data
    #if there is no events of protein, only add a pair is possible
    if(cond1  ){
      case1.out<- case1.00.move.protein.fn(ith.all.data.T= ith.all.data.T, weight.prob=weight.prob)
      new.ith.all.data.T<-  case1.out$new.ith.all.data.T
      new.proporsal.weight.ratio<- case1.out$proporsal.weight.ratio
      propose.type<- case1.out$move.type
      ##########################################################
    } 
    
    #case 2: (0,>=1) or (>=1, 0)  
    #if there is no death/birth events of protein during this interval, only adding or shift weight equally
    ##############################################################################################################
    
    if(cond2){
      case2.out<-  case2.one.jump.type.protein(ith.all.data.T= ith.all.data.T, weight.prob=weight.prob)
      new.ith.all.data.T<-  case2.out$new.ith.all.data.T
      new.proporsal.weight.ratio<- case2.out$proporsal.weight.ratio
      propose.type<- case2.out$move.type
    } 
    
    
    #########################################################################################################
    if( cond3 ){
      case3.out<- case3.ge1move.protein.fn(ith.all.data.T= ith.all.data.T, weight.prob=weight.prob)
      new.ith.all.data.T<-  case3.out$new.ith.all.data.T
      new.proporsal.weight.ratio<- case3.out$proporsal.weight.ratio
      propose.type<- case3.out$move.type
    }
    
    
    #finished cases
    
    #finished the proposal move, let us decide to accept it or not
    ##################################################################################################
    ##################################################################################################
    #MH of protein update
    #old ith density
    #input mrna data into the density of protein birth events
    old.d2.cond.d1.density <- log.density.protein.cond.mrna.jumps.fn( ith.all.data.T = ith.all.data.T, theta= theta, N=N)
    
    #new ith density
    new.d2.cond.d1.density <- log.density.protein.cond.mrna.jumps.fn( ith.all.data.T = new.ith.all.data.T, theta= theta, N=N)
    
    density.log.ratio <-  new.d2.cond.d1.density-  old.d2.cond.d1.density
    
    ##################################################################################################
    
    u.rv<- log( runif(n=1) )
    
    alpha.i<-  density.log.ratio +  new.proporsal.weight.ratio  
    
    if(u.rv < alpha.i ){
      current.ith<- new.ith.all.data.T
      #print(paste("accept the new", k,"th interval data"))
      count= count+1
      
      if( propose.type ==  "add_a_pair"){ events.move$add<- events.move$add +1 }
      if(propose.type == "shift_a_event" ){events.move$shift <-events.move$shift +1 }
      if(propose.type == "delete_a_pair"  ){ events.move$delete <- events.move$delete +1  }
      
      
    }else{
      # print(paste("Not accept the new", k,"th interval data"))
      current.ith<-ith.all.data.T
    }
    
    
    ##################################################################################################
    ##################################################################################################
    # interval.protein.birth.number[k]<- length(which(current.ith$jump.type == "birth.event"))
    ##################################################################################################
    
    updated.all.data.T<- rbind(updated.all.data.T, current.ith[2:nrow(current.ith), ] )
    
  }
  
  #finshed ith interval updates
  acceptance.protein<- count/(chop.number+1)
  # print(paste("accpetance rate of protein for one realization is", acceptance.protein))
  
  #this include the partial observed time
  return(list(updated.all.data.T=updated.all.data.T, 
              acceptance.protein= acceptance.protein,
              # interval.protein.birth.number=interval.protein.birth.number,
              events.move= events.move ))
}
#############################################################################################################################
 
one.mcmc.run.simple.para.acceptance.fn<- function(all.data.T, partial.data.out, weight.prob, para.iter, initial.para.states, N, alp1,alp2, alp3, bet1, bet2,bet3, return.data.index ){
  
  #input the scaled theta, where theta2 is not theta2/N here
  
  one.para.out<- matrix(NA, nrow= (para.iter+1), ncol=3)
  one.para.out[1,]<- initial.para.states
  colnames( one.para.out)<- c("theta1", "theta2", "theta3")
  
  total.events.move.out<- as.data.frame( matrix(NA, nrow= para.iter, ncol=3) )
  colnames(total.events.move.out)<- c("add","shift","delete")
  
  #start each data, one row of parameter corresponding to one of data matrix
  store.all.data.T<- vector("list", (para.iter+1))
  store.all.data.T[[1]]<- all.data.T
  
  acceptance.rate.protein<- rep(NA,l=para.iter)
  protein.birth.number<- rep(NA,l=para.iter)
  #start from true protein birth number
  protein.birth.number[1]<- length(which(all.data.T$jump.type == "birth.event"))
  
  total.interval.protein.birth.number<- matrix(NA, nrow=para.iter, ncol= (nrow(partial.data.out)-1) )
  
  for(i in 2:(para.iter+1)){
    #step 1: Gibbs of paramters 
    #conditional on the [0,T] of full-observed data
    #gibbs update: which is consider the full trajectory not in the ith interval!
    ##################################################################################################
    #this is the info for the summary of data to input the  the full-conditional of parameter 
    full.data.summary.out.T.value<- cond.para.based.full.data.fn(all.data.T = all.data.T)
    
    
    #this is the function of full conditional of parameter out given the full conditional data
    theta.all.can<- full.cond.theta.gamma.fn(full.data.summary.out.T= full.data.summary.out.T.value,
                                             alp1=alp1, alp2=alp2, bet1 = bet1, bet2 = bet2,alp3=alp3, bet3=bet3, N=N)
    
    
    one.para.out[i,]<- theta.all.can
    ##################################################################################################
    
    
    #step 2: Gibbs of mRNA
    #conditional on the [0,T] of protein birth data & N jumps of the total reaction times of mRNA data 
    #because the mRNA jumps is conditional on the protein birth events
    #2.1 input of protein birth data(X2.data)
    ##################################################################################################
    
    mrna.upated.all.data.T <- ful.cond.mrna.jumps.updated.fixed.na.fn(theta= theta.all.can, partial.data.out = partial.data.out, all.data.T = all.data.T,  N=N)
    ##################################################################################################
    #step3: MH of protein  
    #this is the update of all intervals 
    ##################################################################################################
    
    
    protein.all.data.T <- protein.all.data.update.simple.para.acceptance.fn(all.data.T= mrna.upated.all.data.T, partial.data.out=partial.data.out, weight.prob=weight.prob, theta=theta.all.can, N=N )
    
    
    all.data.T<-  protein.all.data.T$updated.all.data.T  #update the current species data
    
    total.events.move.out[i-1, ] <- protein.all.data.T$events.move
    
    #store the acceptance of protein
    acceptance.rate.protein[i-1]<- protein.all.data.T$acceptance.protein
    #each one parameter, store one data set 
    
    #store the ith interval protein birth number
    #total.interval.protein.birth.number[i-1,]<- protein.all.data.T$interval.protein.birth.number
    #store all the data
    
    if( return.data.index == TRUE){
      store.all.data.T[[i]]<- all.data.T
    }
    
    #store the total protein birth number
    # protein.birth.number[i]<- length(which(all.data.T$jump.type == "birth.event"))
    
  }
  
  #  one.sample.mean<- colMeans(one.para.out)
  
  return(list( store.all.data.T=store.all.data.T, 
               one.para.out=one.para.out,
               #one.sample.mean=one.sample.mean, 
               acceptance.rate.protein= acceptance.rate.protein,
               # protein.birth.number= protein.birth.number,
               # total.interval.protein.birth.number= total.interval.protein.birth.number ,
               total.events.move.out= total.events.move.out) )
}
 
given.one.data.mcmc.run.simple.fn<- function( given.result.out, chop.number.times.region, N, true.theta,start.theta,  para.iter,weight.prob, know.reaction.final.index, T.protein.unknown ){
  
  alp1=0.0001
  alp2=0.0001
  bet1=0.0001
  bet2=0.0001
  alp3=0.0001
  bet3=0.0001
  
  
  #fixed one data
  #one.par.mean.out<- matrix(NA, nrow=1, ncol=3)
  #true.posterior.paras<- one.par.mean.out
  #store.par.mean.out.chops<- vector("list",length(chop.number.times.region))
  #theta.ratio.mean.out<- rep(NA, l=length(chop.number.times.region))
  
  
  # for(it in 1:length(chop.number.times.region)){
  #   store.par.mean.out.chops[[it]]<- one.par.mean.out
  # }
  
  
  #start data
  
  for(k in 1:length(chop.number.times.region)){
    
    ################################################################   
    data.p.out <- chop.time.equally.fn(all.data = given.result.out, chop.number.times=chop.number.times.region[k],  know.reaction.final.index=know.reaction.final.index, T.protein.unknown=T.protein.unknown )
    #after chop, out put the partial time data and the true full observed data
    partial.data.out<- data.p.out$partial.all.species.data  
    #this is the partial time observation
    #start true of all species states  and parameter
    all.data.T<- data.p.out$with.T.all.data 
    #############################################################  
    #this is true posterior parameter
    
    true.data.T.summary.out<- cond.para.based.full.data.fn(all.data.T =all.data.T )
    true.posterior.paras<- full.cond.theta.gamma.fn(full.data.summary.out.T = true.data.T.summary.out,
                                                    alp1=alp1, alp2=alp2, bet1 = bet1, bet2 = bet2,alp3=alp3, bet3=bet3, N=N)
    
    
    
    start <- Sys.time()
    
    
    one.mcmc.out<- one.mcmc.run.simple.para.acceptance.fn(para.iter= para.iter,
                                                          all.data.T = all.data.T, partial.data.out=partial.data.out, 
                                                          initial.para.states= start.theta,
                                                          weight.prob=weight.prob, N=N,
                                                          alp1=alp1, alp2=alp2, bet1 = bet1, bet2 = bet2,alp3=alp3, bet3=bet3, return.data.index = FALSE)
    
    end <- Sys.time()
    total.time <- difftime(end, start, units = "auto")
    cat("For simulation time is", format(total.time), "\n","with N=", N)
    
    
    #save the sample mean
    #store.par.mean.out.chops[[k]][j,]<-  one.mcmc.out$one.sample.mean
    
    # theta.ratio.mean.out[k]<- mean( one.mcmc.out$one.para.out[,2]/one.mcmc.out$one.para.out[,3] )
    
    print(paste("finish chop", chop.number.times.region[k], "one data iter" ))
    
    ####################################################################################################  
    #save the plots to a folder  
    para.file <- paste0("para.out-chop.times=",chop.number.times.region[k] , ",N=", N, ".csv")
    write.csv(one.mcmc.out$one.para.out, file = para.file , row.names = FALSE)
    cat("CSV saved to", normalizePath(para.file ), "\n")
    ####################################################################################################  
    acceptance.file<- paste0("acceptance.out-chop.times=", chop.number.times.region[k], ",N=", N, ".csv")
    write.csv(  one.mcmc.out$acceptance.rate.protein , file = acceptance.file , row.names = FALSE)
    cat("CSV saved to", normalizePath(acceptance.file ), "\n")
    ####################################################################################################  
    events.move.file<- paste0("events-move-chop.times=", chop.number.times.region[k], ",N=", N, ".csv")
    write.csv(one.mcmc.out$total.events.move.out, file =  events.move.file , row.names = FALSE)
    cat("CSV saved to", normalizePath( events.move.file ), "\n")
    ####################################################################################################  
    #  data.file<- paste0("data.out-chop.times=", chop.number.times.region[k], ".pdf")
    #   data2.file<- paste0("Protein.out-chop.times=", chop.number.times.region[k], ".pdf")
    #acceptance.file<- paste0("acceptance.out-chop.times=", chop.number.times.region[k], ".pdf")
    #  protein.nb.data.file<- paste0("protein.nb.data-chop.times=", chop.number.times.region[k], ".pdf")
    #   ith.interval.protein.data.file<-paste0("ith.interval.protein.data-chop.times", chop.number.times.region[k], ".pdf")
    
    #theta.ratio.file<- paste0("theta.ratio-chop.times=", chop.number.times.region[k], ".pdf")
    
    ## choose your day-month pattern once
    #today_tag <- format(Sys.time(), "%d%b-%H")     # "18Jun-1123" (adds hour+min)
    
    ## build the full name inside your loop
    #folder.file <- paste0(today_tag, "-One-Data, N=", N)
    
    
    # plotProteinMRNA(total.out =  one.mcmc.out, para.file = para.file, 
    #                 data.file =data.file,data2.file=data2.file,
    #                  folder.file=folder.file, 
    #                  true.posterior.paras=true.posterior.paras, true.theta = true.theta , N=N , k=chop.number.times.region[k], 
    #                  #acceptance.file =acceptance.file,
    #                  protein.nb.data.file= protein.nb.data.file, 
    #                  ith.interval.protein.data.file= ith.interval.protein.data.file,
    #                  theta.ratio.file= theta.ratio.file)
    
    
    
  }
  
  
  # return(list(#store.par.mean.out.chops=store.par.mean.out.chops,
  #   true.posterior.paras=true.posterior.paras, 
  #   theta.ratio.mean.out=theta.ratio.mean.out ))
  
}


true.theta=  c(0.5,80,1)

start.theta<-c( 2,2,2)

chop.number.times.region<- c(7)
  
para.iter= 50000
weight.prob<- c(1/3,1/3,1/3)
N=1000



csv_dir   <- sprintf("/home/pmxyg3/partial-26Jun-n%d" , N)


file_name <- paste0("one_given_data_N=",N, "_",
                    format(Sys.Date(), "%Y%m%d"), ".csv")

## 3.  Full path = dir + file name
full_path <- file.path(csv_dir, file_name)      # safer than paste()
file.exists(full_path)   # should return TRUE

one.data.out<- read.csv(
  file             = full_path,     # absolute path on the server
  header = TRUE,               # first row contains column names
  check.names = FALSE,         # <-- keeps names exactly as-is
  stringsAsFactors = FALSE     # keep character columns as character
)


print("start mcmc")



total.out<- given.one.data.mcmc.run.simple.fn(given.result.out= one.data.out,
                                       true.theta = true.theta,
                                       start.theta= start.theta,
                                       chop.number.times.region=chop.number.times.region,  N=N, 
                                       para.iter=para.iter,
                                       weight.prob = weight.prob, T.protein.unknown=10, 
                                       know.reaction.final.index = FALSE
) 

 