
#### calculate transition counts for each reads ####
#########original############
trans<- function(seq, seqDetail, m=1){  # transition counts
  allcomb<- expand.grid(rep(list(c("A", "C", "T", "G")), m+1))
  end<- dim(seq)[1]
  len<- nchar(seq[1,1])  # length of reads
  dic<- matrix(0, end, 4^(m+1))
  colnames(dic)<- apply(allcomb, 1, paste, collapse="")
  colnames(dic)<- colnames(dic)[order(colnames(dic))]
  ind<- 1
  for(i in 1:end){
    tempSeqDetail<- seqDetail[[i]]
    for(j in 1:(len-m)){
      tuple_char<- tempSeqDetail[j:(j+m)]
      # original
      tuple<- paste(tuple_char, collapse = "")
      idx<- which(colnames(dic) == tuple)
      dic[i, idx] = dic[i, idx] + 1
      
    }
  }
  return(dic)
}

#############original+reverse#########

#### calculate transition counts for each reads ####
# map characters ATCG to its counter-part
map_char<- function(tuple_char){
  for(i in 1:length(tuple_char)){
    tuple_char[i]<- switch(tuple_char[i], T = "A", A = "T", C = "G", G = "C")
  }
  return(tuple_char)
}
trans_4read<- function(seq, seqDetail, m=1){  # transition counts
  allcomb<- expand.grid(rep(list(c("A", "C", "T", "G")), m+1))
  end<- dim(seq)[1]
  len<- nchar(seq[1,1])  # length of reads
  dic<- matrix(0, end, 4^(m+1))
  colnames(dic)<- apply(allcomb, 1, paste, collapse="")
  colnames(dic)<- colnames(dic)[order(colnames(dic))]
  ind<- 1
  for(i in 1:end){
    tempSeqDetail<- seqDetail[[i]]
    for(j in 1:(len-m)){
      tuple_char<- tempSeqDetail[j:(j+m)]
      # original
      tuple<- paste(tuple_char, collapse = "")
      idx<- which(colnames(dic) == tuple)
      dic[i, idx] = dic[i, idx] + 1
      # reversed
      tuple<- paste(rev(tuple_char), collaps = "")
      idx<- which(colnames(dic) == tuple)
      dic[i, idx] = dic[i, idx] + 1
      
      # switch charcters
      tuple_char<- map_char(tuple_char)     
      # original
      tuple<- paste(tuple_char, collapse = "")
      idx<- which(colnames(dic) == tuple)
      dic[i, idx] = dic[i, idx] + 1
      # reversed
      tuple<- paste(rev(tuple_char), collaps = "")
      idx<- which(colnames(dic) == tuple)
      dic[i, idx] = dic[i, idx] + 1
      
    }
  }
  return(dic)
}

                          
#####change counts to probability
prop<- function(x){
  if(sum(x) == 0){
    return(x*0)
  }
  else{
    return(x/sum(x))
  }
}

dicprop<- function(counts){
  dic_prop<- unlist(tapply(counts, rep(1:(2*2^m), each=4), prop))
  return(dic_prop)
}

starprop<- function(counts){
  star_counts<- unlist(tapply(counts, rep(1:(2*2^m), each=4), sum))
  star_prop<-star_counts/sum(star_counts)
  return(t(star_prop))
}

######caculate mode for each read###########
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#######function##############

#1######prior########           #######change#####
prior<-function(z,read){
  p<-vector(mode="list", length=n_sample)
  group<- max(z)
  for(k in 1:n_sample){
    groupCol<- unique(z[z[,k] != 0, k])
    p[[k]]<- matrix(alpha/(alpha+read-1)/(group+1-length(groupCol)), group+1)
    for(j in groupCol){
      p[[k]][j]<-sum(z[,k] == j)/(read-1+alpha) 
    }  
  }
  return(p)
}
#2#####propability for read occurance######
#transition probability#                  ########change#######
trans_p<-function(z,read){
  p<-vector(mode="list", length=n_sample)
  group<- max(z)
  for(k in 1:n_sample){
    p[[k]]<-array(0, dim=c(group+1,4^(m+1)))    
    for (j in 1:(group+1)){
      clusindex<-which(z[,k]==j)
      if(length(clusindex) > 0){
        p[[k]][j,]<-dicprop(colSums(trans_counts_4read[c(clusindex,read),])) 
      } else{
        p[[k]][j,]<-dicprop(trans_counts_4read[read,]) 
      } 
    }
 }
return(p)
}

#starting point propability#            #######change#######
star_p<-function(z,read){
  p<-vector(mode="list", length=n_sample)
  group<- max(z)
  for(k in 1:n_sample){    
    p[[k]]<-array(0, dim=c(group+1,4))
    for(j in 1:(group+1)){
      clusindex<-which(z[,k]==j)
      if(length(clusindex) > 0){
        p[[k]][j,]<-starprop(colSums(trans_counts_4read[c(clusindex,read),]))  
    }else{
      p[[k]][j,]<-starprop(trans_counts_4read[read,])
    }
  }
  }  
  return(p)
}

########SMC posterior ########
startwhere<-function(read,k){
  
  if (seqDetail[[read]][1]=="A"){
    p = star_p(z,read)[[k]][1]
  } else if (seqDetail[[read]][1]=="C"){
    p=star_p(z,read)[[k]][2]
  }else if (seqDetail[[read]][1]=="G"){
    p=star_p(z,read)[[k]][3]
  }else {
    p=star_p(z,read)[[k]][4]
  }
return(p)  
}

                                            #######change#######
post_p_temp<- function(z,read){
  group<- max(z)
  p<- vector(mode="list", length=n_sample)
  trans_p_tmp<- trans_p(z, read)
  for(k in 1:n_sample){
    p[[k]]<- matrix(0, group+1)
    start_p<- startwhere(read,k)
    for(j in 1:(group+1)){
      p[[k]][j]<- w[k]*start_p*prod(trans_p_tmp[[k]][j,]^trans_counts[read,])  
    }
  }
  p<- Reduce("+", p)/length(p)
  return(p)
}
                                         #########change#########
post_p<- function(z,read){
  group<- max(z)
  p<-matrix(0,group+1,n_sample)
    for (k in 1:n_sample){
      for(j in 1:(group+1)){
      p[j,k]<-w[k]*prior(z,read)[[k]][j]*post_p_temp(z,read)[j]
    }
}
#print(p)
output<-rowMeans(p)/sum(rowMeans(p))
return(output)
}

######update weight#####
weight<-function(z,read){
  p<-rep(0,n_sample)
  for (k in 1:n_sample){
    group<-dim(table(z[,k],exclude=0))
    start_p<-startwhere(read,k)
  #  print(k)
   # print(read)
   p[k]<- w[k]*start_p*prod(trans_p(z,read)[[k]][z[read,k],]^trans_counts[read,])
  }
  output<-p/sum(p)
  return(output)
}


args<-c("abundance_species_equal.txt","abundance")
m<-1
set.seed(0)
#####read data######
seq<-read.table(args[1],colClasses = "character")
seqDetail<- apply(seq, 2, strsplit, split = "")$V1
index<-c(1:5,60000:60004,6,100000)
#index<-c(sample(1:100,9),90000)
seqDetail<- seqDetail[index]
######prameters########
trans_counts_4read<- trans_4read(as.matrix(seq[index,]), seqDetail, 1)
trans_counts<-trans(as.matrix(seq[index,]), seqDetail, 1)
n_read<-dim(trans_counts)[1]
n_sample<-100 #number of samples for each read
alpha<-0.5  #alpha prameter
z<-matrix(0,ncol=n_sample,nrow=n_read) # hidden variabke
w<-matrix(0,ncol=n_sample,nrow=n_read) #important weight
z[1,]=1
w<-rep(1/n_sample,n_sample)

######SMC#####
for(i in 2:length(seqDetail)){
  posterior<-post_p(z,i)
  z[i,]<- sample.int(length(posterior), n_sample, replace=T, prob=posterior)
  print(z[i])
  w<-weight(z,i) ##update weight
  print(i)
}

for(i in 1:n_read){
  group[i]<-Mode(z[i,])
}
