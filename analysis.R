library("dplyr")
library(DescTools)
library(patchwork) # To display 2 charts together
source("C://DataCopied/Research/R/createShiny.R")

cluQuality<-function(dt,sim,ND=TRUE){
  temp<-merge(sim,dt,by.y = c("seqID"),by.x=c("bseq"))
  names(temp)[names(temp) == 'cluster'] <- 'clusterB'
  temp<-merge(temp,dt,by.y = c("seqID"),by.x=c("aseq"))
  names(temp)[names(temp) == 'cluster'] <- 'clusterA'
  
  # temp<-temp%>%filter(clusterA==clusterB & aseq!=bseq)
  # clu<-temp%>%group_by(clusterA)%>%
  #   summarise(avg_sim=mean(identity.),
  #             size=length(unique(bseq))+1)
  
  temp<-temp%>%filter(clusterA==clusterB)
  clu<-temp%>%group_by(clusterA)%>%
    summarise(avg_sim=mean(identity.),
              size=length(unique(bseq)))
  
  names(clu)[names(clu) == 'clusterA'] <- 'clus_id'
  clu$clu<-seq.int(nrow(clu))-1
  
  if (ND){
    nodeDistance<-read.csv(paste(path,"/nodeDistance",".txt",sep=""))
    nodeDistance<-nodeDistance%>%group_by(clu)%>%
      summarise(avg_nodeDistance=mean(HD))
    clu<-merge(clu,nodeDistance,by.x="clus_id",by.y="clu")
  }
  clu
}



######################## iterations ##################################
path<-"C://DataCopied/Research/tree/data/toy/toy-k9-w100-s5-s60-l20"
cluQ<-NULL
nodeDistance<-NULL
water<-read.csv(r"(C:\DataCopied\Research\tree\data\toy\toy-waterall.csv)")

for (r in seq(0,9)){
  dt<-plotEnt(paste(path,"/clusters-r",r,".txt",sep=""))
  temp<-cluQuality(dt,water)
  temp$run<-r
  cluQ<-rbind(cluQ,temp)
  
  temp<-read.csv(paste(path,"/nodeDistance-r",r,".txt",sep=""))
  temp$run<-r
  nodeDistance<-rbind(nodeDistance,temp)
  # ggplotly(ggplot(nodeDistance)+
  #            geom_density(aes(x=HD))
  #   # geom_density(aes(x=HD,colour=as.factor(clu)))
  # )
  # mean(nodeDistance$HD)
  # mean(cluQ$avg_sim)
}
meanNodeDistance<-nodeDistance%>%group_by(run)%>%summarise(HD=mean(HD))
meanCluQ<-cluQ%>%group_by(run)%>%summarise(avg_sim=mean(avg_sim))
save.image("C:/DataCopied/Research/tree/data/toy/clusQ.RData")

load("C:/DataCopied/Research/tree/data/toy/clusQ-r.RData")
meanNodeDistance_r<-meanNodeDistance
meanCluQ_r<-meanCluQ
load("C:/DataCopied/Research/tree/data/toy/clusQ.RData")
ggplot()+
  geom_line(data=meanNodeDistance,aes(x=run,y=HD,color="ordered"))+
  geom_line(data=meanNodeDistance_r,aes(x=run,y=HD,color = "random"))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))

ggplot()+
  geom_line(data=meanCluQ,aes(x=run,y=avg_sim,color="ordered"))+
  geom_line(data=meanCluQ_r,aes(x=run,y=avg_sim,color = "random"))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))

ggplot()+
  geom_line(data=meanNodeDistance,aes(x=run,y=HD,color="HD"))+
  geom_line(data=meanCluQ,aes(x=run,y=avg_sim-60,color="avg_sim"))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_y_continuous(
    name = "HD",
    sec.axis = sec_axis(~.+60, name="avg_sim")
  )
######################## toy ##################################
source("C://DataCopied/Research/R/toyFunctions.R")
sigClust<-plotEnt(paste(path,"/sigClust-k9-c48.txt",sep=""))
cluQ_sigClust<-cluQuality(sigClust,water,4,FALSE)

{
dt<-plotEnt(paste(path,"/toy-k9-w100-s5-s80-l18.txt",sep=""))
cluQ_2WMT<-cluQuality(dt,water,4,FALSE)
}

# mean(cluQ$avg_sim)
# avg<- cluQ %>% summarise(across(everything(), list(mean)))
avg<-sumAll(cluQ_2WMT,cluQ_sigClust)
plotTwo(cluQ_2WMT,cluQ_sigClust,3)
plotTwo(normSize(cluQ_2WMT),normSize(cluQ_sigClust),ncol(cluQ_2WMT)+1)
plotTwo(normSize(cluQ_2WMT, 'sd_sim', TRUE),
        normSize(cluQ_sigClust, 'sd_sim', TRUE),
        ncol(cluQ_2WMT)+1)
test<-normSize(cluQ_2WMT, 'sd_sim')
(
  ggplot(cluQ_2WMT)+
  geom_point(aes(x=clu,y=avg_sim,size=size))+
    # geom_point(aes(x=clu,y=avg_nodeDistance,size=size))+
    ylim(50,100)
  )

(
  ggplot(cluQ)+
    geom_line(aes(x=clu,y=avg_sim,colour="avg_sim"))+
    geom_line(aes(x=clu,y=avg_nodeDistance,colour="avg_nodeDistance"))+
    ylim(0,100)+
    
    # scale_y_continuous(
    #   name = "avg_sim",
    #   sec.axis = sec_axis(~.-60, name="avg_nodeDistance")
    # )+
    xlim(0,57)
)

# save.image("C:/DataCopied/Research/tree/data/toy/clusQ-f89e5e9-unbound.RData")

load("C:/DataCopied/Research/tree/data/toy/clusQ-1Ambi.RData")

######################## hierarchy ###############################
hierarchy<-read.csv(paste(path,"/hierarchy.txt",sep=""))
hierarchy<-hierarchy%>%mutate(parent=case_when((rank == 0 & parent == 0 & !child %in% parent) ~ child, TRUE ~ parent))
t<-hierarchy%>%filter(!child %in% parent)

dt<-plotEnt(paste(path,"/toy-k9-w100-s5-s60-l20.txt",sep=""))
dt<-merge(dt,hierarchy,by.x="cluster",by.y="child")


cluQualityParent<-function(dt,sim){
  temp<-merge(sim,dt,by.y = c("seqID"),by.x=c("bseq"))
  temp<-merge(temp,dt,by.y = c("seqID"),by.x=c("aseq"))
  
  temp<-temp%>%filter(parent.x==parent.y & aseq!=bseq)

  # clu<-temp%>%group_by(cluster.y)%>%
  clu<-temp%>%group_by(parent.y)%>%
    summarise(avg_sim=mean(identity.),
              size=length(unique(bseq))+1,
              rank=max(rank.y))
  names(clu)[names(clu) == 'parent.y'] <- 'clus_id'
  clu$clu<-seq.int(nrow(clu))-1
  clu
}
# temp<-cluQualityParent(dt,water)
cluQ<-cluQualityParent(dt,water)

cluQ<-temp%>%group_by(parent.y)%>%
  summarise(avg_sim=min(identity.))

(
  ggplot(cluQ)+
    geom_point(aes(x=clu,y=avg_sim,size=size))
)


load("C:/DataCopied/Research/tree/data/toy/clusQ-1.RData")


######################## staph ################################
source("C://DataCopied/Research/R/staphFunctions.R")

path<-"C://DataCopied/Research/tree/data/staphopia-contigs"
w<-read.csv(paste(path,"/sample-k9-w100-s5-global_sim.txt",sep=""))
water<-waterStaph(paste(path,"/sample_needle.txt",sep=""),
                  paste(path,"/sample_meta.csv",sep=""))


# dt<-read.csv(paste(path,"/sample-k9-w100-s5-s100-l59.txt",sep=""),header=FALSE)
dt<-read.csv(paste(path,"/sample-k9-w100-s5-s90-l60.txt",sep=""),header=FALSE)
colnames(dt)<-c("seqID","cluster")
# cluQ_2WMT<-cluQualityStaph(dt,water)
# cluQ_2WMT<-cluQuality(dt,w,3,FALSE)
cluQ_2WMT<-cluQuality(dt,water,2,FALSE,w,3)

# avg<- cluQ_2WMT %>% summarise(across(everything(), list(mean)))

sigClust<-read.csv(paste(path,"/sigClust-k9-c",nrow(cluQ_2WMT),".txt",sep=""),header=FALSE)
# sigClust<-read.csv(paste(path,"/sigClust.txt",sep=""),header=FALSE)
colnames(sigClust)<-c("seqID","cluster")
# cluQ_sigClust<-cluQualityStaph(sigClust,water,FALSE)
# cluQ_sigClust<-cluQuality(sigClust,w,3,FALSE)
cluQ_sigClust<-cluQuality(sigClust,water,2,FALSE,w,3)

plotTwo(cluQ_2WMT,cluQ_sigClust,3)
plotTwo(normSize(cluQ_2WMT),normSize(cluQ_sigClust),ncol(cluQ_2WMT)+1)
plotTwo(normSize(cluQ_2WMT, 'sd_sim', TRUE),
        normSize(cluQ_sigClust, 'sd_sim', TRUE),
        ncol(cluQ_2WMT)+1)
test<-normSize(cluQ_sigClust, 'sd_sim')

# plotTwo(normSize(cluQ_2WMT,'avg_sim_sig'),normSize(cluQ_sigClust,'avg_sim_sig'),ncol(cluQ_2WMT)+1)

# plotStaph(cluQ_2WMT)


avg<-sumAll(cluQ_2WMT,cluQ_sigClust)

# find nearest neighbours of the outliers
outliers_NN<-function(sim,var, col_i="i", col_j="j"){
  getSubsetBySize<-function(dt,summ,size=-1){
    if (size==-1){
      (dt%>%filter(cluster==(summ%>%slice(which.max(size)))$cluster[1]))$seqID
    }else{
      (dt%>%filter(cluster==(summ%>%slice(which(.data$size == size)))$cluster[1]))$seqID
      
    }
  }
  t_s<-getSubsetBySize(sigClust,cluQ_sigClust)
  t_w<-getSubsetBySize(dt,cluQ_2WMT)
  outliers<-setdiff(t_w,t_s)
  
  t<-sim%>%filter(.data[[col_i]] %in% outliers)
  outliers_best_inter<-t%>%filter(.data[[col_j]] %in% t_s)%>%
    group_by(.data[[col_i]])%>%
    slice(which(min(.data[[var]])<100))
  if (nrow(outliers_best_inter)>0){
    names(outliers_best_inter)[names(outliers_best_inter) == col_j] <- paste("best_inter_",col_j,sep="")
    names(outliers_best_inter)[names(outliers_best_inter) == var] <- paste("best_inter_",var,sep="")
    # colnames(outliers_best_inter)[2:3]<-c("inter_j",paste("best_inter_",var,sep=""))

    true_outliers<-outliers_best_inter%>%pull(col_i)
# 
    NN<-t%>%filter(.data[[col_i]] %in% true_outliers & !.data[[col_j]] %in% t_w)%>%
      group_by(.data[[col_i]])%>%
      slice(which(max(.data[[var]])==.data[[var]]))
    names(NN)[names(NN) == var] <- paste("best_intra_",var,sep="")
    # colnames(NN)[which(colnames(NN)==var)]<-paste("best_intra_",var,sep="")
    
    NN<-merge(NN,dt,by.x=col_j,by.y="seqID")
    NN<-merge(distinct(NN,.data[[col_i]],cluster,.keep_all = TRUE),outliers_best_inter,by=col_i)
    NN<-merge(NN,cluQ_2WMT)
    NN[!grepl('\\.[x|y]$', colnames(NN))]
  }else{
    print("no outliers based on 2WMT signature")
    NULL
  }
}

NN<-outliers_NN(w,"similarity")
# N2<-outliers_NN(water,"identity.","seqID.x","seqID.y")
N2<-outliers_NN(water,"identity.")


plotStaph(cluQ_2WMT,TRUE,seq(3,2))
plotStaph(cluQ_2WMT,FALSE,seq(3,4))

plotStaph(cluQ_sigClust,FALSE)

plotStaph(cluQ_2WMT)+ ggtitle("2W-MT")+
  theme(legend.position="none")+
  plotStaph(cluQ_sigClust)+ ggtitle("sigClust")


# cluQSim<-cluQualityStaph(dt,sim)
# plotStaph(cluQ_2WMT)+ theme(legend.position="none")+
#   plotStaph(cluQ_sigClust)

t<-arrangeByEnt(cluQ_2WMT)
plotTwo(arrangeByEnt(cluQ_2WMT),arrangeByEnt(cluQ_sigClust),2)

temp<-water%>%group_by(seqID.x)%>%
  summarise(size=length(unique(seqID.y)))%>%
  filter(size!=length(dt[,1]))

######################## ecoli ########################
waterEcoli<-function(simFile,skip=0){
  water<-read.csv(simFile,skip = skip)
  # water$clu<-seq.int(nrow(water))-1
  len<-(sqrt(8*length(water[,1])+1)-1)/2
  idx<-NULL
  idy<-NULL
  for (i in 0:(len-1)){
    idx<-c(idx,c(i:(len-1)))
    idy<-c(idy,rep(i,len-i))
  }
  water$aseq=idx
  water$bseq=idy
  # merge(water,meta,by.x="bseq",by.y="name")
  # meta<-read.csv(metaFile,header=FALSE)
  # water<-water%>%mutate(aseq=meta[match(aseq,meta$V2),1],
  #                       bseq=meta[match(bseq,meta$V2),1])
  # water<-mirror(water)
  mirror(water)
}

cluQuality<-function(dt,sim,ND=TRUE){
  temp<-merge(sim,dt,by.y = c("seqID"),by.x=c("bseq"))
  names(temp)[names(temp) == 'cluster'] <- 'clusterB'
  temp<-merge(temp,dt,by.y = c("seqID"),by.x=c("aseq"))
  names(temp)[names(temp) == 'cluster'] <- 'clusterA'
  
  # temp<-temp%>%filter(clusterA==clusterB & aseq!=bseq)
  # clu<-temp%>%group_by(clusterA)%>%
  #   summarise(avg_sim=mean(identity.),
  #             size=length(unique(bseq))+1)
  
  temp<-temp%>%filter(clusterA==clusterB)
  clu<-temp%>%group_by(clusterA)%>%
    summarise(avg_sim=mean(identity.),
              size=length(unique(bseq)))
  
  names(clu)[names(clu) == 'clusterA'] <- 'clus_id'
  clu<-clu%>% arrange(desc(avg_sim),desc(size))
  clu$clu<-seq.int(nrow(clu))-1
  
  if (ND){
    nodeDistance<-read.csv(paste(path,"/nodeDistance",".txt",sep=""))
    nodeDistance<-nodeDistance%>%group_by(clu)%>%
      summarise(avg_nodeDistance=mean(HD))
    clu<-merge(clu,nodeDistance,by.x="clus_id",by.y="clu")
    clu<-clu%>% arrange(clu)
  }
  clu
}


path<-"C://DataCopied/Research/tree/data/ecoli"
water<-waterEcoli(paste(path,"/gene-2500-5000.waterall",sep=""),17)
dt<-read.csv(paste(path,"/gene-2500-5000-k9-w100-s5-s40-l5.txt",sep=""),header=FALSE)
colnames(dt)<-c("seqID","cluster")

cluQ_2WMT<-cluQuality(dt,water)
# avg<- cluQ_2WMT %>% summarise(across(everything(), list(mean)))


sigClust<-read.csv(paste(path,"/sigclust-k9.txt",sep=""),header=FALSE)
colnames(sigClust)<-c("seqID","cluster")
cluQ_sigClust<-cluQuality(sigClust,water,FALSE)
avg<-sumAll(cluQ_2WMT,cluQ_sigClust)


(
  ggplot(cluQ_sigClust)+
  # ggplot(cluQ_2WMT)+
    geom_point(aes(x=clu,y=avg_sim,size=size))+
    # geom_point(aes(x=clu,y=avg_nodeDistance,size=size))+
    ylim(50,100)
)


#############################################################

######################## old ###############################
formatResultMeta<-function(file, metafile){
  meta<-read.csv(metafile, header = F)
  
  dt<-read.csv(file)
  # dt$species_id<-meta$V14
  dt<-dt%>%mutate(species_id=meta[dt$seqID==meta$V1,"V14"])
  ancestors<-unique(dt$ancestor)
  dt<-dt%>%group_by(cluster)%>%
    mutate(ancestor=match(ancestor,ancestors),
           clu_size=length(cluster))
  
  dt<-dt%>%group_by(ancestor)%>%
    mutate(anc_clu_size=length((ancestor)))
  dt
}


# formatResultSigClustMeta<-function(file, metafile){
#   meta<-read.csv(metafile, header = F)
#   
#   dt<-read.csv(file, header = F)
#   colnames(dt)<-c("seqID","cluster")
#   dt$species_id<-meta$V14
#   dt<-dt%>%group_by(cluster)%>%
#     mutate(clu_size=length(cluster))
#   dt
# }
# 
# 
# {
# sigclust<-formatResultSigClustMeta(r"(C:\DataCopied\Research\tree\data\controlled-silva\sigclust-116.csv)",
#                                    metafile)
# ent_clu_sigclust<-sigclust%>%group_by(cluster)%>%
#   summarise(entropy=Entropy(table(species_id),base=exp(1)),
#             size=mean(clu_size))
# mean(ent_clu_sigclust$entropy)
# 
# ggplot()+
#   geom_point(data=ent_clu_sigclust,aes(x=cluster,y=entropy))
# }

# truth<-read.csv(r"(C:\Users\n9417770\OneDrive - Queensland University of Technology\phd\sourceCode\postKTree\controlled-silva\acgt\needle\AAAA02042586.1650.3157.needle)",
#                 skip=18)
# ggplot(truth,aes(x=bseq,y=identity.))+geom_point()

{
  # file<-gsub('"', "", gsub("\\\\", "/", readClipboard()))
  # file<-r"(C:\DataCopied\Research\tree\data\controlled-silva\controlled-silva-k9-w50-s5-s80-l10.txt)"
  # metafile<-r"(C:\DataCopied\Research\tree\data\controlled-silva\controlled-silva-meta.csv)"
  # tree<-formatResultMeta(file,metafile)
  
  file<-("C://DataCopied/Research/tree/data/toy/toy-k9-w100-s5-s60-l20.txt")
  
  # tree<-formatResult("C://DataCopied/Research/tree/data/toy/toy-k9-w100-s5-s50-l30.txt")
  # tree<-formatResult("C://DataCopied/Research/tree/data/toy/toy-k9-w100-s5-s50-l30-random.txt")
# search<-formatResult("C://DataCopied/Research/tree/data/data-s2-l5-search.txt")

singleton<-tree[tree$clu_size==1,]

ent_clu<-tree%>%group_by(cluster)%>%
  summarise(entropy=Entropy(table(species_id),base=exp(1)),
            size=mean(clu_size))


ent_anc<-tree%>%group_by(ancestor)%>%
  summarise(entropy=Entropy(table(species_id),base=exp(1)),
            size=mean(anc_clu_size),
            depth=max(level))

(
  ggplot()+
    geom_point(data=ent_anc,aes(x=ancestor,y=entropy,size=size))
  # geom_point(data=ent_mkm,aes(x=cluster,y=entropy,size=size))+
  #ylim(0,2.5)+ scale_size(limit = c(0,100))
)

ent<-tree%>%group_by(species_id)%>%
  summarise(entropy=Entropy(table(ancestor),base=exp(1)),
            size=length(seqID))
ggplot()+
  geom_point(data=ent,aes(x=species_id,y=entropy))

pure<-ent_anc[ent_anc$entropy==0&ent_anc$size==26,]
mean(ent_anc$entropy)
# mean(ent_clu$entropy)

}



formatSim<-function(needleFile,simFile,metaFile,skip=0){
  needle<-read.csv(needleFile,check.names = F, skip=skip)
  sim<-read.csv(simFile)
  temp<-read.csv(metaFile,header=F)
  colnames(temp)<-c("i","bseq")
  n<-merge(temp,needle)
  colnames(temp)<-c("j","aseq")
  n<-merge(temp,n)
  sim<-left_join(sim,n)
  sim
}



sim1<-formatSim(r"(C:\DataCopied\Research\tree\data\controlled-silva\sample_needle.csv)",
               r"(C:\DataCopied\Research\tree\data\controlled-silva\sample-k9-w100-s5_sim-global.txt)",
               r"(C:\DataCopied\Research\tree\data\controlled-silva\sample-meta.csv)")

sim2<-formatSim(r"(C:\DataCopied\Research\tree\data\controlled-silva\sample_needle.csv)",
               r"(C:\DataCopied\Research\tree\data\controlled-silva\sample-k9-w100-s5_sim.txt)",
               r"(C:\DataCopied\Research\tree\data\controlled-silva\sample-meta.csv)")

sim3<-formatSim(r"(C:\DataCopied\Research\tree\data\controlled-silva\sample_needle.csv)",
                r"(C:\DataCopied\Research\tree\data\controlled-silva\sample_k9_sim.txt)",
                r"(C:\DataCopied\Research\tree\data\controlled-silva\sample-meta.csv)")

# colnames(sim1)[which(colnames(sim1)=="similarity")]<-"similarity_global"
# colnames(sim2)[which(colnames(sim2)=="similarity")]<-"similarity_windows"
# colnames(sim3)[which(colnames(sim3)=="similarity")]<-"similarity_all"
# allsim<-merge(sim1,sim2)
# allsim<-merge(allsim,sim3)
# ggplot(allsim)+
#   geom_point(aes(x=similarity_all,y=similarity_global))+
#   xlim(0,100)+ylim(0,100)+ coord_fixed(ratio = 1)


sim<-NULL
sim<-rbind(sim,sim1%>%mutate(tag="global"))
sim<-rbind(sim,sim2%>%mutate(tag="window"))
sim<-rbind(sim,sim3%>%mutate(tag="all"))
ggplot(sim)+
  geom_point(aes(x=`similarity%`,y=similarity,colour=tag))+
  facet_wrap(~tag)+
  xlim(0,100)+ylim(0,100)+ coord_fixed(ratio = 1)

temp<-sim3[sim3$`similarity%`>=98,]

temp<-sim2[sim2$`similarity`>=50&sim2$`similarity`<=51,]

ggplot(temp)+geom_point(aes(x=`gap%`,y=`similarity%`))

sim2[sim2$bseq=="GU120661.1.1496" & sim2$aseq=="AM085476.1.1485",]

(
ggplot()+
    # geom_point(data=sim1,aes(x=`similarity%`,y=similarity_global,colour="global minimisers"),alpha=0.5)+
    # geom_point(data=sim2,aes(x=`similarity%`,y=similarity_windows,colour="windows minimisers"),alpha=0.5)+
    # geom_point(data=sim3,aes(x=`similarity%`-`gap%`,y=similarity_all,colour="all kmers - gap"),alpha=0.5)+
    # geom_point(data=sim3,aes(x=`similarity%`,y=similarity_all,colour="all kmers"),alpha=0.5)+
    
    
  # geom_point(data=sim1,aes(x=`similarity%`,y=similarity,colour="windows minimisers, s=5"),alpha=0.5)+
  # geom_point(data=sim2,aes(x=`similarity%`,y=similarity,colour="windows minimisers, s=10"),alpha=0.5)+
  # geom_point(data=allkmer,aes(x=`similarity%`,y=similarity,colour="all kmers"),alpha=0.5)+
    xlim(0,100)+ylim(0,100)+ coord_fixed(ratio = 1)+
  xlab("NW similarity")+ylab("Jaccard Similarity with MinimiserSet")
  
)

p<-ggplot()+xlim(0,100)+ylim(0,100)+ coord_fixed(ratio = 1)+
  xlab("NW similarity")+ylab("Jaccard Similarity with MinimiserSet")

plot_grid
####################################################################
library(stringr)
formatSimFolder<-function(metaFile,needleFile,skip,folder, pattern){
  meta<-read.csv(metaFile,header = FALSE)
  needle<-read.csv(needleFile,check.names = F,skip=skip)
  files <- list.files(path=folder, pattern=pattern, full.names=TRUE, recursive=FALSE)
  sim<-NULL
  for (f in files){
    temp<-read.csv(f)
    temp<-temp%>%mutate(k=as.numeric(str_extract(f,"[:digit:]+")))
    sim<-rbind(sim,temp)
  }
  colnames(meta)<-c("i","bseq")
  n<-merge(meta,needle)
  colnames(meta)<-c("j","aseq")
  n<-merge(meta,n)
  sim<-left_join(sim,n)
  sim$k<-as.factor(sim$k)
  sim
}


dt_local<- formatSimFolder(r"(C:\DataCopied\Research\tree\data\controlled-silva\sample-meta.csv)",
                      r"(C:\DataCopied\Research\tree\data\controlled-silva\sample.waterall)",
                      14,
                      r"(C:\DataCopied\Research\tree\data\controlled-silva\)",
                      "sample-.*-all_sim.txt")
dt_local<-dt_local%>%mutate(tag="local")

dt_global<-formatSimFolder(r"(C:\DataCopied\Research\tree\data\controlled-silva\sample-meta.csv)",
                           r"(C:\DataCopied\Research\tree\data\controlled-silva\sample_needle.csv)",
                           0,
                           r"(C:\DataCopied\Research\tree\data\controlled-silva\)",
                           "sample-.*-all_sim.txt")
dt_global<-dt_global%>%mutate(tag="global")

dt<-rbind(dt_global,dt_local)


ggplot(dt[dt$`similarity%`==100,])+
  geom_point(aes(x=length,y=score,colour=tag),alpha=0.5)+
  facet_wrap(~tag)


(
  # ggplot(dt[dt$i<10,])+
  ggplot(dt[dt$k==9,])+
    geom_point(aes(x=score,y=`similarity%`,colour=tag),alpha=0.5)+
    facet_wrap(~tag) +
    ylim(0,100)+ 
    xlab("Alignment Score")+ylab("Jaccard Similarity with all 9-kmers")
)


(
# ggplot(dt[dt$i<10,])+
  ggplot(dt)+
    geom_point(aes(x=`similarity%`,y=similarity,colour=k),alpha=0.5)+
    facet_wrap(~tag) +
    # geom_point(aes(x=`similarity%`,y=similarity,colour=k),alpha=0.5)+
  #   facet_wrap(~k) +
  # geom_smooth(aes(x=`similarity%`,y=similarity,colour=k),method=loess, se=FALSE, linetype='dashed')+
  xlim(0,100)+ylim(0,100)+ coord_fixed(ratio = 1)+
    
  xlab("NW similarity")+ylab("Jaccard Similarity with all kmers")
)

#################
sim_local<-formatSimFolder(r"(C:\DataCopied\Research\tree\data\controlled-silva\sample-meta.csv)",
                           r"(C:\DataCopied\Research\tree\data\controlled-silva\sample.waterall)",
                           14,
                           r"(C:\DataCopied\Research\tree\data\controlled-silva\)",
                           "sample-k9-w100-s5_sim.txt")
sim_local<-sim_local%>%mutate(tag="local")
sim_global<-formatSimFolder(r"(C:\DataCopied\Research\tree\data\controlled-silva\sample-meta.csv)",
                           r"(C:\DataCopied\Research\tree\data\controlled-silva\sample_needle.csv)",
                           0,
                           r"(C:\DataCopied\Research\tree\data\controlled-silva\)",
                           "sample-k9-w100-s5_sim.txt")
sim_global<-sim_global%>%mutate(tag="global")
sim<-rbind(sim_local,sim_global)
ggplot(sim)+
  geom_point(aes(x=`similarity%`,y=similarity),alpha=0.5)+
  facet_wrap(~tag)+
  xlim(0,100)+ylim(0,100)+ coord_fixed(ratio = 1)+
  xlab("alignment similarity")+ylab("Jaccard Similarity with windows minimiser")

#############################
allkmer<-formatSim(r"(C:\DataCopied\Research\tree\data\controlled-silva\sample_needle.csv)",
                r"(C:\DataCopied\Research\tree\data\controlled-silva\sample_k9_sim.txt)")
ggplot(allkmer)+
  geom_point(aes(x=`similarity%`,y=similarity))+
  xlab("NW similarity")+ylab("Jaccard Similarity with MinimiserSet")+
  xlim(0,100)+ylim(0,100)+ coord_fixed(ratio = 1)

##########################################################
silva<-read.csv(r"(C:\DataCopied\Research\tree\data\controlled-silva\sample.histo)",
             header = FALSE, sep=" ")
colnames(silva)<-c("kmerCount","freq")
silva<-silva%>%mutate(scaledkmerCount=kmerCount/max(kmerCount),
                      freq=freq/sum(freq))


staph<-read.csv(r"(C:\DataCopied\Research\tree\data\staphopia-contigs\sample.histo)",
                header = FALSE, sep=" ")
colnames(staph)<-c("kmerCount","freq")
staph<-staph%>%mutate(scaledkmerCount=kmerCount/max(kmerCount),
                      freq=freq/sum(freq))

shiny(
ggplot()+
  geom_bar(data=silva,aes(kmerCount,y=freq,fill="silva"),stat="identity",alpha=0.5)+
  geom_bar(data=staph,aes(kmerCount,y=freq,fill="staph"),stat="identity",alpha=0.5)
)


shiny(
  ggplot()+
    geom_bar(data=silva,aes(scaledkmerCount,y=freq,fill="silva"),stat="identity",alpha=0.5)+
    geom_bar(data=staph,aes(scaledkmerCount,y=freq,fill="staph"),stat="identity",alpha=0.5)
)
##########################################################
sim<-formatSim(r"(C:\DataCopied\Research\tree\data\toy\toy_needle.csv)",
               r"(C:\DataCopied\Research\tree\data\toy\toy-k9-w100-s5_sim.txt)")


sim<-sim%>%mutate(g_i=floor(i/26),
                  g_j=floor(j/26))

g<-0
temp<-sim_local[sim_local$g_i==g&sim_local$g_j==g,]

shiny(
  ggplot(temp)+
  geom_point(aes(x=`similarity%`,y=similarity,colour=as.factor(floor(i/5))))+
  facet_wrap(~as.factor(floor(j/5)))+
  # facet_wrap(~g_i)+
  # facet_wrap(~paste(g_i,g_j))+
  # xlim(0,100)+ylim(0,100)+ coord_fixed(ratio = 1)+
  xlab("NW similarity")+ylab("Jaccard Similarity with BF")
)
formatSim<-function(needleFile,simFile,metaFile,skip=0){
  needle<-read.csv(needleFile,check.names = F, skip=skip)
  needle<-unique(needle)
  sim<-read.csv(simFile)
  temp<-read.csv(metaFile,header=F)
  colnames(temp)<-c("i","bseq")
  n<-merge(temp,needle)
  colnames(temp)<-c("j","aseq")
  n<-merge(temp,n)
  sim<-left_join(sim,n)
  sim
}

sim_local<-formatSim(r"(C:\DataCopied\Research\tree\data\toy\toy.waterall)",
                     r"(C:\DataCopied\Research\tree\data\toy\toy-k9-w100-s5_sim.txt)",
                     r"(C:\DataCopied\Research\tree\data\toy\toy-meta.csv)",
                     14)
sim_local<-sim_local%>%mutate(g_i=floor(i/26),
                  g_j=floor(j/26))





temp_local<-sim_local[sim_local$g_i<1&sim_local$g_j<9,]
ggplot(sim_local)+
  geom_point(aes(x=`similarity%`,y=similarity,colour=g_i==g_j))+
  # facet_wrap(~g_i)+
  # facet_wrap(~paste(g_i,g_j))+
  xlab("Water similarity")+ylab("Jaccard Similarity with BF")+
  xlim(0,100)+ylim(0,100)+ coord_fixed(ratio = 1)






formatSimMeta<-function(needleFile,simFile, metaFile,skip=0){
  needle<-read.csv(needleFile,check.names = F,skip = skip)
  sim<-read.csv(simFile)
  meta<-read.csv(metaFile)
  tempCols<-colnames(meta)
  # temp<-data.frame(i=seq(0,max(sim$i)),aseq=needle$bseq[1:(max(sim$i)+1)])
  # n<-merge(temp,needle)
  # colnames(temp)<-c("j","bseq")
  # n<-merge(temp,n)
  # sim<-left_join(sim,n)
  colnames(meta)<-paste(tempCols, "i", sep = "_")
  sim<-merge(sim,meta,by.x = "i", by.y = "seqID_i");
  
  colnames(meta)<-paste(tempCols, "j", sep = "_")
  sim<-merge(sim,meta,by.x = "j", by.y = "seqID_j");
  colnames(sim)[which(colnames(sim)=="name_i")]<-"bseq"
  colnames(sim)[which(colnames(sim)=="name_j")]<-"aseq"

  left_join(needle,sim)
  
}

allsim<-formatSimMeta(r"(C:\DataCopied\Research\tree\data\staphopia-contigs\sample_needle.txt)",
                    r"(C:\DataCopied\Research\tree\data\staphopia-contigs\sample-k9-all-all_sim.txt)",
                    r"(C:\DataCopied\Research\tree\data\staphopia-contigs\sample_meta.csv)",
                    15)


ggplot(allsim,aes(x=`similarity%`,y=similarity,colour=((mlst_i==mlst_j))))+
  geom_point()+
  xlab("NW similarity")+ylab("Jaccard Similarity with all 9mers")+
  xlim(0,100)+ylim(0,100)+ coord_fixed(ratio = 1)

sim<-formatSimMeta(r"(C:\DataCopied\Research\tree\data\staphopia-contigs\sample_needle.txt)",
               r"(C:\DataCopied\Research\tree\data\staphopia-contigs\sample-k9-w100-s5_sim.txt)",
               r"(C:\DataCopied\Research\tree\data\staphopia-contigs\sample_meta.csv)",
               15)



(
  # ggplot(sim,aes(x=`similarity%`,y=similarity,colour=(paste(mlst_i,mlst_j))))+

  ggplot(sim,aes(x=`similarity%`,y=similarity,colour=((mlst_i==mlst_j))))+
    geom_point()+
  xlab("NW similarity")+ylab("Jaccard Similarity with BF")+
    xlim(0,100)+ylim(0,100)+ coord_fixed(ratio = 1)
  # geom_smooth(method=loess, se=FALSE, col='purple', linetype='dashed')
)

sim<-read.csv("C://DataCopied/Research/tree/data/sim-short.txt")
colnames(sim)<-c("aseq","bseq","window_similarity")
nw<-read.csv("C://DataCopied/Research/tree/data/short/short.csv",check.names = F)
dt<-merge(nw,sim)
ggplot(dt,aes(x=`similarity%`,y=window_similarity))+
  geom_point()+
  geom_smooth(method=loess, se=FALSE, col='purple', linetype='dashed')

############################################################
# Ecoli

formatSimMeta<-function(needleFile,simFile, metaFile,skip=0){
  needle<-read.csv(needleFile,check.names = F,skip = skip)
  sim<-read.csv(simFile)
  meta<-read.csv(metaFile)
  tempCols<-colnames(meta)
  # temp<-data.frame(i=seq(0,max(sim$i)),aseq=needle$bseq[1:(max(sim$i)+1)])
  # n<-merge(temp,needle)
  # colnames(temp)<-c("j","bseq")
  # n<-merge(temp,n)
  # sim<-left_join(sim,n)
  colnames(meta)<-paste(tempCols, "i", sep = "_")
  sim<-merge(sim,meta,by.x = "i", by.y = "seqID_i");
  
  colnames(meta)<-paste(tempCols, "j", sep = "_")
  sim<-merge(sim,meta,by.x = "j", by.y = "seqID_j");
  colnames(sim)[which(colnames(sim)=="name_i")]<-"bseq"
  colnames(sim)[which(colnames(sim)=="name_j")]<-"aseq"
  
  left_join(needle,sim)
}

geneLength<-read.csv(r"(C:\DataCopied\Research\tree\data\ecoli\gene-length.csv)",
                     header = F)


hist(geneLength$V1)

ggplot(geneLength, aes(x=V1)) + 
  geom_histogram(binwidth=50,color="black", fill="white")+xlab("Ecoli geneLength")

ecoli<-formatSimMeta(r"(C:\DataCopied\Research\tree\data\ecoli\gene.needleall)",
                       r"(C:\DataCopied\Research\tree\data\ecoli\gene-all-k9_sim.txt)",
                       r"(C:\DataCopied\Research\tree\data\ecoli\gene-meta.csv)",
                       14)

ecoli_local<-formatSimMeta(r"(C:\DataCopied\Research\tree\data\ecoli\gene.waterall)",
                      r"(C:\DataCopied\Research\tree\data\ecoli\gene-all-k9_sim.txt)",
                      r"(C:\DataCopied\Research\tree\data\ecoli\gene-meta.csv)",
                      14)



ecoli_local_short<-formatSimMeta(r"(C:\DataCopied\Research\tree\data\ecoli\gene-2500-5000.waterall)",
                           r"(C:\DataCopied\Research\tree\data\ecoli\gene-all-k9_sim-2500-5000.txt)",
                           r"(C:\DataCopied\Research\tree\data\ecoli\gene-meta.csv)",
                           0)

temp<-ecoli[ecoli$score>5000&ecoli$score<6000,]

# ggplot(ecoli_local_short,aes(x=`similarity%`,y=similarity))+
ggplot(ecoli_local_short,aes(x=score,y=similarity))+
  geom_point()+
  # xlim(0,100)+ coord_fixed(ratio = 1)+
  ylim(0,100)+xlim(0,43500)+
  ggtitle("geneLength btw 2500-5000")+
  xlab("Water score")+ylab("Water Similarity")
  # xlab("NW similarity")+ylab("Jaccard Simil arity with MinimiserSet")
  # xlab("NW similarity")+ylab("Jaccard Similarity with windowsBF")
  
####################################################

toy<-read.csv(r"(C:\DataCopied\Research\tree\data\toy\toy-waterall.csv)")
hist(toy$identity.)
toy[toy$bseq==15,]

ggplot(toy)+geom_tile(aes(x=aseq,y=bseq,fill=identity.))
ggplot(toy[1:20000,])+geom_point(aes(x=aseq,y=bseq))


ggplot(toy)+
  geom_density(aes(x=identity.,colour=as.factor(floor(aseq/26))),
                 alpha=0.6, position = 'identity')
require(scales)
sig<-read.csv(r"(C:\DataCopied\Research\tree\data\toy\toy-k9-w100-s5-test.sim-global_sim.txt)")
seed<-sig[sig$i%%26==0 & sig$j%%26==0,]
p<-ggplot(seed)+
  # geom_density(aes(x=similarity))
  # geom_bar(aes(x=similarity,y = (..count..)/sum(..count..)))
  # geom_bar(aes(x=similarity,y = ..count..))
  geom_histogram(aes(x=similarity))

ggplot(seed, aes(x = similarity)) +
  geom_histogram(aes(y = stat(density))) +
  scale_y_continuous(labels = scales::percent_format())
q<-ggplot_build(p)

summary(sig$similarity)
sd(toy$similarity)

ggplot(sig)+
  geom_histogram(aes(x=similarity,colour=as.factor(floor(i/26))),
               alpha=0.6, position = 'identity')

ggplot(sig)+
  geom_histogram(aes(x=similarity))+facet_wrap(~floor(i/26))



################################################
file<-"toy_26-k9-all_sim"
file<-"toy_26-k9-w50-s5_chunk-global_sim"
file<-"toy-k9-all_sim"
file<-"toy-k9-w100-s5-global_sim"
path<-paste("C:\\DataCopied\\Research\\tree\\data\\toy\\",file,".txt",sep="")
tall<-mirror(read.csv(path))
# tall<-mirror(read.csv("D:\\phd\\tree\\data\\refseq\\toy_26-k9-all_sim.txt"))
tbin<-mirror(read.csv("D:\\phd\\tree\\data\\refseq\\toy-k9-w100-s5-global_sim.txt"))
# temp<-merge(tall,tbin,by=c("i","j"))%>%filter(similarity.x!=similarity.y)
temp<-tbin%>%group_by(i)%>%summarise(size=length(j))%>%
  filter(size!=length(i))

temp<-check_sim(tall)


sum(tall$similarity)==sum(tbin$similarity)

idx<-c(152,397,652,870,911,1003,1020,1021,1022,1023)

t1<-mirror(read.csv("D:\\phd\\tree\\data\\refseq\\test-all_sim.txt"))
t2<-tbin%>%filter(i%in%idx & j %in%idx)
sum(t1$similarity)==sum(t2$similarity)





