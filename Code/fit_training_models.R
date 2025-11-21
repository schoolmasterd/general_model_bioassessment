#This script fits the model used to predict the presence of taxa based on the
#vegetation data and environmental variables from Coastwide Reference Monitoring System (CRMS).
#This script can take a long time to run due to the cross-validation

#load a library that allows lasso regularization
library(glmnet)

#load the bits of training data that we need
spp_pool<-read.csv("Data/training_spp_pool.csv")
env_names<-read.csv("Data/environmental_var_names.csv")
df<-read.csv("Data/CRMS_training_data.csv")

#grab spp names and site
spp_nms<-spp_pool$Short_Name

#select the rows to use 
#get the environmental variables
env_vars<-env_names$short_names[-(1:2)]
bas_vars<-c("AT","BA","BS","CS","ME","MR","PO","TE","TV")

#remove rows with any missing data
df_comp<-df[complete.cases(df[,c(env_vars,bas_vars)]),]
crms_site<-sapply(strsplit(df_comp$Station_ID,"-"),"[",1)

#select holdout cases. These were selected by stratified random sampling to get
#one site from each habitat type in each region (Chenier, Deltaic) of the Coast.
who<-c("CRMS0556","CRMS2568","CRMS0665","CRMS0146","CRMS0532","CRMS3784","CRMS6302","CRMS4690")
get_cases<-which(crms_site%in%who)

#remove holdout cases
df_train<-df_comp[-get_cases,]
dim(df_train)

#calculate weights for repeated measurement as 1/(number of obs)
wts_tab<-table(df_train$Station_ID)
wts<-1/wts_tab[df_train$Station_ID]

#center the environmental vars and save them for use with the test samples
col_means_mod_vars<-apply(df_train[,env_vars],2,mean,na.rm=T)
col_sd_mod_vars<-apply(df_train[,env_vars],2,sd,na.rm=T)

write.csv(rbind(means=col_means_mod_vars,sd=col_sd_mod_vars),"Output/enviromental_means.csv")

#replace raw env vars with normalized
df_train[,env_vars]<-sweep(sweep(df_train[,env_vars],2,col_means_mod_vars,FUN = "-"),2,col_sd_mod_vars,"/")

#look at occurrences across the training sample
occs<-apply(df_train[,spp_nms],2,sum)

#find ubiquitous spp to hold for possible inclusion later
ubiq<-names(occs[occs/dim(df_train)[1]>=.3])

####fit taxa models using lasso regularization and cross-validation to select "lambda" parameter###
#this step can take an hour or more#

fts<-sapply(spp_nms,function(x)NULL)
lambda_min<-rep(NA,length(spp_nms))
names(lambda_min)<-spp_nms

for(i in spp_nms){
  fmls<-paste0(i,paste0("~",paste0(c(env_vars,bas_vars),collapse = "+")))
  x<-model.matrix(formula(fmls),data = df_train)[,-1]
  y<-df_train[,i]
  cv.lasso<-glmnet::cv.glmnet(x,y,alpha=1,nfolds = 100,intercept=FALSE,family="binomial",
                              weights = as.vector(wts))
  lambda_min[i]<-cv.lasso$lambda.min
  fts[[i]]<-glmnet::glmnet(x,y,alpha=1,family = "binomial", lambda = cv.lasso$lambda,
                           type.logistic = "modified.Newton",maxit = 10^6,weights = as.vector(wts))
}

#save the output for assessing restoration projects
save(fts,file="Output/species_fit_objects.RData")

####Select the taxa to retain in the analysis####
#calc percent deviance
pct_deviance<-sapply(spp_nms,function(x)fts[[x]]$dev.ratio[fts[[x]]$lambda==min(fts[[x]]$lambda)])
hist(unlist(pct_deviance),xlab="Proportion Deviance Explained",main="")
abline(v=.2,lwd=2,lty=2)
sum(unlist(pct_deviance)>.10)

#grab holdout data
df_latest<-df_comp[unlist(get_cases),]
df_latest[,env_vars]<-sweep(sweep(df_latest[,env_vars],2,col_means_mod_vars,FUN = "-"),2,col_sd_mod_vars,"/")
spp_coefs<-sapply(spp_nms,function(x)as.vector(coef(fts[[x]],s=min(fts[[x]]$lambda))))
rownames(spp_coefs)<-c("Intercept",env_vars,bas_vars)

#predict occurrence across all holdout stations
n=dim(df_latest[,env_vars])[1]
foo<-matrix(NA,ncol = 4,nrow=length(spp_nms))
row.names(foo)<-spp_nms
colnames(foo)<-c("act","sim","lci","uci")

for(i in spp_nms){
  ans<-1/(1+exp(-as.matrix(cbind(rep(1,n),df_latest[,c(env_vars,bas_vars)]))%*%(as.matrix(spp_coefs[,i]))))
  bar<-rep(0,5000)
  for(j in 1:5000)bar[j]<-sum(rbinom(n = n,size = 1,prob = ans))
  foo[i,1]<-sum(df_latest[,i])
  foo[i,2:4]<-quantile(bar,c(.5,0.025,0.975))
}

#select taxa with a "accuracy ratio" of between -0.5<x<=0.5
sp_cand<-rownames(foo)[log((foo[,"sim"]+1)/(foo[,"act"]+1))<=0.5&log((foo[,"sim"]+1)/(foo[,"act"]+1))>=-0.5]

#check to see if ubiquitous taxa are included
ubiq%in%sp_cand

#look at coverage of habitats (FIBS) with candidate group
df_fibs<-read.csv("Data/FIBS_data.csv",check.names = F)

#calculate proportion of observations of each candidate spp in each habitat 
get_spp<-df_fibs$`Scientific Name As Currently Recognized`%in%spp_pool$Scientific_Name[match(sp_cand,spp_pool$Short_Name)]
fibs_table<-table(df_fibs$`Scientific Name As Currently Recognized`[get_spp],df_fibs$Community[get_spp])
pct_tab<-sweep(fibs_table,1,STATS = apply(fibs_table,1,sum),FUN = "/")[,c(2,3,1,4)]

#do the same for set of all taxa
get_spp_all<-df_fibs$`Scientific Name As Currently Recognized`%in%spp_pool$Scientific_Name[match(spp_nms,spp_pool$Short_Name)]
fibs_all_table<-table(df_fibs$`Scientific Name As Currently Recognized`[get_spp_all],df_fibs$Community[get_spp_all])
pct_tab_all<-sweep(fibs_all_table,1,STATS = apply(fibs_all_table,1,sum),FUN = "/")[,c(2,3,1,4)]

#create histograms of salinity score where Fresh = 1; Intermediate = 2; Brackish = 3; Saline = 4
grp_all<-hist(as.matrix(pct_tab_all)%*%as.matrix(1:4),freq = F)
grp_sel<-hist(as.matrix(pct_tab)%*%as.matrix(1:4))

#use chi-squared test with (number of categories -1 degree of freedom) to test if these differ
1-pchisq(sum((grp_sel$density-grp_all$density)^2/grp_all$density),length(grp_all$density)-1)

#create ans save the Figure S1
pdf("Output/figure_S1.pdf",height = 5,width = 8)
par(mfrow=c(1,2))
plot(log(foo[sp_cand,"act"]),log(foo[sp_cand,"sim"]),pch=21,bg="grey",bty="l",xlim=c(1,6),ylim=c(1,6),
     xlab="Log Observed Occurences",ylab="Log Predicted Occurences")
arrows(log(foo[sp_cand,"act"]),log(foo[sp_cand,"sim"]),log(foo[sp_cand,"act"]),log(foo[sp_cand,"lci"]),length = 0)
arrows(log(foo[sp_cand,"act"]),log(foo[sp_cand,"sim"]),log(foo[sp_cand,"act"]),log(foo[sp_cand,"uci"]),length = 0)
abline(0,1)
cor_foo<-cor(log(foo[sp_cand,"act"]+1),log(foo[sp_cand,"sim"]+1))
text(2,6,bquote(rho==.(round(cor_foo,3))))
mtext("(a)",3,adj = -.2,padj =-1)
plot(grp_all, col=rgb(0,0,1,1/4), xlim=c(1,4),freq = F,xlab="Salinity Score",main="")  # first histogram
plot(grp_sel, col=rgb(1,0,0,1/4), xlim=c(1,4),freq=F, add=T)  # second
legend("topright",legend=c("All Taxa","Selected Taxa"),fill =c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)),bty='n')
mtext("b)",3,adj = -.2,padj =-1)
dev.off()

#write the coefficients for selected taxa to a file
spp_coefs<-sapply(sp_cand,function(x)as.vector(coef(fts[[x]],s=min(fts[[x]]$lambda))))
rownames(spp_coefs)<-c("Intercept",env_vars,bas_vars)
write.csv(spp_coefs,"Output/species_models.csv")
