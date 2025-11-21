#This script uses the combination of trained model and data from the Restoration
# areas to calculate the MVC

#load the library for latex-style sybmols in plots
library(latex2exp)

#load the data from the training procedure
spp_coef<-read.csv("Output/species_models.csv",row.names = 1)
spp_info<-read.csv("Data/training_spp_pool.csv")
env_adj<-read.csv("Output/enviromental_means.csv",row.names = 1)

#grab the variable for basin
bas_vars<-rownames(spp_coef)[8:16]
df_test<-read.csv("Data/Other_rest_test_data.csv")

#grab the set of species names that are both in the restoration assemblages
#and that we have models for
spp_nms<-names(spp_coef)

#look species absent from restoration assemblage data and add them and set occupancy to zero
missing_spp<-spp_nms[!spp_nms%in%names(df_test)]
miss_mat<-matrix(0,nrow=dim(df_test)[1],ncol=length(missing_spp))
colnames(miss_mat)<-missing_spp
df_test<-data.frame(df_test,miss_mat)

#take a peek
spp_coef[,spp_nms]
df_test[,spp_nms]

#get complete cases for the data
df_no_na<-df_test[complete.cases(df_test),]
env_vars<-names(env_adj)

#set up function to calculate the quantile residual
res_fun<-function(x)(pbinom(x-1,size = 1,spp_prob)+.5*dbinom(x,size = 1,spp_prob))

#theoretical mean and variance of the quantile residual
fun<-function(x,p)qlogis(pbinom(x-1,size = 1,p)+.5*dbinom(x,size = 1,p))
qr_mu<-function(p)abs(fun(1,p))*p+abs(fun(0,p))*(1-p)
qr_var<-function(p)(abs(fun(1,p))-qr_mu(p))^2*p+(abs(fun(0,p))-qr_mu(p))^2*(1-p)

#grab the holdout observations
latest_obs<-tapply(df_no_na$Year,df_no_na$Station_ID,function(x)max(x))
get_cases<-sapply(unique(df_no_na$Station_ID),function(x)which(df_no_na$Station_ID==x&df_no_na$Year==latest_obs[x]))
df_obs<-df_no_na[get_cases,]

#normalize the environmental variables using mean and sd from the training set 
env_means<-as.numeric(env_adj[1,])
env_sd<-as.numeric(env_adj[2,])
names(env_means)<-names(env_sd)<-colnames(env_adj)

#normalize variable 
df_obs[,env_vars]<-sweep(sweep(df_obs[,env_vars],2,env_means,FUN = "-"),2,env_sd,"/")

#set up the variables to hold the output for recording the QR and Upsilon
len<-dim(df_obs)[1]
test_metric<-rep(NA,len)
spp_rich_dist<-rep(NA,len)
names(test_metric)<-df_obs$Station_ID
qres<-matrix(NA,nrow=len,ncol=length(spp_nms))
colnames(qres)<-spp_nms
rownames(qres)<-df_obs$Station_ID

#check to make sure predictors are in the data
env_vars%in%colnames(df_obs)
bas_vars%in%colnames(df_obs)

for(i in 1:len){
  spp_prob <- sapply(spp_nms, function(x)
    1 / (1 + exp(-unlist(c(1, df_obs[i, c(env_vars,bas_vars)])) %*% spp_coef[c("Intercept",env_vars,bas_vars), x])))
  spp_prob[spp_prob<1e-10]<-1e-10
  spp_prob[spp_prob==1]<-(1-1e-10)
  qres[i,spp_nms]<-qlogis(res_fun(unlist(df_obs[i, spp_nms])))
  test_ans <- sum(abs(qres[i,spp_nms]))
  test_metric[i]<-(test_ans-sum(qr_mu(spp_prob)))/sqrt(sum(qr_var(spp_prob)))
  spp_rich_dist[i]<-(sum(df_obs[i, spp_nms])-sum(spp_prob))/sqrt(sum(spp_prob*(1-spp_prob)))
}

plt_details<-sapply(names(test_metric),function(x)NULL)

for(i in 1:length(test_metric)){
  tmp<-qres[i,order(qres[i,])]
  plt_details[[i]]<-list(year=df_obs$Year[i],test_metric=round(test_metric[i],3),spp_rich_metric=round(spp_rich_dist[i],3),spp_present=
                           spp_nms[which(df_obs[i,spp_nms]==1)],
                         surprising_absences=tmp[which(tmp<(-2.944439))], surprising_presences=tmp[which(tmp>(2.944439))])
}

# create results graphs for manuscript

pa_a<-grep("PA-A",names(test_metric))
pa_b<-grep("PA-B",names(test_metric))
pa_c<-grep("PA-C",names(test_metric))
proj_ind<-rep(0,length(test_metric))
proj_ind[pa_a]<-1
proj_ind[pa_b]<-2
proj_ind[pa_c]<-3

#test the correlation between richness and MVC
cor.test(spp_rich_dist,test_metric)

#store the restoration site mean and sd
rest_site_mean<-c(mean(test_metric[pa_a]),mean(test_metric[pa_b]),mean(test_metric[pa_c]))
rest_site_sd<-c(sd(test_metric[pa_a]),sd(test_metric[pa_b]),sd(test_metric[pa_c]))

#create and save Figure 4 of main text
clrs<-c("lightgrey", "#708090", "black")
pdf("Output/figure_4.pdf",width = 10,height = 6)
par(mfrow=c(1,2),mai=c(.8,.5,.2,.2),oma=c(1.5,1.5,1.5,1.5))
plot(1:3,rest_site_mean,
     pch=21,bg=clrs,cex=1.25,bty="l",ylim = c(-1,8),xlim=c(0,4),
     ylab=latex2exp::TeX(r'(Distance from Reference $(\Upsilon)$)'),xlab = "Restoration Project",xaxt='n',cex.lab=1.25)
arrows(1:3,rest_site_mean,1:3,(rest_site_mean+rest_site_sd),length = 0,lwd=1.5)
arrows(1:3,rest_site_mean,1:3,(rest_site_mean-rest_site_sd),length = 0,lwd=1.5)
points(1:3,rest_site_mean,  pch=21,bg=clrs,cex=1.25)
axis(side=1,at = 1:3,labels =c("PA-A","PA-B","PA-C"),cex.lab=1.25)
mtext("(a)",side = 3,adj=-.15)
plot(spp_rich_dist,test_metric,pch=21,bg=clrs[proj_ind],
     ylab="",
     xlab = "Transformed Species Richness",bty="l",cex.lab=1.25,cex=1.25)
legend("bottomright",legend = c("PA-A","PA-B","PA-C"),pch=21,pt.bg = clrs,bty="n")
mtext("(b)",side = 3,adj=-.15)
mtext(latex2exp::TeX(r'(Distance from Reference $(\Upsilon)$)'),side = 2,outer = T,cex=1.25)
dev.off()

#collect unexpected taxa by project area
pa_a_ab<-table(unlist(sapply(plt_details[pa_a],function(x)names(x$surprising_absences))))/length(pa_a)
pa_b_ab<-table(unlist(sapply(plt_details[pa_b],function(x)names(x$surprising_absences))))/length(pa_b)
pa_c_ab<-table(unlist(sapply(plt_details[pa_c],function(x)names(x$surprising_absences))))/length(pa_c)
pa_c_ab<-c(0)
pa_a_pr<-table(unlist(sapply(plt_details[pa_a],function(x)names(x$surprising_presences))))/length(pa_a)
pa_b_pr<-table(unlist(sapply(plt_details[pa_b],function(x)names(x$surprising_presences))))/length(pa_b)
pa_c_pr<-table(unlist(sapply(plt_details[pa_c],function(x)names(x$surprising_presences))))/length(pa_c)

#common surprisingly present taxa
c(names(pa_a_pr),names(pa_b_pr))[duplicated(c(names(pa_a_pr),names(pa_b_pr)))]

#create and save Figure 5 of main text
pdf(file = "Output/figure_5.pdf",height = 8,width = 10)
par(mar=c(5,8,2,2),mfrow=c(1,3),oma=c(2,4,1,1))
barplot(pa_a_pr[order(pa_a_pr,decreasing = T)],horiz = T,
        names.arg =gsub(pattern = "_"," ",names(pa_a_pr[order(pa_a_pr,decreasing = T)])),
        las=1,cex.names = 1,xlab="Proportion of Stations",xlim=c(0,1),cex.axis = 1.25,
        cex.lab=1.25,col="black")
mtext("(a)",side = 3,adj = -.5)
mtext("Project Area A",side = 3)

barplot(pa_b_pr[order(pa_b_pr,decreasing = T)],horiz = T,
        names.arg =gsub(pattern = "_"," ",names(pa_b_pr[order(pa_b_pr,decreasing = T)])),
        las=1,cex.names = 1,xlab="Proportion of Stations",xlim=c(0,1),cex.axis = 1.25,
        cex.lab=1.25,col="black")
mtext("(b)",side = 3,adj = -.5)
mtext("Project Area B",side = 3)
barplot(pa_c_pr[order(pa_c_pr,decreasing = T)],horiz = T,
        names.arg =gsub(pattern = "_"," ",names(pa_c_pr[order(pa_c_pr,decreasing = T)])),
        las=1,cex.names = 1,xlab="Proportion of Stations",xlim=c(0,.5),cex.axis = 1.25,
        cex.lab=1.25,col="black")
mtext("Project Area C",side = 3)
mtext("(c)",side = 3,adj = -.5)
dev.off()

