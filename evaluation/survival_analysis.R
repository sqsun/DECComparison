###########################
## survival analysis
## an example


## once you have cell proportion matrix, you can clustering the mixture sample into different group 
## using k-means clustering methods or other clustering methods, such as hierarchical clustering(hclust)

num_clust <- 3 # for example clustering 3 groups
group <- kmeans(cell_composition, num_clust, nstart=10000, iter.max=10000000)$cluster

##survival
library(survival)
s <- Surv(time, censor)
chisq_val <- survdiff(s ~ group)$chisq
logrank <- 1 - pchisq(chisq_val, num_clust-1)



#############################################
## survival analysis example 
load(file = 'example_survival_group.rdata')
head(cxb6_expr)
load(file = 'example_meta.rdata')
# 0=right censored, 1=event at time, 2=left censored, 3=interval censored.
meta=meta[meta$month<120,]
head(meta)

dat=merge(cxb6_expr,meta,by='gsm')
head(dat)
table(dat$study)


library(survival)
library(survminer)
#dat=dat[dat$study==1,]
sfit <- survfit(Surv(month, event)~group, data=dat)
sfit
summary(sfit)
ggsurvplot(
  sfit, risk.table = TRUE, ggtheme = theme_bw(),
  pval = TRUE, pval.coord = c(0, 0.03)
)

ggsurvplot(sfit, conf.int=F, pval=TRUE)

ggsurvplot(sfit,palette = c("#E7B800", "#2E9FDF"),
           risk.table =TRUE,pval =TRUE,
           conf.int =TRUE,xlab ="Time in months", 
           ggtheme =theme_light(), 
           ncensor.plot = TRUE)
ggsurvplot(
  sfit, risk.table = TRUE, ggtheme = theme_bw(),
  pval = TRUE, pval.coord = c(0, 0.03)
) 

gse_list=unique(dat$study)

i=gse_list[13]
sub_dat=dat[dat$study==i,] 
sfit <- survfit(Surv(month, event)~group, data=sub_dat)
sfit
summary(sfit)
file=paste0('survival_study_',sub_dat[1,6],'.pdf')
p=ggsurvplot(sfit, conf.int=F, pval=TRUE)

p
