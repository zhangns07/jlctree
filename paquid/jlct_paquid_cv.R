library(data.table)
library(JM)
library(lcmm)
library(ipred)
library(prodlim)
library(jlctree)
library("NormPsy")
source('../simulation/util.R')

#====================
# Construct data
paquid$normMMSE <- normMMSE(paquid$MMSE)
paquid$age65 <- (paquid$age - 65) / 10
head(paquid)
paquidS <- paquid[paquid$agedem > paquid$age_init & paquid$age <= paquid$agedem, ]

paquidS2 <- data.table(paquidS)
paquidS2$age <- paquidS2[,{if(.N==1){age_init} else {c(age_init[1], age[c(1:(.N-1))])}},by=ID][,V1]

temp <- subset(paquidS2, select=c(ID, age_init, agedem, dem, male)) # baseline
temp <- unique(temp)
data <- tmerge(temp, temp, id=ID, tstart=age_init,tstop =agedem, death = event(agedem, dem)) #set range
data <- tmerge(data, paquidS2, id=ID,
               age65 = tdc(age, age65), CEP= tdc(age, CEP),
               normMMSE=tdc(age, normMMSE),
               BVRT=tdc(age, BVRT),
               IST=tdc(age, IST),
               HIER=tdc(age, HIER),
               CESD=tdc(age, CESD))

data <- data.table(data)
data <- data[!is.na(normMMSE) & !is.na(BVRT) & !is.na(IST) & !is.na(HIER) & !is.na(CESD)]

data_tmp <- data.table(data)
data_tmp <- data_tmp[,list(age65_one=rep(age65[1],.N),
                           BVRT_one=rep(BVRT[1],.N),
                           IST_one =rep(IST[1],.N),
                           HIER_one=rep(HIER[1],.N),
                           CESD_one=rep(CESD[1],.N)),by =ID]
data$age65_one <- data_tmp$age65_one
data$BVRT_one <- data_tmp$BVRT_one 
data$IST_one <-  data_tmp$IST_one 
data$HIER_one <- data_tmp$HIER_one
data$CESD_one <- data_tmp$CESD_one

#====================
# Setup cross validation: cross validate by subject
ndata <- nrow(data)
nfolds <- 10; n_per_fold <- ceiling(ndata/nfolds) * c(1:nfolds)
N_by_id <- data[,.N,by=id]
set.seed(0)
id_shuffle <- sample(seq_len(nrow(N_by_id)),replace = FALSE); N_by_id <- N_by_id[id_shuffle,]
N_by_id$cumN <- cumsum(N_by_id$N)
N_by_id$cvid<- apply(array(N_by_id$cumN),1,function(x){sum(x>n_per_fold)+1})
data <- merge(data,N_by_id[,list(id,cvid)],by='id')
cvid <- data$cvid


#====================
# Chose model
# inv: only use invariant covariates
# var: use varing covariates
# var1_var: use varing covariates but only take first obs for split; use varing for survival
# var_var1: use varing covariates for split; use varying but only take first obs for survival
# var1_var1: use varing covariates but only take first obs for both split and survival
for(treemodel in c('inv','var','var_var1')){
    for (greedy in c(TRUE)){
        for (STOP in c(2.71, 3.84, 6.63)){
            resultdir <- paste0('RESULTS_',nfolds,'folds_greedy_',greedy,'_stop_',STOP,'/paquid_tree_',treemodel,'_moretimevar')
            dir.create(resultdir,recursive=TRUE)

            #====================
            # Trees
            tree_time <- rep(0,nfolds)
            for (k in seq_len(nfolds)){
                train_data <- data[cvid!=k,]
                test_data <- data[cvid==k,]
                ntrain <- nrow(train_data)
                ntest <- nrow(test_data)

                tik <- Sys.time()
                if(treemodel == 'inv'){
                    tree <- jlctree(survival=Surv(tstart, tstop, death)~CEP+male,
                                    classmb=~CEP+male,
                                    fixed=normMMSE~CEP+poly(age65, degree=2, raw=TRUE),
                                    random=~poly(age65, degree=2, raw=TRUE),
                                    subject='id',data=train_data,
                                    control=list(cp=-Inf),
                                    parms=list(min.nevents=4,stable=TRUE,stop.thre=STOP))
                } else if (treemodel == 'var'){
                    # split: time var
                    tree <- jlctree(survival=Surv(tstart, tstop, death)~CEP+male+BVRT+IST+HIER+CESD+poly(age_init, degree=2, raw=TRUE),
                                    classmb=~CEP+male+age65+BVRT+IST+HIER+CESD,
                                    fixed=normMMSE~CEP+poly(age65, degree=2, raw=TRUE),
                                    random=~poly(age65, degree=2, raw=TRUE),
                                    subject='id',data=train_data,
                                    control=list(cp=-Inf),
                                    parms=list(min.nevents=4,stable=TRUE,stop.thre=STOP))
                } else if (treemodel == 'var1_var1'){
                    tree <- jlctree(survival=Surv(tstart, tstop, death)~CEP+male+BVRT_one+IST_one+HIER_one+
                                    CESD_one+poly(age_init, degree=2, raw=TRUE),
                                classmb=~CEP+male+age65_one+BVRT_one+IST_one+HIER_one+CESD_one,
                                fixed=normMMSE~CEP+poly(age65, degree=2, raw=TRUE),
                                random=~poly(age65, degree=2, raw=TRUE),
                                subject='id',data=train_data,
                                control=list(cp=-Inf),
                                parms=list(min.nevents=4,stable=TRUE,stop.thre=STOP))
                } else if (treemodel == 'var1_var'){
                    tree <- jlctree(survival=Surv(tstart, tstop, death)~CEP+male+BVRT+IST+HIER+CESD+poly(age_init, degree=2, raw=TRUE),
                                    classmb=~CEP+male+age65_one+BVRT_one+IST_one+HIER_one+CESD_one,
                                    fixed=normMMSE~CEP+poly(age65, degree=2, raw=TRUE),
                                    random=~poly(age65, degree=2, raw=TRUE),
                                    subject='id',data=train_data,
                                    control=list(cp=-Inf),
                                    parms=list(min.nevents=4,stable=TRUE,stop.thre=STOP))
                } else if (treemodel == 'var_var1'){
                    if(greedy){control <- rpart.control(cp=-Inf)} else {control <- rpart.control()}
                    tree <- jlctree(survival=Surv(tstart, tstop, death)~CEP+male+BVRT_one+IST_one+HIER_one+
                                    CESD_one+poly(age_init, degree=2, raw=TRUE),
                                classmb=~CEP+male+age65+BVRT+IST+HIER+CESD,
                                fixed=normMMSE~CEP+poly(age65, degree=2, raw=TRUE),
                                random=~poly(age65, degree=2, raw=TRUE),
                                subject='id',data=train_data,
                                control=list(cp=-Inf),
                                parms=list(min.nevents=4,stable=TRUE,stop.thre=STOP))
                }

                tok <- Sys.time()

                tree_time[k] <- round(difftime(tok,tik,units='secs'),4)
                filename <- paste0(resultdir,'/tree_fold_',k,'.RData')
                save(tree, file = filename)
            }
            cat(tree_time,sep=',',file=paste0(resultdir,'/tree_time'))


            #--------------------
            # Compare biomarker prediction: MSE
            PRED_tree <- rep(0,ndata)
            MSE_tree_pred <- rep(0,nfolds)
            for (k in seq_len(nfolds)){
                train_data <- data[cvid!=k,]
                test_data <- data[cvid==k,]
                ntrain <- nrow(train_data)
                ntest <- nrow(test_data)

                filename <- paste0(resultdir,'/tree_fold_',k,'.RData')
                load(filename)

                tree2 <- tree$tree
                tree2$frame[grepl('leaf',tree2$frame$var),]$yval <- sort(unique(tree2$where))
                test_class <- predict(tree2,test_data)
                test_data$node <- as.factor(test_class)

                PRED_tree[cvid==k] <- predict(tree$lmmmodel, newdata=test_data,allow.new.levels=TRUE)
                MSE_tree_pred[k] <- sqrt(mean((test_data$normMMSE - PRED_tree[cvid==k])^2))
            }

            MSE_tree <- sqrt( mean((data$normMMSE - PRED_tree)^2))
            cat(c(MSE_tree, MSE_tree_pred),sep=',',file=paste0(resultdir,'/MSE_tree'))

            #--------------------
            # Compare survival prediction
            evaltimes <- seq(min(data$tstart), max(data$tstop), length.out=100)
            IBS_tree_pred_psesub <- IBS_tree_pred <- matrix(0,nrow=nfolds,ncol=3)
            colnames(IBS_tree_pred_psesub) <- colnames(IBS_tree_pred) <-  c('diffh_diffs','diffh','diffs')
            SHAT <- array(0,c(nrow(data),100,3))

            for (k in seq_len(nfolds)){
                train_data <- data.table(data[cvid!=k,])
                test_data <- data.table(data[cvid==k,])
                ntrain <- nrow(train_data)
                ntest <- nrow(test_data)

                filename <- paste0(resultdir,'/tree_fold_',k,'.RData')
                load(filename)

                tree2 <- tree$tree
                tree2$frame[grepl('leaf',tree2$frame$var),]$yval <- sort(unique(tree2$where))
                test_class <- predict(tree2,test_data)
                test_data$node <- as.factor(test_class)

                SHAT_k <- array(0,c(nrow(test_data),100,3))
                for (j in unique(test_data$id)){
                    tmpidx <- test_data$ID == j 
                    tmpdata <- test_data[tmpidx,]
                    ntmp <- nrow(tmpdata)
                    tmpdata[ntmp,'tstop'] <- Inf
                    start.time <- tmpdata[1,tstart]

                    SHAT_k[tmpidx,,1] <- tryCatch({
                        tmp_survfit1 <- survfit(tree$coxphmodel_diffh_diffs,newdata=tmpdata,individual=TRUE)
                        tmp_survfit1$time <- tmp_survfit1$time + start.time
                        rep(getsurv(tmp_survfit1,evaltimes),each=nrow(tmpdata))
                    }, error = function(e){0})
                    SHAT_k[tmpidx,,2] <- tryCatch({
                        tmp_survfit2 <- survfit(tree$coxphmodel_diffh,newdata=tmpdata,individual=TRUE)
                        tmp_survfit2$time <- tmp_survfit2$time + start.time
                        rep(getsurv(tmp_survfit2,evaltimes),each=nrow(tmpdata))
                    }, error = function(e){0})
                    SHAT_k[tmpidx,,3] <- tryCatch({
                        tmp_survfit3 <- survfit(tree$coxphmodel_diffs,newdata=tmpdata,individual=TRUE)
                        tmp_survfit3$time <- tmp_survfit3$time + start.time
                        rep(getsurv(tmp_survfit3,evaltimes),each=nrow(tmpdata))
                    }, error = function(e){0})
                }

                SHAT[cvid==k,,] <- SHAT_k


                torm <- c(which(evaltimes < min(test_data$agedem)), which(evaltimes > max(test_data$agedem)))
                test_data$rownames <- c(1:nrow(test_data))
                rowtokeep <- test_data[,rownames[1],by=ID][,V1]

                SHAT_k <- SHAT_k[,-torm,]
                brier_surv <- Surv(test_data$agedem,test_data$dem)

                IBS_tree_pred[k,] <- apply(SHAT_k[rowtokeep,,],3,function(x){ sbrier(brier_surv[rowtokeep], t(x), evaltimes[-torm]) })
                IBS_tree_pred_psesub[k,] <- apply(SHAT_k,3,function(x){ sbrier(brier_surv, t(x), evaltimes[-torm])})
            }

            write.table(IBS_tree_pred,file=paste0(resultdir,'/IBS_tree_pred.csv'),quote=F,col.names=T,sep=",",row.names=F)
            write.table(IBS_tree_pred_psesub,file=paste0(resultdir,'/IBS_tree_pred_psesub.csv'),quote=F,col.names=T,sep=",",row.names=F)


            if(1==0){ # overall IBS
                torm <- c(which(evaltimes < min(data$agedem)), which(evaltimes > max(data$agedem)))
                data$rownames <- c(1:nrow(data))
                rowtokeep <- data[,rownames[1],by=ID][,V1]
                SHAT <- SHAT[,-torm,]
                brier_surv <- Surv(data$agedem,data$dem)

                apply(SHAT[rowtokeep,,],3,function(x){ sbrier(brier_surv[rowtokeep], t(x), evaltimes[-torm]) })
                apply(SHAT[,,],3,function(x){ sbrier(brier_surv, t(x), evaltimes[-torm]) })
            }

            # Five models:
            # 1. time-invariant survival and splits, tree and lcmm
            # 2. time-varying survival and splits but only first obs, tree and lcmm 
            # 3. time-varying survival and splits, tree 
        }
    }
}

#====================
# Read result
nfolds  <- 10; greedy  <- TRUE; 
treemodel <- 'var'
treemodel <- 'var_var1'
folder  <- paste0('RESULTS_',nfolds,'folds_greedy_',greedy,'_stop_6.63/paquid_tree_',treemodel,'_moretimevar')
folder  <- paste0('RESULTS_',nfolds,'folds_greedy_',greedy,'/paquid_tree_',treemodel,'_moretimevar')
folder  <- paste0('RESULTS_',nfolds,'folds_greedy_',greedy,'_stop_2.71/paquid_tree_',treemodel,'_moretimevar')

filename  <- paste0(folder,'/IBS_tree_pred.csv')
ret  <- fread(filename)
colMeans(ret)

filename  <- paste0(folder,'/MSE_tree')
ret  <- fread(filename)
mean(unlist(ret[,c(1:11),with=F]))

filename  <- paste0(folder,'/tree_time')
ret  <- fread(filename)
mean(unlist(ret))



# ====================
# Fit an overall model and plot
data$age  <- round(10*data$age65+65)
data$age_one  <-  round(10*data$age65_one+65)
save(data,file='data.rdata')
load('data.rdata')
# 'inv'
set.seed(0)
tree_inv <- jlctree(survival=Surv(tstart, tstop, death)~CEP+male,
                classmb=~CEP+male,
                fixed=normMMSE~CEP+poly(age65, degree=2, raw=TRUE),
                random=~poly(age65, degree=2, raw=TRUE),
                greedy=TRUE, maxng=6,
                subject='id',data=data,
                parms=list(min_nevent=2,stop_thre=3.84,stable=TRUE))
tree_inv$tree

# 'var1_var1'
tree_mid <- jlctree(survival=Surv(tstart, tstop, death)~CEP+male+
                BVRT_one+IST_one+HIER_one+CESD_one+poly(age_init, degree=2, raw=TRUE),
                classmb=~CEP+male+age_one+BVRT_one+IST_one+HIER_one+CESD_one,
                fixed=normMMSE~CEP+poly(age65, degree=2, raw=TRUE),
                random=~poly(age65, degree=2, raw=TRUE),
                greedy=TRUE,
                subject='id',data=data,
                parms=list(min_nevent=5,stop_thre=3.84,stable=TRUE))
tree_mid$tree

# 'var_var1'
tree_varsinvt <- jlctree(survival=Surv(tstart, tstop, death)~CEP+male+
                BVRT_one+IST_one+HIER_one+CESD_one+poly(age_init, degree=2, raw=TRUE),
                classmb=~CEP+male+age+BVRT+IST+HIER+CESD,
                fixed=normMMSE~CEP+poly(age65, degree=2, raw=TRUE),
                random=~poly(age65, degree=2, raw=TRUE),
                subject='id',data=data,
                control=list(cp=-Inf),
                parms=list(min.nevents=5,stop.thre=3.84,stable=TRUE))
tree_varsinvt$tree


# 'var'
tree_var <- jlctree(survival=Surv(tstart, tstop, death)~CEP+male+
                    BVRT+IST+HIER+CESD+poly(age_init, degree=2, raw=TRUE),
                classmb=~CEP+male+age+BVRT+IST+HIER+CESD,
                fixed=normMMSE~CEP+poly(age65, degree=2, raw=TRUE),
                random=~poly(age65, degree=2, raw=TRUE),
                subject='id',data=data,
                #control=list(cp--Inf)
                parms=list(min.nevents=5,stop.thre=6.63,stable=TRUE))
tree_var$tree


library(rpart.plot)
if(1==1){
pdf('paquid_tree_inv.pdf')
rpart.plot(tree_inv$tree,cex=2)
dev.off()

pdf('paquid_tree_var1_var1.pdf')
rpart.plot(tree_mid$tree,cex=1.5)
dev.off()

pdf('paquid_tree_var_var1.pdf')
rpart.plot(tree_varsinvt$tree,cex=2)
dev.off()


pdf('paquid_tree_var.pdf')
rpart.plot(tree_var$tree,cex=2)
dev.off()

}

#====================
# See whether IST is significiant in survival model
# for var_var
data$age  <- round(10*data$age65+65)
data$age_one  <-  round(10*data$age65_one+65)
tree_var <- jlctree(survival=Surv(tstart, tstop, death)~CEP+male+
                    BVRT+IST+HIER+CESD+poly(age_init, degree=2, raw=TRUE),
                classmb=~CEP+male+age+BVRT+IST+HIER+CESD,
                fixed=normMMSE~CEP+poly(age65, degree=2, raw=TRUE),
                random=~poly(age65, degree=2, raw=TRUE),
                subject='id',data=data,control=list(cp=-Inf),
                parms=list(min.nevents=5,stop.thre=3.84, stable=TRUE))


coef <- get_cox_coef(tree_var$coxphmodel_diffh_diffs)
se <- get_cox_se(tree_var$coxphmodel_diffh_diffs)
coef/se

#  2) age>=89.5 156  0.8837669  0.8837669 *
#  3) age< 89.5 1828  6.0951790  6.0951790
#      6) age>=81.5 636  2.5066020  2.5066020 *
#      7) age< 81.5 1192  2.7919420  2.7919420 *
#node4:IST                                     -1.55 0.122
#node5:IST                                     -2.28 0.023

tree_varsinvt <- jlctree(survival=Surv(tstart, tstop, death)~CEP+male+
                         BVRT_one+IST_one+HIER_one+CESD_one+poly(age_init, degree=2, raw=TRUE),
                     classmb=~CEP+male+age+BVRT+IST+HIER+CESD,
                     fixed=normMMSE~CEP+poly(age65, degree=2, raw=TRUE),
                     random=~poly(age65, degree=2, raw=TRUE),
                     subject='id',data=data,
                     control=list(cp=-Inf),
                     parms=list(min.nevents=5,stop.thre=3.84,stable=TRUE))
coef <- get_cox_coef(tree_varsinvt$coxphmodel_diffh_diffs)
se <- get_cox_se(tree_varsinvt$coxphmodel_diffh_diffs)
coef/se

#  2) age>=84.5 507  3.03677100  3.03677100 *
#  3) age< 84.5 1477  4.11212700  4.11212700
#      6) IST< 19.5 101  0.02539347  0.02539347 *
#      7) IST>=19.5 1376  2.59193800  2.59193800 *

#node4:IST_one                                 -0.77 0.443
#node5:IST_one                                 -1.80 0.071


# Conclusion: 
# for var model, IST is significant in groups where age < 89
# for vars_invt mdoel, IST_one is also siginficant in age < 84.5 groups.


