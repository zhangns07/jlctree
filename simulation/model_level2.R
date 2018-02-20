# model1: time varying, section 5 in biostats supp
# 1. y ~ betay * x + (1|id) + noise
# 2. Surv ~ coxph(0*biomarker+beta_x*X1+beta_z*Z) where the first 0 is for biomarker.
# Censor rate parameter:
# I,II, exponential, 0.5: lam_D=0.25
# I,II, exponential, 0.2: lam_D=0.07

# I,II, gompertz, 0.5: lam_D=0.55
# I,II, gompertz, 0.2: lam_D=0.15

# I,II, weibull, 0.5: lam_D=0.53
# I,II, weibull, 0.2: lam_D=0.13


library(optparse)
library(data.table)
library(plyr)
source('../0.init.R')
source('0.gen_survival.R')
option_list <- list(
                    make_option(c("-N", "--Nsub"), type="numeric", default=100, help=""),
                    make_option(c("-c", "--censor"), type="numeric", default=0, help=""),
                    make_option(c("-d", "--dist"), type="character", default='exponential', help=""),
                    make_option(c("-s", "--stop_thre"), type="numeric", default=0, help=""),
                    make_option(c("-t", "--test"), type="character", default='lrt', help="Test statistics, rsq, lrt or wald."),
                    make_option(c("-i", "--type"), type="character", default='I', help="I or II, time varying dichotomous type")
                    )

#make_option(c("-l", "--left"), type="logical", default=TRUE, help="Wether left truncated"),
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

FLAGS <- opt
Nsub <- FLAGS$Nsub
censor_rate <- FLAGS$censor
dist <- FLAGS$dist
stop_thre <- FLAGS$stop_thre
type <- FLAGS$type

sd_ranef <- 2
sd_e <- 1
Nsim <- 200
parms <- switch(dist,
                exponential=list(lambda=0.1, shape=0, beta=0.8, betaz=1.4),
                weibull=list(lambda=0.3, shape=0.8, beta=0.9, betaz=1.6),
                gompertz=list(lambda=0.2, shape=0.1, beta=1.2, betaz=2.0))

lam_D <- switch(dist,
                exponential=c(0.07,0.25), #corresponds to censor rate 0.2 and 0.5 
                weibull=c(0.13,0.53),
                gompertz=c(0.15,0.55))

survs <- survs_v3
survlist <- list(eval=surve, split=survs, init=survi)

opt2 <- FLAGS
opt2$help <- NULL
opt2$betay <- 1
opt2$model <- 'condind'
basefilename <- paste0(paste0(names(opt2),"_",opt2),collapse="_")

beta_y <- rep(opt2$betay, 5)
nsplit_RET <- rep(0, Nsim)
nsplit_prune_RET <- matrix(0,nrow=Nsim,ncol=4)
val_root_RET <- rep(0,Nsim)

for (sim in c(1:Nsim)){
    set.seed(sim)

    # X1 and time
    X1 <- as.numeric(runif(Nsub)>0.5)
    X3 <- as.numeric(runif(Nsub)>0.5)
    ntimes <- switch(type,'I'=1,'II'=3)
    times <- matrix(runif(ntimes*Nsub,min=0.6,max=6),ncol=ntimes)
    if (ntimes>1){ times <- t(apply(times,1,sort)) }
    time_T <- gen_vary_survival(Nsub, X1, times, type, dist, parms)

    if (censor_rate==0){
        time_C <- Inf
    } else{
        time_C <- rexp(Nsub, rate=ifelse(censor_rate==0.2, lam_D[1],lam_D[2]))
    }

    delta <- as.numeric(time_C >= time_T)
    time_Y <- pmin(time_T, time_C)

    # convert to pseudo object
    num_pseudo <- 1 + apply(array(seq_len(Nsub)),1,function(i){
                                sum(times[i,] < time_Y[i]) })
    all_Z <- switch(FLAGS$type,'I'=c(0,1), 'II'=c(0,1,0,1))
    LTRC_data <- ldply(array(seq_len(Nsub)),function(i){
                           num_i <- num_pseudo[i]
                           tmp_time <- pmin(c(0,times[i,],Inf),time_Y[i])[1:(num_i+1)]
                           ret <- cbind(sub=i, #y = rep(ranef[i],num_i), 
                                        X1=rep(X1[i],num_i),X2=all_Z[1:num_i],X3=rep(X3[i],num_i),
                                        time_L=tmp_time[1:num_i],time_Y=tmp_time[-1],
                                        delta=c(rep(0,num_i-1),delta[i]))})

    Nsub_pseudo <- nrow(LTRC_data)
    LTRC_data$X4 <- round(runif(Nsub_pseudo),1)
    LTRC_data$X5 <- sample(c(1:5),Nsub_pseudo,replace=TRUE)

    # biomarker
    ranef <- matrix(rnorm(Nsub*6, sd=sd_ranef),nrow=Nsub)
    ranefs <- ldply(LTRC_data$sub,function(i){ranef[i,]})
    LTRC_data$y <- rowSums(cbind(1,LTRC_data[,c('X1','X2','X3','X4','X5')]) * ranefs) 
    LTRC_data$y <- LTRC_data$y + rnorm(Nsub_pseudo, sd=sd_e)

    data <- data.frame(LTRC_data)
    cond_ind_tree <- rpart(model.frame
                           (cbind(time_L,time_Y,delta,y,X1,X2,X3,X4,X5)~X1+X2+X3+X4+X5,
                            data = data), control=rpart.control(minsplit=5),
                           method=survlist, 
                           parms=list(LTRC=1, test_stat=FLAGS$test, stop_thre=stop_thre, min_nevent=6))

    nsplit <- max(cond_ind_tree$cptable[,'nsplit'])
    nsplit_RET[sim] <- nsplit

    xfit <- xpred.rpart(cond_ind_tree, xval=5)
    xerror <- apply(xfit,2,mean)
    xstd <- apply(xfit,2,sd)/sqrt(nrow(data))
    cptable <- cbind(cond_ind_tree$cptable, xerror,xstd)

    cventry <- which.min(cptable[, "xerror"])
    xerrorcv <- cptable[cventry, "xerror"]

    for (kse in c(0:3)){
        sexerrorcv <- xerrorcv + kse*cptable[cventry, "xstd"] 
        cpcvse <- cptable[which.max(cptable[, "xerror"] <= sexerrorcv), "CP"]
        cond_ind_tree_prune <- prune(cond_ind_tree, cp=cpcvse)
        nsplit_prune <- max(cond_ind_tree_prune$cptable[,'nsplit'])
        nsplit_prune_RET[sim,(kse+1)] <- nsplit_prune
    }

    f1 <- Surv(time_L, time_Y, delta)~X1+X2+X3+X4+X5
    f2 <- Surv(time_L, time_Y, delta)~y+X1+X2+X3+X4+X5

    if (FLAGS$test == 'lrt'){
        val_root <- 2*tryCatch(coxph(f2,data)$loglik[2]- coxph(f1,data)$loglik[2], error=function(e) NA)
    } else if (FLAGS$test == 'wald'){
        val_root <- tryCatch({model <- coxph(f2,data);(model$coefficients[1])^2/ model$var[1,1]}, error=function(e) NA)
    }

    val_root_RET[sim] <- val_root

    cat('sim: ',sim,'\n')

}

for (kse in c(0:3)){
        filename <- paste0('simret_model2/',basefilename,'_s_',stop_thre,'_k_',kse,'.csv')
    write.table(cbind(nsplit_RET, nsplit_prune_RET[,(kse+1)]), file=filename, sep=',',col.names=TRUE, quote=FALSE)
}

filename <- paste0('simret_model2/',basefilename,'_root.csv')
write.table(val_root_RET,  file=filename, sep=',',col.names=TRUE, quote=FALSE)

if (1==0){
    # find the 'right' lam_D for different distribution and censor rate
    Nsub <- 1000; type='II'; dist='exponential'
    parms <- switch(dist,
                    exponential=list(lambda=0.1, shape=0, beta=0.8, betaz=1.4),
                    weibull=list(lambda=0.3, shape=0.8, beta=0.9, betaz=1.6),
                    gompertz=list(lambda=0.2, shape=0.1, beta=1.2, betaz=2.0))

    X1 <- as.numeric(runif(Nsub)>0.5)
    ntimes <- switch(type,'I'=1,'II'=3)
    times <- matrix(runif(ntimes*Nsub,min=0.6,max=6),ncol=ntimes)
    if (ntimes>1){ times <- t(apply(times,1,sort)) }
    time_T <- gen_vary_survival(Nsub, X1, times, type, dist, parms)


    sum(rexp(Nsub, rate=0.2)<time_T)/Nsub
    1==1

    # I,II, exponential, 0.5: lam_D=0.25
    # I,II, exponential, 0.2: lam_D=0.07

    # I,II, gompertz, 0.5: lam_D=0.55
    # I,II, gompertz, 0.2: lam_D=0.15

    # I,II, weibull, 0.5: lam_D=0.53
    # I,II, weibull, 0.2: lam_D=0.13

}

if (1==0){

    nsplit_RET <- rep(0,200)
    for(sim in c(1:200)){
        set.seed(sim)

        # X1 and time
        X1 <- as.numeric(runif(Nsub)>0.5)
        X3 <- as.numeric(runif(Nsub)>0.5)
        ntimes <- switch(type,'I'=1,'II'=3)
        times <- matrix(runif(ntimes*Nsub,min=0.6,max=6),ncol=ntimes)
        if (ntimes>1){ times <- t(apply(times,1,sort)) }
        time_T <- gen_vary_survival(Nsub, X1, times, type, dist, parms)

        if (censor_rate==0){
            time_C <- Inf
        } else{
            time_C <- rexp(Nsub, rate=ifelse(censor_rate==0.2, lam_D[1],lam_D[2]))
        }

        delta <- as.numeric(time_C >= time_T)
        time_Y <- pmin(time_T, time_C)

        # convert to pseudo object
        num_pseudo <- 1 + apply(array(seq_len(Nsub)),1,function(i){
                                    sum(times[i,] < time_Y[i]) })

        all_Z <- switch(FLAGS$type,'I'=c(0,1), 'II'=c(0,1,0,1))
        LTRC_data <- ldply(array(seq_len(Nsub)),function(i){
                               num_i <- num_pseudo[i]
                               tmp_time <- pmin(c(0,times[i,],Inf),time_Y[i])[1:(num_i+1)]
                               ret <- cbind(sub=i, #y = rep(ranef[i],num_i), 
                                            X1=rep(X1[i],num_i),X2=all_Z[1:num_i],X3=rep(X3[i],num_i),
                                            time_L=tmp_time[1:num_i],time_Y=tmp_time[-1],
                                            delta=c(rep(0,num_i-1),delta[i]))})

        Nsub_pseudo <- nrow(LTRC_data)
        LTRC_data$X4 <- round(runif(Nsub_pseudo),1)
        LTRC_data$X5 <- sample(c(1:5),Nsub_pseudo,replace=TRUE)

        # biomarker
        #    ranef <- matrix(rnorm(Nsub*6, sd=sd_ranef),nrow=Nsub)
        #    ranefs <- ldply(LTRC_data$sub,function(i){ranef[i,]})
        #    LTRC_data$y <- rowSums(cbind(1,LTRC_data[,c('X1','X2','X3','X4','X5')]) * ranefs) 
        #    LTRC_data$y <- LTRC_data$y + rnorm(Nsub_pseudo, sd=sd_e)
        #    data <- copy(LTRC_data)
        ranef <- matrix(rnorm(Nsub,sd=0.2),ncol=1)
        ranefs <- unlist(ldply(LTRC_data$sub,function(i){ranef[i,]}))
        LTRC_data$y <- rowSums(LTRC_data[,c('X1','X2')] ) + ranefs
        LTRC_data$y <- LTRC_data$y + rnorm(Nsub_pseudo, sd=0.1)
        data <- copy(LTRC_data)

        cond_ind_tree <- rpart(model.frame
                               (cbind(time_L,time_Y,delta,y,X3,X4,X5)~X1+X2+X3+X4+X5,
                                data = data), control=rpart.control(minsplit=5),
                               method=survlist, 
                               parms=list(LTRC=1, test_stat='lrt', stop_thre=3.84, min_nevent=6))

        nsplit <- max(cond_ind_tree$cptable[,'nsplit'])
        nsplit_RET[sim] <- nsplit
        cat('sim: ',sim,'\n')

    }


    f1 <- Surv(time_L, time_Y, delta)~X3+X4+X5
    f2 <- Surv(time_L, time_Y, delta)~y+X3+X4+X5
    PARMS <- expand.grid(c(0,1,NA),c(0,1,NA))
    RET <- c(0,9)
    for(i in seq_len(nrow(PARMS))){
        subdata <- data
        if (!is.na(PARMS[i,1])){
            subdata <-subset(subdata,X1==PARMS[i,1]) 
        }
        if (!is.na(PARMS[i,2])){
            subdata <-subset(subdata,X2==PARMS[i,2]) 
        }
        RET[i] <- 2*(coxph(f2,subdata)$loglik[2]- coxph(f1,subdata)$loglik[2])
    }
    cbind(PARMS,RET)


    data <- copy(LTRC_data)
    data$z <- 1
    f1 <- Surv(time_L, time_Y, delta)~X3+X4+X5
    f2 <- Surv(time_L, time_Y, delta)~y+X3+X4+X5
    f3 <- Surv(time_L, time_Y, delta)~z+X3+X4+X5
    coxph(f1,subset(data,X1==0&X2==0))
    coxph(f2,subset(data,X1==0&X2==0))
    coxph(f3,subset(data,X1==0&X2==0))

    coxph(f1,data)$loglik
    coxph(f2,data)$loglik
    coxph(f3,data)$loglik

    m1 <- coxph(f1,data)
}
