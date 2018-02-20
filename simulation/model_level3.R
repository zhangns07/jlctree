# model3: non-null 
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
                    make_option(c("-t", "--test"), type="character", default='lrt', help="Test statistics, rsq, lrt or wald."))

#make_option(c("-i", "--type"), type="character", default='I', help="I or II, time varying dichotomous type")
#make_option(c("-l", "--left"), type="logical", default=TRUE, help="Wether left truncated"),
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

FLAGS <- opt
Nsub <- FLAGS$Nsub
censor_rate <- FLAGS$censor
dist <- FLAGS$dist
stop_thre <- FLAGS$stop_thre

sd_ranef <- 0.2
#Nsim <- 200
Nsim <- 5
beta_y <- rep(0,5)
parms <- switch(dist,
                exponential=list(lambda=0.1),
                weibulld=list(alp=0.9, beta=1),
                weibulli=list(alp=3, beta=2))


# By average ratios
slopes <- switch(dist,
                 exponential=matrix(c(0,0,0,0.56,0.56,0.09,0.92,0.92,0.15,1.46,1.46,0.24),ncol=3,byrow=TRUE),
                 weibulld=matrix(c(-1.17,-1.17,-0.19,-0.66,-0.66,-0.11,-0.55,-0.55,-0.09,0,0,0),ncol=3,byrow=TRUE),
                 weibulli=matrix(c(-3.22,-3.22,-0.54,-2.26,-2.26,-0.38, -1.53,-1.53,-0.26, 0,0,0),ncol=3,byrow=TRUE))

lam_D <- switch(dist,
                exponential=list('0.2'=c(0.026,0.048,0.086,0.156),
                                 '0.5'=c(0.114,0.232,0.396,0.74)),
                weibulld=list('0.2'=c(0.024,0.072,0.084,0.238),
                              '0.5'=c(0.13,0.318,0.384,0.956)),
                weibulli=list('0.2'=c(0.024,0.04,0.076,0.19),
                              '0.5'=c(0.084,0.168,0.224,0.64)))

survs <- survs_v3
survlist <- list(eval=surve, split=survs, init=survi)

opt2 <- FLAGS
opt2$help <- NULL
opt2$model <- 'model3'
basefilename <- paste0(paste0(names(opt2),"_",opt2),collapse="_")

nsplit_RET <- rep(0, Nsim)
nsplit_prune_RET <- matrix(0,nrow=Nsim,ncol=4)
val_root_RET <- rep(0,Nsim)

RET <- matrix(0,ncol=7,nrow=Nsim*5)
colnames(RET) <- c('sim','k','nsplit','ISE','IBS_pred','MSE_b','MSE_beta')
RET_iter <- 1

for (sim in c(1:Nsim)){
    set.seed(sim)

    # X1 - X5
    X1 <- as.numeric(runif(2*Nsub)>0.5)
    X2 <- as.numeric(runif(2*Nsub)>0.5)
    X3 <- as.numeric(runif(2*Nsub)>0.5)
    X4 <- round(runif(2*Nsub),1)
    X5 <- sample(c(1:5),2*Nsub,replace=TRUE)
    X <- cbind(X1,X2,X3,X4,X5)
    g <- 2*(X1)+X2+1

    ebx <- exp(rowSums(slopes[g,] * cbind(X3,X4,X5)))
    time_T <- gen_model3_survival(ebx, dist, parms)
    time_L <- runif(2*Nsub, min=0, max=1)

    time_tokeep <- time_L < time_T
    time_L <- time_L[time_tokeep][1:Nsub]
    time_T <- time_T[time_tokeep][1:Nsub]
    X <- X[time_tokeep,][1:Nsub,]
    g <- g[time_tokeep][1:Nsub]
    ebx <- ebx[time_tokeep][1:Nsub]

    if (censor_rate==0){ 
        time_C <- Inf
    } else{
        time_C <- time_L + rexp(Nsub,lam_D[[censor_rate]][g])
    }

    delta <- as.numeric(time_C >= time_T)
    time_Y <- pmin(time_T, time_C)


#    # biomarker
#    y <- X[,'X1']+X[,'X2'] + rnorm(Nsub, sd=sd_ranef) 
#    data <- data.frame(cbind(time_L, time_Y, delta, X, y))
    num_measure <-  1+rpois(Nsub,lambda=1)
    LTRC_data <- ldply(array(seq_len(Nsub)),function(i){
                           tmp_time <- c(time_L[i],sort(runif(num_measure[i], min=time_L[i],max=time_Y[i])),time_Y[i])
                           num_i <- num_measure[i]+1
                           ret <- cbind(ID=i, X1=rep(X1[i],num_i),
                                        X2=rep(X2[i],num_i), X3=rep(X3[i],num_i),
                                        X4=rep(X4[i],num_i), X5=rep(X5[i],num_i),
                                        time_L=tmp_time[1:num_i],time_Y=tmp_time[-1],
                                        delta=c(rep(0,num_i-1),delta[i])) })

    # biomarker
    ranef <- rnorm(Nsub, sd=sd_ranef)
    ranefs <- ranef[LTRC_data$ID]
    y <- LTRC_data$X1 + LTRC_data$X2 + ranefs + rnorm(nrow(LTRC_data), sd=sd_e) 
    data <- cbind(LTRC_data,y)

    cond_ind_tree <- rpart(model.frame
                           (cbind(time_L,time_Y,delta,y,X3,X4,X5)~X1+X2+X3+X4+X5,
                            data = data), control=rpart.control(minsplit=5),
                           method=survlist, 
                           parms=list(LTRC=1, test_stat=FLAGS$test, stop_thre=stop_thre, min_nevent=6))

    nsplit <- max(cond_ind_tree$cptable[,'nsplit'])
    RET[RET_iter,] <- c(sim,k=-Inf, nsplit,eval_tree_pred(data,dist, slopes, parms, cond_ind_tree$where))
    RET_iter <- RET_iter+1

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
        if (nsplit_prune == nsplit){
            RET[RET_iter,] <- RET[RET_iter-1,]
            RET[RET_iter,2] <- kse
        } else {
            RET[RET_iter,] <- c(sim,k=kse, nsplit_prune, eval_tree_pred(data,dist, slopes, parms, cond_ind_tree_prune$where))
        }
        RET_iter <- RET_iter+1
    }

    cat('sim: ',sim,'\n')
}


filename <- paste0('simret_model3/',basefilename,'.csv')
write.table(RET, file=filename, sep=',',col.names=TRUE, quote=FALSE)

if(1==0){

    dist='weibulli'
    Nsub=2000
    slopes <- switch(dist,
                     exponential=matrix(c(0,0,0,0.56,0.56,0.09,0.92,0.92,0.15,1.46,1.46,0.24),ncol=3,byrow=TRUE),
                     weibulld=matrix(c(-1.17,-1.17,-0.19,-0.66,-0.66,-0.11,-0.55,-0.55,-0.09,0,0,0),ncol=3,byrow=TRUE),
                     weibulli=matrix(c(-3.22,-3.22,-0.54,-2.26,-2.26,-0.38, -1.53,-1.53,-0.26, 0,0,0),ncol=3,byrow=TRUE))


    parms <- switch(dist,
                    exponential=list(lambda=0.1),
                    weibulld=list(alp=0.9, beta=1),
                    weibulli=list(alp=3, beta=2),
                    bathtub=list(a=0.01,b=1,c=5))


    # X1 - X5
    X1 <- as.numeric(runif(2*Nsub)>0.5)
    X2 <- as.numeric(runif(2*Nsub)>0.5)
    X3 <- as.numeric(runif(2*Nsub)>0.5)
    X4 <- round(runif(2*Nsub),1)
    X5 <- sample(c(1:5),2*Nsub,replace=TRUE)
    X <- cbind(X1,X2,X3,X4,X5)
    g <- 2*(X1)+X2+1

    ebx <- exp(rowSums(slopes[g,] * cbind(X3,X4,X5)))
    time_T <- gen_model3_survival(ebx, dist, parms)
    time_L <- runif(2*Nsub, min=0, max=1)

    time_tokeep <- time_L < time_T
    time_L <- time_L[time_tokeep][1:Nsub]
    time_T <- time_T[time_tokeep][1:Nsub]
    X <- X[time_tokeep,][1:Nsub,]
    g <- g[time_tokeep][1:Nsub]
    ebx <- ebx[time_tokeep][1:Nsub]

    RET <- matrix(0,nrow=2,ncol=4)
    for (G in c(1:4)){
        tmpidx <- g==G
        lams <- c(1:500)/500
        rates <- unlist(lapply(lams,function(x){ sum(time_L[tmpidx] + rexp(sum(tmpidx),x) < time_T[tmpidx])/sum(tmpidx)}))
        RET[,G] <- c(max(lams[which(rates < 0.2)]), max(lams[which(rates < 0.5)]))
    }

    cat(paste0(RET[1,],collapse=','),'\n')
    cat(paste0(RET[2,],collapse=','),'\n')

#----------
# Figure out what censor parameters to use.
lam_D <- switch(dist,
                exponential=list('0.2'=c(0.026,0.048,0.086,0.156),
                                 '0.5'=c(0.114,0.232,0.396,0.74)),
                weibulld=list('0.2'=c(0.024,0.072,0.084,0.238),
                              '0.5'=c(0.13,0.318,0.384,0.956)),
                weibulli=list('0.2'=c(0.024,0.04,0.076,0.19),
                              '0.5'=c(0.084,0.168,0.224,0.64)))


# X1 - X5
X1 <- as.numeric(runif(2*Nsub)>0.5)
X2 <- as.numeric(runif(2*Nsub)>0.5)
X3 <- as.numeric(runif(2*Nsub)>0.5)
X4 <- round(runif(2*Nsub),1)
X5 <- sample(c(1:5),2*Nsub,replace=TRUE)
X <- cbind(X1,X2,X3,X4,X5)
g <- 2*(X1)+X2+1

Nsub=2000

dist='weibulld'
slopes <- matrix(c(-1.17,-1.17,-0.19,-0.66,-0.66,-0.11,-0.55,-0.55,-0.09,0,0,0),ncol=3,byrow=TRUE)
parms <- list(alp=0.9, beta=1)
ebx <- exp(rowSums(slopes[g,] * cbind(X3,X4,X5)))
time_WD <- gen_model3_survival(ebx, dist, parms)

dist='exponential'
slopes <- matrix(c(0,0,0,0.56,0.56,0.09,0.92,0.92,0.15,1.46,1.46,0.24),ncol=3,byrow=TRUE)
parms <- list(lambda=0.1)
ebx <- exp(rowSums(slopes[g,] * cbind(X3,X4,X5)))
time_EXP <- gen_model3_survival(ebx, dist, parms)

dist='weibulli'
slopes <- matrix(c(-3.22,-3.22,-0.54,-2.26,-2.26,-0.38, -1.53,-1.53,-0.26, 0,0,0),ncol=3,byrow=TRUE)
parms <- list(alp=3, beta=2)
ebx <- exp(rowSums(slopes[g,] * cbind(X3,X4,X5)))
time_WI <- gen_model3_survival(ebx, dist, parms)


node_idx <- cond_ind_tree$where

eval_tree_pred<- function
(data,dist, slopes, parms, node_idx){
    Nsub_pseudo  <- nrow(data)
    g <- 2*(data$X1)+data$X2+1
    ebx <- exp(rowSums(slopes[g,] * cbind(data$X3,data$X4,data$X5)))

    # ---------- Survival Prediction 
    # get true survival prob
    pred_parms <- matrix(0,nrow=Nsub,ncol=3)
    true_parms <- slopes[g,]

    # get predicted survival prob
    d_idx <- node_idx
    formula <-Surv(time_L,time_Y,delta) ~ X3+X4+X5

    times_order <- order(data$time_Y)
    times <- data$time_Y[times_order]
    ntimes <- length(times)
    ISE <- 0
#    SURVS <- matrix(0,nrow=ntimes,ncol=nrow(data))

    for (i in unique(d_idx)){
        subset_id <- d_idx== i
        tmpdata <- data[subset_id,]
        bo<-0
        while(bo!=10){
            mod  <- try(coxph(formula, tmpdata),silent=TRUE)
            if (class(mod)=="try-error") {
                bo <- bo+1
                tmpdata <- tmpdata[sample(c(1:nrow(tmpdata)),replace = TRUE),]
            } else break 
        }

        for (x in c(1:nrow(tmpdata))){
            Shat <- getsurv(survfit(mod,newdata=tmpdata[x,]),times)
#            org_idx <- seq_along(subset_id)[subset_id][x]
#            SURVS[,org_idx] <- Shat

            tmpg <- 2*(tmpdata[x,'X1'])+tmpdata[x,'X2']+1
            tmpebx <- exp(sum(slopes[tmpg,] * tmpdata[x,c('X3','X4','X5')]))
            if (dist=='exponential'){
                Strue <- exp(-tmpebx*times*parms$lambda)
            } else if (dist=='weibulld' | dist=='weibulli'){
                Strue <- exp(-(times/parms$beta)^(parms$alp) * tmpebx)
            }

            scores <- (Shat - Strue)^2
            SE <- sum(0.5*(scores[-1]+scores[-ntimes]) * diff(times)) / diff(range(times))
            ISE <- ISE + SE
        }

        pred_parms[subset_id,] <- rep(mod$coefficients,each=sum(subset_id))
    }

#    surv <- Surv(data$time_Y, data$delta)
#    IBS_pred <- sbrier(surv, SURVS)[1,1]
    ISE <- ISE/Nsub

    pred_parms[is.na(pred_parms)] <- 0
    MSE_b <- mean(rowSums((true_parms - pred_parms)^2))

    # ---------- Biomarker Prediction
    ymod <- lmer(y ~ d_idx + (1|ID),data=data)
    pred_y<- predict(ymod)
    MSE_y <- mean((pred_y - data$y)^2)

    return(c(ISE=ISE,MSE_b=MSE_b,MSE_y=MSE_y))
}
}
