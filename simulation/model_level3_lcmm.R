# model3: non-null 
library(optparse)
library(data.table)
library(plyr)
library(lcmm)
source('../0.init.R')
source('0.gen_survival.R')
option_list <- list(
                    make_option(c("-N", "--Nsub"), type="numeric", default=100, help=""),
                    make_option(c("-c", "--censor"), type="numeric", default=0, help=""),
                    make_option(c("-d", "--dist"), type="character", default='exponential', help=""))

#make_option(c("-t", "--test"), type="character", default='lrt', help="Test statistics, rsq, lrt or wald."))
#make_option(c("-s", "--stop_thre"), type="numeric", default=0, help=""),
#make_option(c("-i", "--type"), type="character", default='I', help="I or II, time varying dichotomous type")
#make_option(c("-l", "--left"), type="logical", default=TRUE, help="Wether left truncated"),
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

FLAGS <- opt
Nsub <- FLAGS$Nsub
censor_rate <- FLAGS$censor
dist <- FLAGS$dist


sd_ranef <- 0.2
sd_e <- 0.1
Nsim <- 200
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

opt2 <- FLAGS
opt2$help <- NULL
opt2$model <- 'model3'
opt2$alg<-'lcmm'
basefilename <- paste0(paste0(names(opt2),"_",opt2),collapse="_")

RET <- matrix(0,ncol=4,nrow=Nsim)
colnames(RET) <- c('sim','ISE','MSE_b','MSE_beta')
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

    # biomarker
    y <- X[,'X1']+X[,'X2'] + rnorm(Nsub, sd=sd_ranef) 
    data <- data.frame(cbind(time_L, time_Y, delta, X, y))
    data$ID <- as.numeric(rownames(data))

    m1 <- Jointlcmm(fixed=y~X1+X2+X3+X4+X5,
                    subject='ID',
                    survival = Surv(time_L,time_Y,delta)~X3+X4+X5,
                    hazard="Weibull",hazardtype="Specific",ng=1,data=data)


    m4 <- Jointlcmm(fixed=y~X1+X2+X3+X4+X5,mixture=~1,
                    classmb=~X1+X2+X3+X4+X5,subject='ID',
                    survival = Surv(time_L,time_Y,delta)~mixture(X3+X4+X5),
                    hazard="Weibull",hazardtype="Specific",ng=4,data=data,B=m1)

    if (m4$conv %in% c(4,5)){
        m4 <- Jointlcmm(fixed=y~X1+X2+X3+X4+X5,mixture=~1,
                        classmb=~X1+X2+X3+X4+X5,subject='ID',
                        survival = Surv(time_L,time_Y,delta)~mixture(X3+X4+X5),
                        hazard="Weibull",hazardtype="Specific",ng=4,data=data,maxiter=30)
    }

    if (m4$conv %in% c(4,5)){
        m4 <- Jointlcmm(fixed=y~X1+X2+X3+X4+X5,mixture=~1,
                        classmb=~X1+X2+X3+X4+X5,subject='ID',
                        survival = Surv(time_L,time_Y,delta)~mixture(X3+X4+X5),
                        hazard="Weibull",hazardtype="Specific",ng=4,data=data,B=random(m1),maxiter=30)
    }

    if (m4$conv %in% c(4,5)){
        RET[RET_iter,] <- c(sim, 0, 0, 0)
    } else {
        RET[RET_iter,] <- c(sim, eval_lcmm_pred(data,dist,slopes,parms,m4))
    }
    RET_iter <- RET_iter+1

    cat('sim: ',sim,'\n')

    if (sim %% 10 ==0){
        filename <- paste0('simret_model3/',basefilename,'.csv')
        write.table(RET, file=filename, sep=',',col.names=TRUE, quote=FALSE)
    }
}

filename <- paste0('simret_model3/',basefilename,'.csv')
write.table(RET, file=filename, sep=',',col.names=TRUE, quote=FALSE)


if(1==0){
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

    m1 <- Jointlcmm(fixed=y~X1+X2+X3+X4+X5,
                    subject='ID',
                    survival = Surv(time_L,time_Y,delta)~X3+X4+X5,
                    hazard="Weibull",hazardtype="Specific",ng=1,data=data)


    m4 <- Jointlcmm(fixed=y~X1+X2+X3+X4+X5,mixture=~1,
                    classmb=~X1+X2+X3+X4+X5,subject='ID',
                    survival = Surv(time_L,time_Y,delta)~mixture(X3+X4+X5),
                    hazard="Weibull",hazardtype="Specific",ng=4,data=data,B=m1)

    eval_lcmm_pred(data,dist,slopes,parms,m4)

eval_lcmm_pred<- function
(data, dist, slopes, parms, mod){
    Nsub_pseudo<- nrow(data)
    g <- 2*(data$X1)+data$X2+1
    true_parms <- slopes[g,]
    ebx <- exp(rowSums(slopes[g,] * cbind(data$X3,data$X4,data$X5)))

    #  coeff
    coefs <- mod$best
    coefstart <- 27
    pred_slopes <- matrix(coefs[coefstart:(coefstart+11)],nrow=4,ncol=3)

    predclass <- (mod$pprob$class)[data$ID]
    pred_parms <- pred_slopes[predclass,]
    MSE_b <- mean(rowSums((true_parms - pred_parms)^2))

    # Shat
    times <- mod$predSurv[,1]
    ntimes <- length(times)
    ISE <- 0

    for (x in c(1:nrow(data))){
        tmpclass <- predclass[x]
        tmpebx <- exp(sum(pred_slopes[tmpclass,] * data[x,c('X3','X4','X5')]))
        Shat <- exp(-tmpebx*mod$predSurv[,paste0('event1.CumRiskFct',tmpclass)])

        tmpg <- g[x]
        tmpebx <- exp(sum(slopes[tmpg,] * data[x,c('X3','X4','X5')]))
        if (dist=='exponential'){
            Strue <- exp(-tmpebx*times*parms$lambda)
        } else if (dist=='weibulld' | dist=='weibulli'){
            Strue <- exp(-(times/parms$beta)^(parms$alp) * tmpebx)
        }

        scores <- (Shat - Strue)^2
        SE <- sum(0.5*(scores[-1]+scores[-ntimes]) * diff(times)) / diff(range(times))
        ISE <- ISE + SE
    }

    ISE <- ISE/Nsub_pseudo

    # ---------- Biomarker Prediction
    pred_y <- (mod$pred$pred_ss)
    MSE_beta <- mean((pred_y - data$y)^2)

    return(c(ISE=ISE,MSE_b=MSE_b,MSE_beta=MSE_beta))
}


#      ISE     MSE_b  MSE_beta
#0.0649192 1.4870936 0.0181876



}
