library(survival)
library(rpart)
library(lme4)
library(ipred)
library(prodlim)


eval_tree_pred<- function
(data,dist, slopes, parms, node_idx, g=NULL){
    Nsub_pseudo  <- nrow(data)
    if (is.null(g)){
        g <- 2*(data$X1)+data$X2+1
    } 
    ebx <- exp(rowSums(slopes[g,] * cbind(data$X3,data$X4,data$X5)))

    # ---------- Survival Prediction 
    # get true survival prob
    pred_parms <- matrix(0,nrow=Nsub_pseudo,ncol=3)
    true_parms <- slopes[g,]

    # get predicted survival prob
    d_idx <- node_idx
    formula <-Surv(time_L,time_Y,delta) ~ X3+X4+X5

    times <- sort(data$time_Y)
    timessub <- seq(from=min(times),to=max(times), length.out = 100)
    ntimes <- length(timessub)
    ISE <- 0

    # purity
    tmptable <- table(node_idx, g)
    purity <- sum(apply(tmptable,1,function(tb){ if(sum(tb>0)==1){sum(tb)} else {0}}))/nrow(data)

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

        KM <- FALSE
        if (class(mod)=="try-error") {
            mod <- coxph(Surv(time_L,time_Y,delta) ~ 1, tmpdata)
            KM <- TRUE
        }

        for (x in c(1:nrow(tmpdata))){
            if(!KM){
                Shat <- getsurv(survfit(mod,newdata=tmpdata[x,]),timessub)
            } else{
                Shat <- getsurv(survfit(mod),timessub)
            }

            #tmpg <- 2*(tmpdata[x,'X1'])+tmpdata[x,'X2']+1
            tmpg <- g[subset_id][x]
            tmpebx <- exp(sum(slopes[tmpg,] * tmpdata[x,c('X3','X4','X5')]))
            if (dist=='exponential'){
                Strue <- exp(-tmpebx*timessub*parms$lambda)
            } else if (dist=='weibulld' | dist=='weibulli'){
                Strue <- exp(-(timessub/parms$beta)^(parms$alp) * tmpebx)
            } else if (dist=='lognormal'){
                tmpparms <- parms[tmpg,]
                Strue <- (1-pnorm((log(timessub)-tmpparms[1])/tmpparms[2]))^tmpebx
            }


            scores <- (Shat - Strue)^2
            SE <- sum(0.5*(scores[-1]+scores[-ntimes]) * diff(timessub)) / diff(range(timessub))
            ISE <- ISE + SE
        }

        if (!KM){
            pred_parms[subset_id,] <- rep(mod$coefficients,each=sum(subset_id))
        } else {
            pred_parms[subset_id,] <- rep(c(0,0,0),each=sum(subset_id))
        }
    }

    ISE <- ISE/Nsub_pseudo

    pred_parms[is.na(pred_parms)] <- 0
    MSE_b <- mean(rowSums((true_parms - pred_parms)^2))

    # ---------- Biomarker Prediction
    if(length(unique(data$ID))== nrow(data)){
        ymod <- lm(y ~ d_idx, data=data)
    } else {
        ymod <- lmer(y ~ d_idx + (1|ID),data=data)
    }
    pred_y<- predict(ymod)
    MSE_y <- mean((pred_y - data$y)^2)

    return(round(c(ISE=ISE,MSE_b=MSE_b,MSE_y=MSE_y,purity=purity),4))
}

eval_lcmm_pred<- function
(data, dist, slopes, parms, mod,g=NULL,inter=TRUE){
    Nsub_pseudo<- nrow(data)
    if(is.null(g)){
        g <- 2*(data$X1)+data$X2+1
    }
    true_parms <- slopes[g,]
    ebx <- exp(rowSums(slopes[g,] * cbind(data$X3,data$X4,data$X5)))

    #  coeff
    nclasses <- ncol(mod$pprob)-2
    coefs <- mod$best
    coefstart <- max(which(grepl('Weibull',names(coefs))))+1
    pred_slopes <- matrix(coefs[coefstart:(coefstart+nclasses*3-1)],nrow=nclasses,ncol=3)

    # predclass
    predclass <- (mod$pprob$class)[data$ID]
    #coefend <- min(which(grepl('Weibull',names(coefs))))-1
    #coefs_multilogit <- matrix(coefs[1:coefend],nrow=nclasses-1)
    #tmpX <- cbind(1,data[,c('X1','X2','X3','X4','X5')])
    #if(inter){
    #    tmpX <- cbind(tmpX,data[,'X1']*data[,'X2'])
    #}

    #tmp1 <- as.matrix(tmpX) %*% t(coefs_multilogit)
    #predclass <- apply(tmp1, 1,function(x){ which.max(exp(c(x,0)))})


    pred_parms <- pred_slopes[predclass,]
    MSE_b <- mean(rowSums((true_parms - pred_parms)^2))

    # purity
    tmptable <- table(predclass, g)
    purity <- sum(apply(tmptable,1,function(tb){ if(sum(tb>0)==1){sum(tb)} else {0}}))/nrow(data)

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
        } else if (dist=='lognormal'){
            tmpparms <- parms[tmpg,]
            Strue <- (1-pnorm((log(times)-tmpparms[1])/tmpparms[2]))^tmpebx
        }

        scores <- (Shat - Strue)^2
        SE <- sum(0.5*(scores[-1]+scores[-ntimes]) * diff(times)) / diff(range(times))
        ISE <- ISE + SE
    }

    ISE <- ISE/Nsub_pseudo

    # ---------- Biomarker Prediction
    pred_y <- (mod$pred$pred_ss)
    MSE_y <- mean((pred_y - data$y)^2)

    return(round(c(ISE=ISE,MSE_b=MSE_b,MSE_y=MSE_y, purity=purity),4))
}


get_parms <- function(dist){
    parms <- switch(dist,
                    exponential=list(lambda=0.1),
                    weibulld=list(alp=0.9, beta=1),
                    weibulli=list(alp=3, beta=2),
                    lognormal=matrix(c(2,0.3,1.7,0.2,1.3,0.3,0.5,0.5),byrow=TRUE,ncol=2))

    slopes <- switch(dist,
                     exponential=matrix(c(0,0,0,0.56,0.56,0.09,0.92,0.92,0.15,1.46,1.46,0.24),ncol=3,byrow=TRUE),
                     weibulld=matrix(c(-1.17,-1.17,-0.19,-0.66,-0.66,-0.11,-0.55,-0.55,-0.09,0,0,0),ncol=3,byrow=TRUE),
                     weibulli=matrix(c(-3.22,-3.22,-0.54,-2.26,-2.26,-0.38, -1.53,-1.53,-0.26, 0,0,0),ncol=3,byrow=TRUE),
                     lognormal=matrix(rep(0,12),ncol=3))

    lam_D <- switch(dist,
                    exponential=list('0.2'=c(0.026,0.048,0.086,0.156),
                                     '0.5'=c(0.114,0.232,0.396,0.74)),
                    weibulld=list('0.2'=c(0.024,0.072,0.084,0.238),
                                  '0.5'=c(0.13,0.318,0.384,0.956)),
                    weibulli=list('0.2'=c(0.024,0.04,0.076,0.19),
                                  '0.5'=c(0.084,0.168,0.224,0.64)),
                    lognormal= list('0.2'=c(0.03,0.05,0.08, 0.15),
                                    '0.5'=c(0.1,0.12,0.2,0.4)))


    return (list(parms=parms, slopes=slopes, lam_D=lam_D))
}



get_latent_class <- function(X1,X2,struct,member,seed=0,majprob=NULL){

    if(is.null(majprob)){
        Wmul <- 1
    } else {
        Wlist <- switch(struct,
                        tree=list('0.25' = 0, '0.5'=1, '0.7'=2, '0.85'=4.5, '0.99'=50), 
                        linear=list('0.25' = 0, '0.5'=1.5, '0.7'=3.8, '0.85'=9, '0.99'=150))
        Wmul <- Wlist[[as.character(majprob)]]
    }

    if (struct == 'tree'){
        W <- matrix(c( -1,-1, -1,1,1,-1,1,1),ncol=2,byrow=TRUE)
        if (member == 'partition'){
            g <- apply(cbind(X1,X2),1,function(x){
                           weights <- W %*% (x-0.5)
                           which.max(weights)})
        } else if (member == 'multinomial'){
            W <- W*Wmul
            g <- apply(cbind(X1,X2),1,function(x){
                           weights <- exp(2*W %*% (x-0.5))
                           weights <- weights/sum(weights)
                           cumweights <- cumsum(weights)
                           which(runif(1)  < cumweights)[1] })
        }
    } else if (struct == 'XOR'){# {(0,0), (1,1)} -> g = 4; {(1,0),(0,1)} -> g=1
        if (member == 'partition'){
            g <- (2*X1-1) * (2*X2-1)
            g <- 1.5*g + 2.5
        } else if (member == 'multinomial'){
            W <- c(1,-1)
            X <- (2*X1-1)*(2*X2-1)
            g <- apply(array(X),1,function(x){
                           weights <- exp(x * W)
                           weights <- weights/sum(weights)
                           cumweights <- cumsum(weights)
                           which(runif(1)  < cumweights)[1] })
            g <- ifelse(g==2,1,4)
        }

    } else if (struct == 'linear'){
        W <- matrix(c(0.8,-0.6,0.9,0.5,-0.8,0.6,0.5,0.9),byrow=TRUE,ncol=2) # random sample from sphere
        if (member == 'partition'){
            g <- apply(cbind(X1,X2),1,function(x){
                           weights <- W %*% (x-2)
                           which.max(weights)})
        } else if (member == 'multinomial'){
            W <- W*Wmul
            g <- apply(cbind(X1,X2),1,function(x){
                           weights <- exp(W %*% (x-2))
                           weights <- weights/sum(weights)
                           cumweights <- cumsum(weights)
                           which(runif(1)  < cumweights)[1] })
        }
    } else {
        if (struct=='nonlinear'){
            f1 <-  X1^2+X2^2 < 0.75^2
            f2 <- (X1-0)^2+(X2-1)^2 < 0.75^2
            g <- ifelse(f1, ifelse(f2, 4,1),ifelse(f2, 2,3))
        } else if (struct == 'asym'){
            g <- ifelse(X1>0.75, 1, ifelse(X2<0.33, 2, ifelse(X2<0.67,3,4)))
        } else if (struct == 'quad'){
            g <- ifelse(X1<0.25,1,ifelse(X1<0.5,4,ifelse(X1<0.75,2,3)))
        }

        if(is.null(majprob)){majprob <- 0.7}
        if (member == 'multinomial'){
            g_multi <- rep(0, length(g))
            for (i in seq_along(X1)){
                gi <- g[i]
                weights <- rep((1-majprob)/3,4); weights[gi] <- majprob
                #weights <- rep(0.1,4); weights[gi] <- 0.7
                cumweights <- cumsum(weights)
                g_multi[i] <- which(runif(1)  < cumweights)[1]
            }
            g <- g_multi
        }
    }
    return (g)
}

gen_data <- function(FLAGS, PARMS, seed){

    Nsub <- FLAGS$Nsub
    censor_rate <- FLAGS$censor
    dist <- FLAGS$dist
    parms <- PARMS$parms; slopes <- PARMS$slopes; lam_D <- PARMS$lam_D
    sd_ranef <- 0.2
    sd_e <- 0.1

    set.seed(seed)
    # X1 - X5
    if (FLAGS$continuous){
        if(FLAGS$struct == 'linear'){
            X1 <- round(runif(2*Nsub, min=1,max=3),2)
            X2 <- round(runif(2*Nsub, min=1,max=3),2)
        } else {
            X1 <- round(runif(2*Nsub),2)
            X2 <- round(runif(2*Nsub),2)
        }
    } else {
        if(FLAGS$struct == 'linear'){
            X1 <- sample(c(1:3),2*Nsub,replace=TRUE)
            X2 <- sample(c(1:3),2*Nsub,replace=TRUE)
        } else if (FLAGS$struct == 'nonlinear'){
            stop("Nonlinear must have continuous X1 and X2.")
        } else {
            X1 <- as.numeric(runif(2*Nsub)>0.5)
            X2 <- as.numeric(runif(2*Nsub)>0.5)
        }
    }
    X3 <- as.numeric(runif(2*Nsub)>0.5)
    X4 <- round(runif(2*Nsub),1)
    X5 <- sample(c(1:5),2*Nsub,replace=TRUE)
    X <- cbind(X1,X2,X3,X4,X5)

    g <- get_latent_class(X1,X2,FLAGS$struct, FLAGS$member, seed=seed, FLAGS$majprob)
    if(FLAGS$dist == 'lognormal'){ 
        ebx <- rep(1 , 2*Nsub)
        tmp_parms <- parms[g,]
    } else{
        ebx <- exp(rowSums(slopes[g,] * cbind(X3,X4,X5))) 
        tmp_parms <- parms
    }

    time_T <- gen_model3_survival(ebx, dist, tmp_parms)
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
                           if (FLAGS$alg == 'jlctree'){
                               tmp_time_L <- tmp_time[1:num_i]
                               tmp_time_Y <- tmp_time[-1]
                               tmp_delta <- c(rep(0,num_i-1),delta[i])
                           } else if (FLAGS$alg == 'jlcmm'){
                               tmp_time_L <- rep(time_L[i],num_i)
                               tmp_time_Y <- rep(time_Y[i],num_i)
                               tmp_delta <- rep(delta[i],num_i)
                           }

                           ret <- cbind(ID=i, X1=rep(X[i,1],num_i),
                                        X2=rep(X[i,2],num_i), X3=rep(X[i,3],num_i),
                                        X4=rep(X[i,4],num_i), X5=rep(X[i,5],num_i),
                                        time_L=tmp_time_L, time_Y=tmp_time_Y, delta=tmp_delta)})


    ranef <- rnorm(Nsub, sd=sd_ranef)
    ranefs <- ranef[LTRC_data$ID]
    fixef <- c(0,1,1,2)
    pseudo_g <- g[LTRC_data$ID]
    y <- fixef[pseudo_g] + ranefs + rnorm(nrow(LTRC_data), sd=sd_e) 
    data <- cbind(LTRC_data,y)

    if(!is.null(FLAGS$extra)){
        if( FLAGS$extra ){ # add extra unrelated predictors
        X6 <- (round(runif(Nsub),2))[LTRC_data$ID]
        X7 <- (round(runif(Nsub),2))[LTRC_data$ID]
        X8 <- (as.numeric(runif(Nsub)>0.5))[LTRC_data$ID]
        X9 <- (round(runif(Nsub),1))[LTRC_data$ID]
        X10 <- (sample(c(1:5),Nsub,replace=TRUE))[LTRC_data$ID]
        data <- cbind(data, X6,X7,X8,X9,X10)
    }}
    return(list(data=data,pseudo_g=pseudo_g))
}

predict_class <- function(obj, newdata){

    if(class(obj) == 'rpart'){
        obj2<- obj
        obj2$frame[grepl('leaf',obj2$frame$var),]$yval <- sort(unique(obj2$where))
        test_class <- predict(obj2,newdata)
    } else if (class(obj) == 'Jointlcmm'){
        nclasses <- ncol(obj$pprob)-2
        coefs <- obj$best
        coefend <- min(which(grepl('Weibull',names(coefs))))-1
        coefs_multilogit <- matrix(coefs[1:coefend],nrow=nclasses-1)
        tmpX <- cbind(1,newdata[,paste0('X',c(1:5))])
        if ('X6' %in% colnames(newdata)){ tmpX <- cbind(tmpX,newdata[,paste0('X',c(6:10))]) }
        if (ncol(coefs_multilogit) %in% c(7,12)){
            tmpX <- cbind(tmpX,newdata[,'X1']*newdata[,'X2'])
        } 


        linearval <- as.matrix(tmpX) %*% t(coefs_multilogit)
        test_class <- t(apply(linearval, 1,function(x){ exps=exp(c(x,0)); exps/sum(exps)}))
    }

    return (test_class)
}

get_tree_ISE <- function(mod, subdata, subg, evaltimes, dist, slopes, parms,  KM){

    ISE <- 0
    for (x in c(1:nrow(subdata))){
        if(!KM){
            Shat <- getsurv(survfit(mod,newdata=subdata[x,]),evaltimes)
        } else{
            Shat <- getsurv(survfit(mod),evaltimes)
        }

        tmpg <- subg[x]
        tmpebx <- exp(sum(slopes[tmpg,] * subdata[x,c('X3','X4','X5')]))

        if (dist=='exponential'){
            Strue <- exp(-tmpebx*evaltimes*parms$lambda)
        } else if (dist=='weibulld' | dist=='weibulli'){
            Strue <- exp(-(evaltimes/parms$beta)^(parms$alp) * tmpebx)
        } else if (dist=='lognormal'){
            tmpparms <- parms[tmpg,]
            Strue <- (1-pnorm((log(evaltimes)-tmpparms[1])/tmpparms[2]))^tmpebx
        }

        scores <- (Shat - Strue)^2
        ntimes <- length(evaltimes)
        SE <- sum(0.5*(scores[-1]+scores[-ntimes]) * diff(evaltimes)) / diff(range(evaltimes))
        ISE <- ISE+SE
    }

    return (ISE)
}

eval_tree_pred_inout <- function
(data, data_test, dist, slopes, parms, 
 idx, idx_test, 
 g, g_test){

    Nobs <- nrow(data); Nobs_test <- nrow(data_test)
    ebx <- exp(rowSums(slopes[g,] * data[,c('X3','X4','X5')]))
    ebx_test <- exp(rowSums(slopes[g_test,] * data_test[,c('X3','X4','X5')]))

    uniqd <- unique(idx)

    # ---------- Survival Prediction 
    # get true survival prob
    pred_parms <- matrix(0,nrow=Nobs,ncol=3)
    true_parms <- slopes[g,]

    # get predicted survival prob
    formula <-Surv(time_L,time_Y,delta) ~ X3+X4+X5

    ntimes <- 100
    evaltimes <- seq(from=min(data$time_Y),to=max(data$time_Y),length.out=ntimes)
    evaltimes_test <- seq(from=min(data_test$time_Y),to=max(data_test$time_Y),length.out=ntimes)
    ISE <- 0; ISE_test <- 0

    for (i in uniqd){
        # subsets
        sid <- idx == i; sdata <- data[sid,]; sg <- g[sid]
        sid_test <- idx_test==i; sdata_test <- data_test[sid_test,]; sg_test <- g_test[sid_test]

        KM <- FALSE; err <-0
        while(err!=10){
            mod  <- try(coxph(formula, sdata),silent=TRUE)
            if (class(mod)=="try-error") {
                err <- err +1;  sdata <- sdata[sample(c(1:nrow(sdata)),replace = TRUE),]
            } else break 
        }
        if (class(mod)=="try-error") { mod <- coxph(Surv(time_L,time_Y,delta) ~ 1, sdata); KM <- TRUE }
        orgsdata <- data[sid,]; #if sdata has been reordered, then survfit(mod) will have trouble.

        ISE <- ISE + get_tree_ISE(mod, orgsdata, sg, evaltimes, dist, slopes, parms, KM)
        if(nrow(sdata_test)>0){
            ISE_test <- ISE_test + get_tree_ISE(mod, sdata_test, sg_test, evaltimes_test, dist, slopes,parms, KM)
        }

        if (!KM){
            pred_parms[sid,] <- rep(mod$coefficients,each=nrow(sdata))
        } else {
            pred_parms[sid,] <- rep(c(0,0,0),each=nrow(sdata))
        }
    }
    ISE <- ISE/Nobs; ISE_test <- ISE_test/Nobs_test

    pred_parms[is.na(pred_parms)] <- 0
    MSE_b <- mean(rowSums((true_parms - pred_parms)^2))

    # ---------- Biomarker Prediction
    data$idx <- factor(idx)
    data_test$idx <- factor(idx_test); data_test$ID <- 0

    if(length(unique(idx))==1){
        #        ymod <- lm(y ~ idx, data=data)
        ymod <- lmer(y ~ X1+X2+X3+X4+X5+(1|ID), data=data)
    } else {
        #        ymod <- lmer(y ~ idx + (1|ID),data=data)
        ymod <- lmer(y ~ X1+X2+X3+X4+X5+(1|idx) + (1|ID),data=data)
    }
    predy<- predict(ymod); predy_test <- predict(ymod,newdata=data_test, allow.new.levels=TRUE)
    MSE_y <- mean((predy - data$y)^2); MSE_y_test <- mean((predy_test - data_test$y)^2)

    # purity
    tmptable <- table(idx, g)
    purity <- sum(apply(tmptable,1,function(tb){ if(sum(tb>0)==1){sum(tb)} else {0}}))/Nobs

    tmptable <- table(idx_test, g_test)
    purity_test <- sum(apply(tmptable,1,function(tb){ if(sum(tb>0)==1){sum(tb)} else {0}}))/Nobs_test


    return(round(c(ISE=ISE,MSE_b=MSE_b,MSE_y=MSE_y,purity=purity,
                   ISE_test=ISE_test, MSE_y_test = MSE_y_test,purity_test=purity_test),4))
}

get_lcmm_ISE <- function(mod, data, g, predclass, dist, slopes, parms){
    times <- mod$predSurv[,1]
    ntimes <- length(times)
    ISE <- 0

    nclasses <- ncol(mod$pprob)-2
    coefs <- mod$best
    coefstart <- max(which(grepl('Weibull',names(coefs))))+1
    pred_slopes <- matrix(coefs[coefstart:(coefstart+nclasses*3-1)],nrow=nclasses,ncol=3)

    if (length(dim(predclass))==0){ avg <- FALSE} else {avg <- TRUE}

    for (x in c(1:nrow(data))){
        if (!avg){
            tmpc <- predclass[x]
            tmpebx <- exp(sum(pred_slopes[tmpc,] * data[x,c('X3','X4','X5')]))
            Shat <- exp(-tmpebx*mod$predSurv[,paste0('event1.CumRiskFct',tmpc)])
        } else {
            Shat_raw <- matrix(0,ncol=nclasses,nrow=ntimes)
            for (tmpc in c(1:nclasses)){
                tmpebx <- exp(sum(pred_slopes[tmpc,] * data[x,c('X3','X4','X5')]))
                Shat_raw[,tmpc] <- exp(-tmpebx*mod$predSurv[,paste0('event1.CumRiskFct',tmpc)])
            }
            Shat<- c(Shat_raw %*% predclass[x,])
        }

        tmpg <- g[x]
        tmpebx <- exp(sum(slopes[tmpg,] * data[x,c('X3','X4','X5')]))
        if (dist=='exponential'){
            Strue <- exp(-tmpebx*times*parms$lambda)
        } else if (dist=='weibulld' | dist=='weibulli'){
            Strue <- exp(-(times/parms$beta)^(parms$alp) * tmpebx)
        } else if (dist=='lognormal'){
            tmpparms <- parms[tmpg,]
            Strue <- (1-pnorm((log(times)-tmpparms[1])/tmpparms[2]))^tmpebx
        }

        scores <- (Shat - Strue)^2
        SE <- sum(0.5*(scores[-1]+scores[-ntimes]) * diff(times)) / diff(range(times))
        ISE <- ISE + SE
    }
    ISE <- ISE / nrow(data)
    return (ISE)
}

eval_lcmm_pred_inout <- function
(data, data_test, dist, slopes, parms, mod,
 g, g_test){

    Nobs <- nrow(data)
    true_parms <- slopes[g,]

    #  coeff
    nclasses <- ncol(mod$pprob)-2
    coefs <- mod$best
    coefstart <- max(which(grepl('Weibull',names(coefs))))+1
    pred_slopes <- matrix(coefs[coefstart:(coefstart+nclasses*3-1)],nrow=nclasses)

    # insample predclass
    predclass_in <- (mod$pprob$class)[data$ID]
    pred_parms <- pred_slopes[predclass_in,]
    MSE_b <- mean(rowSums((true_parms - pred_parms)^2))

    # out of sample predclass: a vector of probabilities
    predclass_test <- predict_class(mod, data_test)
    predclass_test_max <- apply(predclass_test,1,which.max)


    # ---------- Purity
    tmptable <- table(predclass_in, g)
    purity <- sum(apply(tmptable,1,function(tb){ if(sum(tb>0)==1){sum(tb)} else {0}}))/Nobs

    tmptable <- table(predclass_test_max, g_test)
    purity_test <- sum(apply(tmptable,1,function(tb){ if(sum(tb>0)==1){sum(tb)} else {0}}))/nrow(data_test)

    # ---------- ISE
    ISE <- get_lcmm_ISE(mod, data, pseudo_g, predclass_in, dist, slopes, parms)
    ISE_test_max <- get_lcmm_ISE(mod, data_test, pseudo_g_test, predclass_test_max, dist, slopes, parms)
    ISE_test_avg <- get_lcmm_ISE(mod, data_test, pseudo_g_test, predclass_test, dist, slopes, parms)

    # ---------- Biomarker Prediction
    predy <- (mod$pred$pred_ss)
    predy_test_raw <- predictY(mod,newdata=data_test)$pred
    predy_test_max <- apply(array(seq_along(predclass_test_max)),1,function(x){predy_test_raw[x,predclass_test_max[x]]})
    predy_test_avg <- rowSums(predy_test_raw * predclass_test)

    MSE_y <- mean((predy - data$y)^2)
    MSE_y_test_max <- mean((predy_test_max - data_test$y)^2)
    MSE_y_test_avg <- mean((predy_test_avg - data_test$y)^2)

    return(round(c(ISE=ISE,MSE_b=MSE_b,MSE_y=MSE_y, purity=purity,
                   ISE_test_max=ISE_test_max,MSE_y_test_max=MSE_y_test_max,
                   ISE_test_avg=ISE_test_avg,MSE_y_test_avg=MSE_y_test_avg, purity_test=purity_test),4))
}


gen_data_timevar <- function(FLAGS, PARMS, seed, survvar=FALSE){

    Nsub <- FLAGS$Nsub
    censor_rate <- FLAGS$censor
    dist <- FLAGS$dist
    parms <- PARMS$parms; slopes <- PARMS$slopes; lam_D <- PARMS$lam_D
    sd_ranef <- 0.2
    sd_e <- 0.1

    set.seed(seed)
    # X1 - X5
    if (FLAGS$continuous){
        if(FLAGS$struct == 'linear'){
            X1 <- round(runif(2*Nsub, min=1,max=3),2)
            X2 <- round(runif(2*Nsub, min=1,max=3),2)
            X1_next <- round(pmax(pmin(X1 + runif(2*Nsub, min=-0.6, max=0.6),3),1),2)
            X2_next <- round(pmax(pmin(X2 + runif(2*Nsub, min=-0.6, max=0.6),3),1),2)
        } else {
            X1 <- round(runif(2*Nsub),2)
            X2 <- round(runif(2*Nsub),2)
            X1_next <- round(pmax(pmin(X1 + runif(2*Nsub, min=-0.3, max=0.3),1),0),2)
            X2_next <- round(pmax(pmin(X2 + runif(2*Nsub, min=-0.3, max=0.3),1),0),2)
        }
    } 
    X3 <- as.numeric(runif(2*Nsub)>0.5)
    X4 <- round(runif(2*Nsub),1)
    X5 <- sample(c(1:5),2*Nsub,replace=TRUE)

    if (survvar){
        X3_next <- as.numeric(runif(2*Nsub)>0.5)
        X4_next <- round(pmax(pmin(X4 + runif(2*Nsub, min=-0.3, max=0.3),1),0),1)
        X5_next <- pmax(pmin(X5+sample(c(0,-1,1),2*Nsub,replace=TRUE),5),1)
    } else {
        X3_next <- X3; X4_next <- X4; X5_next <- X5
    }

    if (FLAGS$majprob==1){
        g <- get_latent_class(X1,X2,FLAGS$struct, 'partition', seed=seed, FLAGS$majprob)
        g_next <- get_latent_class(X1_next,X2_next,FLAGS$struct, 'partition', seed=seed, FLAGS$majprob)
    } else {
        g <- get_latent_class(X1,X2,FLAGS$struct, 'multinomial', seed=seed, FLAGS$majprob)
        g_next <- get_latent_class(X1_next,X2_next,FLAGS$struct, 'multinomial', seed=seed, FLAGS$majprob)
    }

    X <- cbind(X1,X2,X3,X4,X5, X1_next, X2_next, X3_next, X4_next, X5_next)

    ebx1 <- exp(rowSums(slopes[g,] * cbind(X3,X4,X5))) 
    ebx2 <- exp(rowSums(slopes[g_next,] * cbind(X3_next,X4_next,X5_next))) 
    ebx <- cbind(ebx1,ebx2)
    changepoint <- runif(2*Nsub, min=1,max=3)
    time_T <- gen_model4_survival(ebx,dist,parms,changepoint)
    time_L <- runif(2*Nsub, min=0, max=1)

    time_tokeep <- time_L < time_T
    time_L <- time_L[time_tokeep][1:Nsub]
    time_T <- time_T[time_tokeep][1:Nsub]
    X <- X[time_tokeep,][1:Nsub,]
    g <- g[time_tokeep][1:Nsub]
    g_next <- g_next[time_tokeep][1:Nsub]
    ebx <- ebx[time_tokeep,][1:Nsub,]
    changepoint <- changepoint[time_tokeep][1:Nsub]

    if (censor_rate==0){ 
        time_C <- Inf
    } else{
        time_C <- time_L + rexp(Nsub,lam_D[[censor_rate]][g])
    }

    delta <- as.numeric(time_C >= time_T)
    time_Y <- pmin(time_T, time_C)

    LTRC_data <- ldply(array(seq_len(Nsub)),function(i){
                           if(changepoint[i] < time_Y[i]){
                               tmp_time <- c(time_L[i],changepoint[i],time_Y[i])
                               num_i <- 2
                           } else{
                               tmp_time <- c(time_L[i],time_Y[i])
                               num_i <- 1
                           }
                           if (FLAGS$alg == 'jlctree'){
                               tmp_time_L <- tmp_time[1:num_i]
                               tmp_time_Y <- tmp_time[-1]
                               tmp_delta <- c(rep(0,num_i-1),delta[i])
                               ret <- cbind(ID=i, 
                                            X1=c(X[i,1], X[i,6])[1:num_i], X2=c(X[i,2],X[i,7])[1:num_i],
                                            X3=c(X[i,3], X[i,8])[1:num_i], X4=c(X[i,4],X[i,9])[1:num_i],
                                            X5=c(X[i,5], X[i,10])[1:num_i])

                           } else if (FLAGS$alg == 'jlcmm'){
                               tmp_time_L <- rep(time_L[i],num_i)
                               tmp_time_Y <- rep(time_Y[i],num_i)
                               tmp_delta <- rep(delta[i],num_i)
                               tmp_changepoint <- rep(changepoint[i], num_i)
                               ret <- cbind(ID=i, changepoint=tmp_changepoint,
                                            X1=rep(X[i,1], num_i), X2=rep(X[i,2],num_i),
                                            X3=rep(X[i,3], num_i), X4=rep(X[i,4],num_i),X5=rep(X[i,5],num_i),
                                            X3_true=c(X[i,3], X[i,8])[1:num_i], X4_true=c(X[i,4],X[i,9])[1:num_i],
                                            X5_true=c(X[i,5], X[i,10])[1:num_i])

                           }

                           ret <- cbind(ret, g = c(g[i],g_next[i])[1:num_i],
                                        time_L=tmp_time_L, time_Y=tmp_time_Y, delta=tmp_delta)})

    ranef <- rnorm(Nsub, sd=sd_ranef)
    ranefs <- ranef[LTRC_data$ID]
    fixef <- c(0,1,1,2)
    pseudo_g <- LTRC_data$g
    y <- fixef[pseudo_g] + ranefs + rnorm(nrow(LTRC_data), sd=sd_e) 
    data <- cbind(LTRC_data,y)

    return(list(data=data,pseudo_g=pseudo_g))
}


eval_tree_pred_inout_timevar <- function
(data, data_test, dist, slopes, parms, 
 idx, idx_test, 
 g, g_test,
 firstobs=FALSE){

    Nobs <- nrow(data); Nobs_test <- nrow(data_test)
    #ebx <- exp(rowSums(slopes[g,] * data[,c('X3','X4','X5')]))
    #ebx_test <- exp(rowSums(slopes[g_test,] * data_test[,c('X3','X4','X5')]))

    uniqd <- unique(idx)

    if(firstobs){ # only use first observation in fitting.
        data_tmp <- data.table(data)
        data_tmp <- data_tmp[,list(X1=X1[1],X2=X2[1],X3=X3[1],X4=X4[1],X5=X5[1],g=g[1],time_L,time_Y,delta,y),by =ID]
        data_touse <- data_tmp

        data_tmp <- data.table(data_test)
        data_tmp <- data_tmp[,list(X1=X1[1],X2=X2[1],X3=X3[1],X4=X4[1],X5=X5[1],g=g[1],time_L,time_Y,delta,y),by =ID]
        data_test_touse <- data_tmp
    } else {
        data_touse <- data
        data_test_touse <- data_test
    }

    # ---------- Survival Prediction 
    # get true survival prob
    pred_parms <- matrix(0,nrow=Nobs,ncol=3)
    true_parms <- slopes[g,]

    # get predicted survival prob
    formula <-Surv(time_L,time_Y,delta) ~ X3+X4+X5

    for (i in uniqd){
        # subsets
        sid <- idx == i; sdata <- data_touse[sid,]; sg <- g[sid]
        sid_test <- idx_test==i; sdata_test <- data_test_touse[sid_test,]; sg_test <- g_test[sid_test]

        KM <- FALSE; err <-0
        while(err!=10){
            mod  <- try(coxph(formula, sdata,model=TRUE),silent=TRUE)
            if (class(mod)=="try-error") {
                err <- err +1;  sdata <- sdata[sample(c(1:nrow(sdata)),replace = TRUE),]
            } else break 
        }
        if (class(mod)=="try-error") { mod <- coxph(Surv(time_L,time_Y,delta) ~ 1, sdata); KM <- TRUE }
        assign(paste0('mod',i),mod)
        assign(paste0('KM',i),KM)

        if (!KM){
            pred_parms[sid,] <- rep(mod$coefficients,each=nrow(sdata))
        } else {
            pred_parms[sid,] <- rep(c(0,0,0),each=nrow(sdata))
        }
    }

    # ---------- Parameter Prediction
    pred_parms[is.na(pred_parms)] <- 0
    MSE_b <- mean(rowSums((true_parms - pred_parms)^2))


    # ---------- ISE in sample
    ntimes <- 100
    evaltimes <- seq(from=min(data$time_Y),to=max(data$time_Y),length.out=ntimes)
    ISE <- 0

    for (x in unique(data$ID)){
        sdata <- subset(data,ID==x)
        sg <- subset(g, data$ID==x)
        si <- subset(idx, data$ID==x)

        changepoint <- sdata$time_Y[1]
        changepoint_timeidx <-  which(changepoint < evaltimes)[1]

        if(nrow(sdata)==1){
            mod <- get(paste0('mod',si));
            KM <- get(paste0('KM',si));
            ISE <- ISE + get_tree_ISE(mod, sdata, sg, evaltimes, dist, slopes,parms,KM)
        } else {
            # Shat
            SHAT <- matrix(0,nrow=2,ncol=length(evaltimes))
            prob <- rep(0,2)
            for (i in c(1:2)){
                mod <- get(paste0('mod',si[i])); KM <- get(paste0('KM',si[i]));
                if(!KM){ 
                    SHAT[i,] <- getsurv(survfit(mod,newdata=sdata[i,]),evaltimes)
                    prob[i] <- getsurv(survfit(mod,newdata=sdata[i,]),changepoint)
                } else{ 
                    SHAT[i,]<- getsurv(survfit(mod),evaltimes) 
                    prob[i] <- getsurv(survfit(mod),changepoint)
                }
            }
            Shat <- c(SHAT[1,1:(changepoint_timeidx-1)], SHAT[2,(changepoint_timeidx:ntimes)]*prob[1]/prob[2])

            # True S
            STRUE <- matrix(0,nrow=2,ncol=length(evaltimes))
            prob <- rep(0,2)
            for(i in c(1:2)){
                tmpebx <- exp(sum(slopes[sg[i],] * sdata[i,c('X3','X4','X5')]))
                if (dist=='exponential'){
                    STRUE[i,] <- exp(-tmpebx*evaltimes*parms$lambda)
                    prob[i] <- exp(-tmpebx*changepoint*parms$lambda)
                } else if (dist=='weibulld' | dist=='weibulli'){
                    STRUE[i,] <- exp(-(evaltimes/parms$beta)^(parms$alp) * tmpebx)
                    prob[i] <- exp(-(changepoint/parms$beta)^(parms$alp) * tmpebx)
                } 
            }

            Strue<- c(STRUE[1,1:(changepoint_timeidx-1)], 
                      STRUE[2,(changepoint_timeidx:ntimes)]*prob[1]/prob[2])


            scores <- (Shat - Strue)^2
            ISE <- ISE + sum(0.5*(scores[-1]+scores[-ntimes]) * diff(evaltimes)) / diff(range(evaltimes))
        }
    }

    ISE <- ISE / length(unique(data$ID))

    # ---------- ISE out of sample
    evaltimes_test <- seq(from=min(data_test$time_Y),to=max(data_test$time_Y),length.out=ntimes)
    ISE_test <- 0

    for (x in unique(data_test$ID)){
        sdata <- subset(data_test,ID==x)
        sg <- subset(g_test, data_test$ID==x)
        si <- subset(idx_test, data_test$ID==x)

        changepoint <- sdata$time_Y[1]
        changepoint_timeidx <-  which(changepoint < evaltimes_test)[1]

        if(nrow(sdata)==1){
            mod <- get(paste0('mod',si));
            KM <- get(paste0('KM',si));
            ISE_test <- ISE_test + get_tree_ISE(mod, sdata, sg, evaltimes_test, dist, slopes,parms,KM)
        } else {
            SHAT <- matrix(0,nrow=2,ncol=ntimes)
            prob <- rep(0,2)
            for (i in c(1:2)){
                mod <- get(paste0('mod',si[i])); KM <- get(paste0('KM',si[i]));
                if(!KM){ 
                    SHAT[i,] <- getsurv(survfit(mod,newdata=sdata[i,]),evaltimes_test)
                    prob[i] <- getsurv(survfit(mod,newdata=sdata[i,]),changepoint)
                } else{ 
                    SHAT[i,]<- getsurv(survfit(mod),evaltimes_test) 
                    prob[i] <- getsurv(survfit(mod),changepoint)
                }
            }
            Shat <- c(SHAT[1,1:(changepoint_timeidx-1)], SHAT[2,(changepoint_timeidx:ntimes)]*prob[1]/prob[2])

            # True S
            STRUE <- matrix(0,nrow=2,ncol=ntimes)
            prob <- rep(0,2)
            for(i in c(1:2)){
                tmpebx <- exp(sum(slopes[sg[i],] * sdata[i,c('X3','X4','X5')]))
                if (dist=='exponential'){
                    STRUE[i,] <- exp(-tmpebx*evaltimes_test*parms$lambda)
                    prob[i] <- exp(-tmpebx*changepoint*parms$lambda)
                } else if (dist=='weibulld' | dist=='weibulli'){
                    STRUE[i,] <- exp(-(evaltimes_test/parms$beta)^(parms$alp) * tmpebx)
                    prob[i] <- exp(-(changepoint/parms$beta)^(parms$alp) * tmpebx)
                } 
            }

            Strue<- c(STRUE[1,1:(changepoint_timeidx-1)], 
                      STRUE[2,(changepoint_timeidx:ntimes)]*prob[1]/prob[2])

            scores <- (Shat - Strue)^2
            ISE_test <- ISE_test + sum(0.5*(scores[-1]+scores[-ntimes]) * diff(evaltimes_test)) / diff(range(evaltimes_test))
        }
    }

    ISE_test <- ISE_test / length(unique(data_test$ID))

    # ---------- Biomarker Prediction
    data$idx <- factor(idx)
    data_test$idx <- factor(idx_test); data_test$ID <- 0

    if(length(unique(idx))==1){
        #        ymod <- lm(y ~ idx, data=data)
        ymod <- lmer(y ~ X1+X2+X3+X4+X5+(1|ID), data=data)
    } else {
        #        ymod <- lmer(y ~ idx + (1|ID),data=data)
        ymod <- lmer(y ~ X1+X2+X3+X4+X5+(1|idx) + (1|ID),data=data)
    }
    predy<- predict(ymod); predy_test <- predict(ymod,newdata=data_test, allow.new.levels=TRUE)
    MSE_y <- mean((predy - data$y)^2); MSE_y_test <- mean((predy_test - data_test$y)^2)

    # purity
    tmptable <- table(idx, g)
    purity <- sum(apply(tmptable,1,function(tb){ if(sum(tb>0)==1){sum(tb)} else {0}}))/Nobs

    tmptable <- table(idx_test, g_test)
    purity_test <- sum(apply(tmptable,1,function(tb){ if(sum(tb>0)==1){sum(tb)} else {0}}))/Nobs_test

    return(round(c(ISE=ISE,MSE_b=MSE_b,MSE_y=MSE_y,purity=purity,
                   ISE_test=ISE_test, MSE_y_test = MSE_y_test,purity_test=purity_test),4))
}



get_lcmm_ISE_timevar <- function(mod, data, g, predclass, dist, slopes, parms){
    times <- mod$predSurv[,1]
    ntimes <- length(times)
    ISE <- 0

    nclasses <- ncol(mod$pprob)-2
    coefs <- mod$best
    coefstart <- max(which(grepl('Weibull',names(coefs))))+1
    pred_slopes <- matrix(coefs[coefstart:(coefstart+nclasses*3-1)],nrow=nclasses,ncol=3)

    if (length(dim(predclass))==0){ avg <- FALSE} else {avg <- TRUE}

    for (x in unique(data$ID)){
        sdata <- subset(data,ID==x)
        sg <- subset(g, data$ID==x)
        si <- subset(predclass, data$ID==x)
        changepoint <- sdata$changepoint[1]
        changepoint_timeidx <-  which(changepoint < times)[1]
        ns <- nrow(sdata)

        if (!avg){
            tmpc <- si[1]
            tmpebx <- exp(sum(pred_slopes[tmpc,] * sdata[1,c('X3','X4','X5')]))
            Shat <- exp(-tmpebx*mod$predSurv[,paste0('event1.CumRiskFct',tmpc)])
        } else {
            Shat_raw <- matrix(0,ncol=nclasses,nrow=ntimes)
            for (tmpc in c(1:nclasses)){
                tmpebx <- exp(sum(pred_slopes[tmpc,] * sdata[1,c('X3','X4','X5')]))
                Shat_raw[,tmpc] <- exp(-tmpebx*mod$predSurv[,paste0('event1.CumRiskFct',tmpc)])
            }
            Shat<- c(Shat_raw %*% si[1,])
        }


        STRUE <- matrix(0,nrow=2,ncol=ntimes)
        prob <- rep(0,2)
        for(i in c(1:ns)){
            tmpebx <- exp(sum(slopes[sg[i],] * sdata[i,c('X3_true','X4_true','X5_true')]))
            if (dist=='exponential'){
                STRUE[i,] <- exp(-tmpebx*times*parms$lambda)
                prob[i] <- exp(-tmpebx*changepoint*parms$lambda)
            } else if (dist=='weibulld' | dist=='weibulli'){
                STRUE[i,] <- exp(-(times/parms$beta)^(parms$alp) * tmpebx)
                prob[i] <- exp(-(changepoint/parms$beta)^(parms$alp) * tmpebx)
            } 
        }
        if (ns==1){
            Strue <- STRUE[1,]
        } else {
            Strue<- c(STRUE[1,1:(changepoint_timeidx-1)], 
                      STRUE[2,(changepoint_timeidx:ntimes)]*prob[1]/prob[2])
        }
        scores <- (Shat - Strue)^2
        SE <- sum(0.5*(scores[-1]+scores[-ntimes]) * diff(times)) / diff(range(times))
        ISE <- ISE + SE
    }
    ISE <- ISE / length(unique(data$ID))
    return (ISE)
}

eval_lcmm_pred_inout_timevar <- function
(data, data_test, dist, slopes, parms, mod,
 g, g_test){

    Nobs <- nrow(data)
    true_parms <- slopes[g,]

    #  coeff
    nclasses <- ncol(mod$pprob)-2
    coefs <- mod$best
    coefstart <- max(which(grepl('Weibull',names(coefs))))+1
    pred_slopes <- matrix(coefs[coefstart:(coefstart+nclasses*3-1)],nrow=nclasses)

    # insample predclass
    predclass_in <- (mod$pprob$class)[data$ID]
    pred_parms <- pred_slopes[predclass_in,]
    MSE_b <- mean(rowSums((true_parms - pred_parms)^2))

    # out of sample predclass: a vector of probabilities
    predclass_test <- predict_class(mod, data_test)
    predclass_test_max <- apply(predclass_test,1,which.max)

    # ---------- Purity
    tmptable <- table(predclass_in, g)
    purity <- sum(apply(tmptable,1,function(tb){ if(sum(tb>0)==1){sum(tb)} else {0}}))/Nobs

    tmptable <- table(predclass_test_max, g_test)
    purity_test <- sum(apply(tmptable,1,function(tb){ if(sum(tb>0)==1){sum(tb)} else {0}}))/nrow(data_test)

    # ---------- ISE
    ISE <- get_lcmm_ISE_timevar(mod, data, pseudo_g, predclass_in, dist, slopes, parms)
    ISE_test_max <- get_lcmm_ISE_timevar(mod, data_test, pseudo_g_test, predclass_test_max, dist, slopes, parms)
    ISE_test_avg <- get_lcmm_ISE_timevar(mod, data_test, pseudo_g_test, predclass_test, dist, slopes, parms)

    # ---------- Biomarker Prediction
    predy <- (mod$pred$pred_ss)
    predy_test_raw <- predictY(mod,newdata=data_test)$pred
    predy_test_max <- apply(array(seq_along(predclass_test_max)),1,function(x){predy_test_raw[x,predclass_test_max[x]]})
    predy_test_avg <- rowSums(predy_test_raw * predclass_test)

    MSE_y <- mean((predy - data$y)^2)
    MSE_y_test_max <- mean((predy_test_max - data_test$y)^2)
    MSE_y_test_avg <- mean((predy_test_avg - data_test$y)^2)

    return(round(c(ISE=ISE,MSE_b=MSE_b,MSE_y=MSE_y, purity=purity,
                   ISE_test_max=ISE_test_max,MSE_y_test_max=MSE_y_test_max,
                   ISE_test_avg=ISE_test_avg,MSE_y_test_avg=MSE_y_test_avg, purity_test=purity_test),4))
}






get_pred_curv <- function(obj, newdata, evaltimes){
    neval <- length(evaltimes)
    ret <- matrix(0,nrow=neval, ncol=nrow(newdata))
    newdata_id <- unique(newdata$id)

    if(class(obj) == 'survreg'){
        pct <- 1:999/1000
        for (j in newdata_id){
            tmpidx <- newdata$id == j
            tmpdata <-newdata[tmpidx,]
            ntmp <- sum(tmpidx)
            tmpdata$tstop[ntmp] <- max(evaltimes)

            ptime <- predict(obj, newdata=tmpdata, type='quantile',p=pct)
            tmpsurv <- matrix(0,nrow=neval,ncol=ntmp)
            if(ntmp==1){
                fakeobj <- list(surv = 1-pct, time=ptime);class(fakeobj) <- 'survfit'
                tmpsurv[,1] <- getsurv(fakeobj, evaltimes)
                ret[,tmpidx] <- tmpsurv
            } else {
                for (z in seq_len(ntmp)){
                    fakeobj <- list(surv = 1-pct, time=ptime[z,]);class(fakeobj) <- 'survfit'
                    tmpsurv[,z] <- getsurv(fakeobj, evaltimes)
                }
                final_curv <- agg_curv(tmpsurv, evaltimes, tmpdata$tstop, evaltimes)
                ret[,tmpidx] <- rep(final_curv, ntmp)
            }
        }
    } else if (class(obj) == 'coxph'){
        # obj is fit with strata
        # use individual=TRUE to automatically aggregate across pseudo-subjects.
        for (j in newdata_id){
            tmpidx <- newdata$id == j
            tmpdata <- newdata[tmpidx,]
            ntmp <- sum(tmpidx)
            tmpdata$tstop[ntmp] <- max(evaltimes)

            if(ntmp==1){ survfitmod <- survfit(obj, newdata = tmpdata,id='class')
            } else { 
                survfitmod <- survfit(obj, newdata = tmpdata,individual=TRUE) 
                survfitmod$time <- survfitmod$time+tmpdata[1,tstart]
            }
            final_curv <- getsurv(survfitmod, evaltimes)
            ret[, tmpidx] <- rep( final_curv, ntmp)
        }
    } else if (class(obj)=='coxph.null'){
        final_curv <- getsurv(survfit(obj),evaltimes)
        ret <- rep(final_curv, nrow(newdata))
    }
    return (ret)
}


agg_curv <- function(surv, time, tstops, evaltimes){
    ntmp <- length(tstops)

    if (ntmp>1){
        breaks <- c(-Inf,tstops[-ntmp],Inf)
        cut_ret <- as.numeric(cut(time, breaks))
        surv_ret <- c()
        prevprob <- 1
        currtimeidx <- 0
        for(z in seq_len(ntmp)){
            if(any(cut_ret==z)){
                surv_ret <- c(surv_ret, surv[cut_ret==z,z]*prevprob)
                currtimeidx <- max(which(cut_ret==z))
                if (z < ntmp & currtimeidx < length(time)){
                    if (surv[currtimeidx+1, z+1] < 1e-4){
                        prevprob <- 0
                    } else {prevprob <- prevprob*surv[currtimeidx+1,z]/surv[currtimeidx+1,z+1] }
                }
            }
        }
        fake_obj <- list(surv=surv_ret, time=time); class(fake_obj) <- 'survfit'
        return(getsurv(fake_obj,evaltimes))
    } else {
        return(surv)
    }
}


eval_tree_pred_inout_jlctree <- function
(data, data_test, dist, slopes, parms, 
 idx, idx_test, 
 g, g_test, jlctree, timevar){

    Nobs <- nrow(data); Nobs_test <- nrow(data_test)
#    ebx <- exp(rowSums(slopes[g,] * data[,c('X3','X4','X5')]))
#    ebx_test <- exp(rowSums(slopes[g_test,] * data_test[,c('X3','X4','X5')]))

    # ---------- Coxph Parameter 
    # get true survival prob
    pred_parms <- matrix(0,nrow=Nobs,ncol=3)
    true_parms <- slopes[g,]

    # get pred survival prob
    pred_parms <- array(0, c(Nobs, ncol(true_parms), 3))
    if(!is.null(jlctree$coxph_model_diffh_diffs)){
        pred_parms_raw <- get_cox_coef(jlctree$coxph_model_diffh_diffs)
        pred_parms[,,1] <- pred_parms_raw[as.character(idx),]
    }
    if(!is.null(jlctree$coxph_model_diffh)){
        pred_parms_raw <- get_cox_coef(jlctree$coxph_model_diffh)
        pred_parms[,,2] <- rep(pred_parms_raw, each = Nobs)
    }
    if(!is.null(jlctree$coxph_model_diff)){
        pred_parms_raw <- get_cox_coef(jlctree$coxph_model_diffs)
        pred_parms[,,3] <- pred_parms_raw[as.character(idx),]
    }
    MSE_b <- apply(pred_parms,3,function(x){mean(rowSums((true_parms-x)^2))})

    # ---------- Survival Prediction 
    ISE <- c(get_tree_ISE_jlctree(data,idx,pseudo_g, jlctree$coxph_model_diffh, PARMS$slopes, PARMS$parms, FLAGS$dist, timevar),
             get_tree_ISE_jlctree(data,idx,pseudo_g, jlctree$coxph_model_diffh_diffs, PARMS$slopes, PARMS$parms, FLAGS$dist, timevar),
             get_tree_ISE_jlctree(data,idx,pseudo_g, jlctree$coxph_model_diffs, PARMS$slopes, PARMS$parms, FLAGS$dist, timevar))

    ISE_test <- c(get_tree_ISE_jlctree(data_test,idx_test,pseudo_g_test,jlctree$coxph_model_diffh, PARMS$slopes, PARMS$parms, FLAGS$dist, timevar),
                  get_tree_ISE_jlctree(data_test,idx_test,pseudo_g_test,jlctree$coxph_model_diffh_diffs, PARMS$slopes, PARMS$parms, FLAGS$dist, timevar),
                  get_tree_ISE_jlctree(data_test,idx_test,pseudo_g_test,jlctree$coxph_model_diffs, PARMS$slopes, PARMS$parms, FLAGS$dist, timevar))


    # ---------- Biomarker Prediction
    predy <- predict(jlctree$lmm_model)
    data_test$node <- as.factor(idx_test); data_test$ID <- 0
    predy_test <- predict(jlctree$lmm_model,newdata=data_test,allow.new.levels=TRUE)
    MSE_y <- mean((predy - data$y)^2); MSE_y_test <- mean((predy_test - data_test$y)^2)

    # purity
    tmptable <- table(idx, g)
    purity <- sum(apply(tmptable,1,function(tb){ if(sum(tb>0)==1){sum(tb)} else {0}}))/Nobs

    tmptable <- table(idx_test, g_test)
    purity_test <- sum(apply(tmptable,1,function(tb){ if(sum(tb>0)==1){sum(tb)} else {0}}))/Nobs_test


    return(round(c(ISE=ISE,MSE_b=MSE_b,MSE_y=MSE_y,purity=purity,
                   ISE_test=ISE_test, MSE_y_test = MSE_y_test,purity_test=purity_test),4))
}


get_cox_coef <- function(coxph_model){
    if(is.null(coxph_model$xlevels$node)){
        pred_parms <- matrix(coef(coxph_model),nrow=1)
        # when there is one class, for consistency name it '1'
        rownames(pred_parms) <- '1'
    } else {
        all_coef <- coef(coxph_model)
        nodes <- paste0('node',coxph_model$xlevels$node)
        all_coef[is.na(all_coef)] <- 0

        baselevel <- nodes[1]
        base_parms <- all_coef[!grepl('node',names(all_coef))]
        pred_parms <- base_parms
        for (d in c(2:length(nodes))){
            tmp_parms <- all_coef[grepl(nodes[d],names(all_coef))]
            tmp_parms <- base_parms + tmp_parms[-1]+tmp_parms[1]
            pred_parms <- rbind(pred_parms, tmp_parms)
        }
        rownames(pred_parms) <- coxph_model$xlevels$node
    }
    return (pred_parms)
}



get_tree_ISE_jlctree <- function
(data, idx, g, coxph_model, slopes, parms, dist, timevar, weighted=TRUE)
{
    if(is.null(coxph_model)){
        ISE <- 0
    } else {
        data$node <- as.factor(idx)
        SE <- rep(0,length(unique(data$ID)))

        ntimes <- 100
        evaltimes <- seq(from=min(data$time_Y),to=max(data$time_Y),length.out=ntimes)

        for(j in unique(data$ID)){
            tmpdata <- subset(data, ID==j);ntmp <- nrow(tmpdata)
            tmpg <- subset(g, data$ID==j)
            tmpdata[ntmp,'time_Y'] <- Inf
            survmod <- survfit(coxph_model, newdata=tmpdata,individual=TRUE)
            survmod$time <- survmod$time + tmpdata[1,'time_L']

            pred_surv <- getsurv(survmod, evaltimes)
            true_surv <- get_true_surv(slopes[tmpg,,drop=FALSE], parms, dist, tmpdata, evaltimes, timevar)

            scores <- (pred_surv - true_surv)^2
            ntimes <- length(evaltimes)
            SE[j] <- sum(0.5*(scores[-1]+scores[-ntimes]) * diff(evaltimes)) / diff(range(evaltimes))
        }

        if(weighted){
            weights <- table(data$ID)
            ISE <- sum(weights*SE)/sum(weights)
        } else {
            ISE <- mean(SE)
        }
    }
    return(ISE)
}

get_true_surv <- function(slopes, parms, dist, data, evaltimes,timevar)
{
    ntimes <- length(evaltimes)
    if(timevar & nrow(data)>1){
        changepoint <- data$time_Y[1]
        changepoint_timeidx <-  which(changepoint < evaltimes)[1]

        # True S
        STRUE <- matrix(0,nrow=2,ncol=length(evaltimes))
        prob <- rep(0,2)
        for(i in c(1:2)){
            tmpebx <- exp(sum(slopes[i,] * data[i,c('X3','X4','X5')]))
            if (dist=='exponential'){
                STRUE[i,] <- exp(-tmpebx*evaltimes*parms$lambda)
                prob[i] <- exp(-tmpebx*changepoint*parms$lambda)
            } else if (dist=='weibulld' | dist=='weibulli'){
                STRUE[i,] <- exp(-(evaltimes/parms$beta)^(parms$alp) * tmpebx)
                prob[i] <- exp(-(changepoint/parms$beta)^(parms$alp) * tmpebx)
            } 
        }
        Strue<- c(STRUE[1,1:(changepoint_timeidx-1)], 
                  STRUE[2,(changepoint_timeidx:ntimes)]*prob[1]/prob[2])
    } else {
        tmpebx  <- exp(sum(slopes[1,] * data[1,c('X3','X4','X5')]))

        if (dist=='exponential'){
            Strue <- exp(-tmpebx*evaltimes*parms$lambda)
        } else if (dist=='weibulld' | dist=='weibulli'){
            Strue <- exp(-(evaltimes/parms$beta)^(parms$alp) * tmpebx)
        } else if (dist=='lognormal'){
            Strue <- (1-pnorm((log(evaltimes)-parms[1])/parms[2]))^tmpebx
        }
    }

   return(Strue)
}


