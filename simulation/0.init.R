library(survival)
library(rpart)
library(lme4)
library(ipred)


survi<- function
#  Initialize
(y, offset, parms, wt) {
    if (!missing(offset) && length(offset) > 0)
        warning("offset argument ignored")

    sfun <- function(yval, dev, wt, ylevel, digits ) {
        paste("  mean=", format(signif(yval, digits)),
              ", MSE=" , format(signif(dev/wt, digits)),
              sep = '') }
    environment(sfun) <- .GlobalEnv
    #list(y = data.frame(y), parms = parms, numresp = 1, numy = ncol(y), summary = sfun)
    list(y = y, parms = parms, numresp = 1, numy = ncol(y), summary = sfun)
}


surve <- function
(y, wt, parms){

    if (is.null(parms$stable)){ parms$stable <- FALSE }
    y = data.frame(y)
    if (parms$LTRC){
        colnames(y)[1:4] <- c('start','end','event','biomarker')
        if (parms$test_stat=='rsq'){
            formulay1 <-  formulay2 <- Surv(start, end, event) ~ biomarker
        } else {
            formulay1 <- Surv(start, end, event) ~ . - biomarker
            formulay2 <- Surv(start, end, event) ~ .
        }
    } else {
        colnames(y)[1:3] <- c('time','event','biomarker')
        if (parms$test_stat=='rsq'){
            formulay1 <- formulay2 <- Surv(time, event) ~ biomarker
        } else {
            formulay1 <- Surv(time, event) ~ . - biomarker
            formulay2 <- Surv(time, event) ~ .
        }
    }

    node_val <- get_node_val(formulay1, formulay2, y, test_stat=parms$test_stat, stable=parms$stable)
    list(label = node_val, deviance = node_val)
    # deviance: it should be closely related to the split criteria.
    # label: does not matter, but we use node_val here.
}


survs_v1<- function
# Split
# version 1: if using log likelihood, the root level is treating potential two children 
# groups as same, thus same coefficients.
(y, wt, x, parms, continuous)
{
    y = data.frame(y)
    if (parms$LTRC){
        colnames(y)[1:4] <- c('start','end','event','biomarker')
        if (parms$test_stat=='rsq'){
            formulay1 <-  formulay2 <- Surv(start, end, event) ~ biomarker
        } else {
            formulay1 <- Surv(start, end, event) ~ . - biomarker
            formulay2 <- Surv(start, end, event) ~ .
        }
    } else {
        colnames(y)[1:3] <- c('time','event','biomarker')
        if (parms$test_stat=='rsq'){
            formulay1 <- formulay2 <- Surv(time, event) ~ biomarker
        } else {
            formulay1 <- Surv(time, event) ~ . - biomarker
            formulay2 <- Surv(time, event) ~ .
        }
    }

    root_val <- get_node_val(formulay1, formulay2, y, test_stat= parms$test_stat)

    if (continuous) {
        # continuous x variable: do all the logistic regressions
        n <- nrow(y)
        goodness <- double(n-1)
        direction <- goodness
        for (i in 1:(n-1)) {
            if (x[i] != x[i+1]) {
                result <- tryCatch({
                    left_val <- get_node_val(formulay1, formulay2, y[1:i,], test_stat= parms$test_stat)
                    right_val <- get_node_val(formulay1, formulay2, y[(i+1):n,], test_stat= parms$test_stat)
                    get_split_utility(root_val, left_val, right_val, parms$test_stat)
                }, error = function(e){ c(0, sign(1))})#, warning = function(w){c(0,sign(0))})
                goodness[i] = result[1]; direction[i] = result[2]

            }
        }
    } else {
        # Categorical X variable
        n <- nrow(y)
        ux <- sort(unique(x))
        nx <- length(ux)
        goodness <- double(nx-1)
        direction <- goodness

        for (i in 1:(nx-1)){
            #next_start <- min(which(x > ux[i]))
            next_start <- min(which(x == ux[i+1]))
            result <- tryCatch({
                left_val <- get_node_val(formulay1, formulay2, y[1:(next_start-1),], test_stat= parms$test_stat)
                right_val <- get_node_val(formulay1, formulay2, y[next_start:n,], test_stat= parms$test_stat)
                get_split_utility(root_val, left_val, right_val, parms$test_stat)
            }, error = function(e){ c(0, sign(0))})#, warning = function(w){c(0,sign(0))})
            goodness[i] = result[1]
        }
        names(goodness) <- ux[1:(nx-1)]
        direction <- ux
    }

    list(goodness=goodness, direction=direction)
}

survs_v2<- function
# Split
# version 2: if using log likelihood, the root level is using 
# coefficients for different groups (thus, changes after every split)
(y, wt, x, parms, continuous)
{
    y = data.frame(y)
    if (parms$LTRC){
        colnames(y)[1:4] <- c('start','end','event','biomarker')
        if (parms$test_stat=='rsq'){
            formulay1 <-  formulay2 <- Surv(start, end, event) ~ biomarker
        } else {
            formulay1 <- Surv(start, end, event) ~ . - biomarker
            formulay2 <- Surv(start, end, event) ~ .
        }
    } else {
        colnames(y)[1:3] <- c('time','event','biomarker')
        if (parms$test_stat=='rsq'){
            formulay1 <- formulay2 <- Surv(time, event) ~ biomarker
        } else {
            formulay1 <- Surv(time, event) ~ . - biomarker
            formulay2 <- Surv(time, event) ~ .
        }
    }


    if (parms$test_stat=='rsq'){
        root_val <- get_node_val(formulay1, formulay2, y, test_stat= parms$test_stat)
    } 


    if (continuous) {
        # continuous x variable: do all the logistic regressions
        n <- nrow(y)
        goodness <- double(n-1)
        direction <- goodness

        for (i in 1:(n-1)) {
            if (x[i] != x[i+1]) {
                result <- tryCatch({
                    left_val <- get_node_val(formulay1, formulay2, y[1:i,], test_stat= parms$test_stat)
                    right_val <- get_node_val(formulay1, formulay2, y[(i+1):n,], test_stat= parms$test_stat)

                    if (parms$test_stat=='lrt'){
                        if (parms$LTRC){
                            addcols <- y[,-c(1:4)] * c(rep(1,i), rep(0,n-i))
                        } else {
                            addcols <- y[,-c(1:3)] * c(rep(1,i), rep(0,n-i))
                        }
                        colnames(addcols) <- paste0(colnames(addcols),'_dummy')
                        y2 <- cbind(y, addcols)
                        root_val <- get_node_val(formulay1, formulay2, y2, test_stat= parms$test_stat)
                    }
                    get_split_utility(root_val, left_val, right_val, parms$test_stat)
                }, error = function(e){ c(0, sign(0))})#, warning = function(w){c(0,sign(0))})

                goodness[i] = result[1]; direction[i] = result[2]
            }
        }
    } else {
        # Categorical X variable
        n <- nrow(y)
        ux <- sort(unique(x))
        nx <- length(ux)
        goodness <- double(nx-1)
        direction <- goodness

        for (i in 1:(nx-1)){
            #next_start <- min(which(x > ux[i]))
            next_start <- min(which(x == ux[i+1]))
            result <- tryCatch({
                left_val <- get_node_val(formulay1, formulay2, y[1:(next_start-1),], test_stat= parms$test_stat)
                right_val <- get_node_val(formulay1, formulay2, y[next_start:n,], test_stat= parms$test_stat)
                if (parms$test_stat=='lrt'){
                    if (parms$LTRC){
                        addcols <- y[,-c(1:4)] * c(rep(1,next_start-1), rep(0,n-next_start+1))
                    } else {
                        addcols <- y[,-c(1:3)] * c(rep(1,next_start-1), rep(0,n-next_start+1))
                    }
                    colnames(addcols) <- paste0(colnames(addcols),'_dummy')
                    y2 <- cbind(y2, addcols)
                    root_val <- get_node_val(formulay1, formulay2, y2, test_stat= parms$test_stat)

                }
                get_split_utility(root_val, left_val, right_val, parms$test_stat)
            }, error = function(e){ c(0, sign(0))})#, warning = function(w){c(0,sign(0))})
            goodness[i] = result[1]
        }

        names(goodness) <- ux[1:(nx-1)]
        direction <- ux
 
    }
    list(goodness=goodness, direction=direction)
}

survs_v3<- function
# Split
# version 1: if using log likelihood, the root level is treating potential two children 
# groups as same, thus same coefficients.
(y, wt, x, parms, continuous)
{
    y = data.frame(y)
    if (parms$LTRC){
        colnames(y)[1:4] <- c('start','end','event','biomarker')
        if (parms$test_stat=='rsq'){
            formulay1 <-  formulay2 <- Surv(start, end, event) ~ biomarker
        } else {
            formulay1 <- Surv(start, end, event) ~ . - biomarker
            formulay2 <- Surv(start, end, event) ~ .
        }
    } else {
        colnames(y)[1:3] <- c('time','event','biomarker')
        if (parms$test_stat=='rsq'){
            formulay1 <- formulay2 <- Surv(time, event) ~ biomarker
        } else {
            formulay1 <- Surv(time, event) ~ . - biomarker
            formulay2 <- Surv(time, event) ~ .
        }
    }

    if (is.null(parms$min_nevent)){ parms$min_nevent <- 1 }
    if (is.null(parms$stop_thre)){ parms$stop_thre <- 1 } # expected Chi_1
    if (is.null(parms$stable)){ parms$stable <- FALSE }

    nevents <- sum(y[,'event'])
    root_val <- get_node_val(formulay1, formulay2, y, test_stat=parms$test_stat, stable=parms$stable)

    if (nevents <= parms$min_nevent*2 | root_val < parms$stop_thre){
        if (continuous){
            goodness <- rep(0,nrow(y)-1); direction<-goodness;
        } else{
            ux <- sort(unique(x))
            goodness <- rep(0,length(ux)-1); direction <- ux
        }
    } else {
        if (continuous) {
            # continuous x variable: do all the logistic regressions
            n <- nrow(y)
            goodness <- double(n-1)
            direction <- goodness
            for (i in 1:(n-1)) {
                if (x[i] != x[i+1]) {
                    nevents_l <- sum(y$event[1:i])
                    nevents_r <- sum(y$event[(i+1):n])
                    if (nevents_l <= parms$min_nevent | nevents_r <= parms$min_nevent ){
                        result <- c(0,0)
                    } else{
                        result <- tryCatch({
                            left_val <- get_node_val(formulay1, formulay2, y[1:i,], test_stat= parms$test_stat, stable=parms$stable)
                            right_val <- get_node_val(formulay1, formulay2, y[(i+1):n,], test_stat= parms$test_stat, stable=parms$stable)
                            get_split_utility(root_val, left_val, right_val, parms$test_stat)
                        }, error = function(e){ c(0, sign(1))})#, warning = function(w){c(0,sign(0))})
                    }
                    goodness[i] = result[1]; direction[i] = result[2]
                }
            }
            goodness <- pmax(0,goodness)
        } else {
            # Categorical X variable
            n <- nrow(y)
            ux <- sort(unique(x))
            nx <- length(ux)
            goodness <- double(nx-1)
            direction <- goodness

            for (i in 1:(nx-1)){
                #next_start <- min(which(x > ux[i]))
                next_start <- min(which(x == ux[i+1]))
                nevents_l <- sum(y$event[1:(next_start-1)])
                nevents_r <- sum(y$event[next_start:n])
                if (nevents_l <= parms$min_nevent | nevents_r <= parms$min_nevent ){
                    result <- c(0,0)
                } else{
                    result <- tryCatch({
                        left_val <- get_node_val(formulay1, formulay2, y[1:(next_start-1),], test_stat= parms$test_stat, stable=parms$stable)
                        right_val <- get_node_val(formulay1, formulay2, y[next_start:n,], test_stat= parms$test_stat, stable=parms$stable)
                        get_split_utility(root_val, left_val, right_val, parms$test_stat)
                    }, error = function(e){ c(0, sign(1))})#, warning = function(w){c(0,sign(0))})
                }
                goodness[i] = result[1]
            }
            names(goodness) <- ux[1:(nx-1)]
            goodness <- pmax(0,goodness)
            direction <- ux
        }
    }
    list(goodness=goodness, direction=direction)
}


get_node_val_OLD <- function
(f1,f2,data, 
 loglik = c(1,0)[1]){
    if (!loglik ){
        ret <- get_rsquare(f1, data)
    } else if (loglik){
        ret <- get_loglik_diff(f1, f2, data)
    }
    return (ret)
}

get_node_val<- function
(f1,f2,data, 
 test_stat=c('rsq','lrt','wald')[2],
 stable=FALSE){
    if (test_stat == 'rsq'){
        ret <- get_rsquare(f1, data)
    } else if (test_stat=='lrt'){
        ret <- get_loglik_diff(f1, f2, data, stable)
    } else if (test_stat=='wald'){
        ret <- get_wald(f2, data)
    }
    return (ret)
}


get_wald <- function(f, data){
    ret <- tryCatch({
        bo<-0
        while(bo!=10){
            coxml <- coxph(f,data)
            if (is.na(coxml$coefficients[1])){
                bo <- bo+1
                data <- data[sample(c(1:nrow(data)),replace = TRUE),]
            } else  break 
        }
        (coxml$coefficients[1])^2/ coxml$var[1,1]
    }, error = function(e){Inf})
}

get_rsquare <- function(f, data){
    coxml <- coxph(f, data)
    logtest <- -2 * (coxml$loglik[1] - coxml$loglik[2])
    rsq <-  1-exp(-logtest/coxml$nevent)
    return(rsq)
}


get_loglik_diff <- function(f1, f2, data, stable=FALSE){
    ret <- tryCatch({
        bo<-0
        while(bo!=10){
            coxml1 = coxph(f1, data)
            if (coxml1$loglik[2] < coxml1$loglik[1]){
                bo <- bo+1
                data <- data[sample(c(1:nrow(data)),replace = TRUE),]
            } else  break 
        }

        bo<-0
        while(bo!=10){
            coxml2 = coxph(f2, data)
            if (coxml2$loglik[2] < coxml2$loglik[1]){
                bo <- bo+1
                data <- data[sample(c(1:nrow(data)),replace = TRUE),]
            } else  break 
        }

        loglik_diff <- 2*(coxml2$loglik[2] - coxml1$loglik[2] )

        if(stable){
            if (max(c(diag(coxml1$var), diag(coxml2$var))) > 1e5){
                loglik_diff <- Inf
            }
        }

        max(0,loglik_diff)
    }, error = function(e){Inf}) # Inf, thus the parent node will not split to this childnode

    return(ret)
}

get_split_utility <- function
(var_root, var_left, var_right, 
 test_stat=c('rsq','lrt','wald')[1]){
    if (test_stat=='rsq'){
        c(var_root- (var_left+var_right)/2,
          sign(var_left- var_right))
    } else {
        c( 
          var_root - (var_left+var_right),
          #var_root - (var_left+var_right)/2,
          sign(var_left- var_right))
    }
}



independence_score <- function(gidx, newdata, LTRC){
    group_uq <- sort(unique(gidx))
    ind_score <- matrix(0, nrow=length(group_uq),ncol=2)
    newdata = data.frame(newdata)
    if (LTRC){
        colnames(newdata)[1:4] <- c('start','end','event','biomarker')
        formulay1 <- Surv(start, end, event) ~ . - biomarker
        formulay2 <- Surv(start, end, event) ~ .
    } else{
        colnames(newdata)[1:3] <- c('time','event','biomarker')
        formulay1 <- Surv(time, event) ~ . - biomarker
        formulay2 <- Surv(time, event) ~ .
    }

    for (i in seq_along(group_uq)){
        g <- group_uq[i]
        data_sub <- newdata[gidx==g,]
        coxml1 <- coxph(formulay1, data_sub)
        coxml2 <- coxph(formulay2, data_sub)
        ind_score[i,1] <- nrow(data_sub)
        ind_score[i,2] <- coxml2$loglik[2] - coxml1$loglik[2] 
    }

    return(ind_score)
}


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



get_latent_class <- function(X1,X2,struct,member,seed=0){
    if (struct == 'tree'){
        W <- matrix(c( -1,-1, -1,1,1,-1,1,1),ncol=2,byrow=TRUE)
        if (member == 'partition'){
            g <- apply(cbind(X1,X2),1,function(x){
                           weights <- W %*% (x-0.5)
                           which.max(weights)})
        } else if (member == 'multinomial'){
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

        if (member == 'multinomial'){
            g_multi <- rep(0, length(g))
            for (i in seq_along(X1)){
                gi <- g[i]
                weights <- rep(0.1,4); weights[gi] <- 0.7
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

    g <- get_latent_class(X1,X2,FLAGS$struct, FLAGS$member, seed=seed)
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
        tmpX <- cbind(1,newdata[,c('X1','X2','X3','X4','X5')])
        if (ncol(coefs_multilogit)==7){
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
    pred_slopes <- matrix(coefs[coefstart:(coefstart+nclasses*3-1)],nrow=nclasses,ncol=3)

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



