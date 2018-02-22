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

    node_val <- get_node_val(formulay1, formulay2, y, test_stat= parms$test_stat)
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

    nevents <- sum(y[,'event'])
    root_val <- get_node_val(formulay1, formulay2, y, test_stat= parms$test_stat)

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
                            left_val <- get_node_val(formulay1, formulay2, y[1:i,], test_stat= parms$test_stat)
                            right_val <- get_node_val(formulay1, formulay2, y[(i+1):n,], test_stat= parms$test_stat)
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
                        left_val <- get_node_val(formulay1, formulay2, y[1:(next_start-1),], test_stat= parms$test_stat)
                        right_val <- get_node_val(formulay1, formulay2, y[next_start:n,], test_stat= parms$test_stat)
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
 test_stat = c('rsq','lrt','wald')[2]){
    if (test_stat == 'rsq'){
        ret <- get_rsquare(f1, data)
    } else if (test_stat=='lrt'){
        ret <- get_loglik_diff(f1, f2, data)
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


get_loglik_diff <- function(f1, f2, data){
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

#        if (loglik_diff < threshold){
#            RET <- rep(0,5)
#            for (i in seq_len(5)){
#                reorder <- sample(seq_len(nrow(data)),replace = FALSE)
#                RET[i]<-get_loglik_diff(f1, f2, data[reorder,])
#            }
#            loglik_diff <- median(RET)
#        }
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


# This function is wrong.
#get_IBS <- function(surv, pred){
#    times <- surv[,1]
#    N <- length(times)
#    scores <- apply(array(times),1,function(x){sbrier(surv, pred, x)})
#    torder <- order(times)
#    scores <- scores[torder]; times <- times[torder]
#
#    IBS <- sum(0.5*(scores[-1]+scores[-N]) * (times[-1]-times[-N] )) / (max(times)-min(times))
#    return (IBS)
#}
#

eval_tree_pred_OLD<- function # for y not longituindal
(data,dist, slopes, parms, node_idx){
    Nsub <- nrow(data)
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
    SURVS <- matrix(0,nrow=ntimes,ncol=nrow(data))

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
            org_idx <- seq_along(subset_id)[subset_id][x]
            SURVS[,org_idx] <- Shat

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

    surv <- Surv(data$time_Y, data$delta)
    IBS_pred <- sbrier(surv, SURVS)[1,1]
    ISE <- ISE/Nsub

    pred_parms[is.na(pred_parms)] <- 0
    MSE_b <- mean(rowSums((true_parms - pred_parms)^2))

    # ---------- Biomarker Prediction
    pred_y<- rep(0,Nsub)
    pred_ef <- matrix(0,ncol=2, nrow=Nsub)
    true_ef <- cbind(data$X1+data$X2, data$y - data$X1-data$X2)

    d_idx <- node_idx
    for (i in unique(d_idx)){
        subset_id <- d_idx== i
        tmpdata <- data[subset_id,]
        alp_d <- mean(tmpdata$y)
        beta_i  <- tmpdata$y - alp_d
        pred_ef[subset_id,] <- cbind(rep(alp_d,sum(subset_id)), beta_i)
    }
    MSE_beta <- mean(rowSums((true_ef- pred_ef)^2))

    return(c(ISE=ISE,IBS_pred=IBS_pred,MSE_b=MSE_b,MSE_beta=MSE_beta))
}

eval_lcmm_pred_OLD <- function
(data, dist, slopes, parms, mod){
    Nsub <- nrow(data)
    g <- 2*(data$X1)+data$X2+1
    true_parms <- slopes[g,]
    ebx <- exp(rowSums(slopes[g,] * cbind(data$X3,data$X4,data$X5)))

    #  coeff
    coefs <- mod$best
    coefstart <- 27
    pred_slopes <- matrix(coefs[coefstart:(coefstart+11)],nrow=4,ncol=3)

    predclass <- mod$pprob$class
    pred_parms <- pred_slopes[predclass,]
    MSE_b <- mean(rowSums((true_parms - pred_parms)^2))

    # Shat
    times <- mod$predSurv[,1]
    ntimes <- length(times)
    ISE <- 0

    for (x in c(1:nrow(data))){
        tmpclass <- mod$pprob$class[x]
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

    ISE <- ISE/Nsub

    # ---------- Biomarker Prediction
    pred_y <- mod$pred$pred_ss
    MSE_beta <- mean((pred_y - data$y)^2)

    return(c(ISE=ISE,MSE_b=MSE_b,MSE_beta=MSE_beta))
}


eval_tree_pred<- function
(data,dist, slopes, parms, node_idx){
    Nsub_pseudo  <- nrow(data)
    g <- 2*(data$X1)+data$X2+1
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

            tmpg <- 2*(tmpdata[x,'X1'])+tmpdata[x,'X2']+1
            tmpebx <- exp(sum(slopes[tmpg,] * tmpdata[x,c('X3','X4','X5')]))
            if (dist=='exponential'){
                Strue <- exp(-tmpebx*timessub*parms$lambda)
            } else if (dist=='weibulld' | dist=='weibulli'){
                Strue <- exp(-(timessub/parms$beta)^(parms$alp) * tmpebx)
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
(data, dist, slopes, parms, mod){
    Nsub_pseudo<- nrow(data)
    g <- 2*(data$X1)+data$X2+1
    true_parms <- slopes[g,]
    ebx <- exp(rowSums(slopes[g,] * cbind(data$X3,data$X4,data$X5)))


    #  coeff
    nclasses <- ncol(mod$pprob)-2
    coefs <- mod$best
    coefstart <- max(which(grepl('Weibull',names(coefs))))+1
    #coefstart <- 27
    pred_slopes <- matrix(coefs[coefstart:(coefstart+nclasses*3-1)],nrow=nclasses,ncol=3)

    predclass <- (mod$pprob$class)[data$ID]
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
