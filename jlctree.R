library(survival)
library(rpart)


survi <- function
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
#  Evaludate
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


survs <- function
# Split 
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
    if (is.null(parms$split_add)){ parms$split_add <- 20} # add parms$split_add to goodness in splitting utility

    nevents <- sum(y[,'event'])
    root_val <- get_node_val(formulay1, formulay2, y, test_stat=parms$test_stat, stable=parms$stable)

    if (nevents <= parms$min_nevent*2 | root_val < parms$stop_thre){
        if (continuous){
            goodness <- rep(-Inf,nrow(y)-1); direction<-goodness;
        } else{
            ux <- sort(unique(x))
            goodness <- rep(-Inf,length(ux)-1); direction <- ux
        }
    } else {
        if (continuous) {
            # continuous x variable
            n <- nrow(y)
            goodness <- rep(-Inf,n-1)
            direction <- double(n-1)
            for (i in 1:(n-1)) {
                if (x[i] != x[i+1]) {
                    nevents_l <- sum(y$event[1:i])
                    nevents_r <- sum(y$event[(i+1):n])
                    if (nevents_l <= parms$min_nevent | nevents_r <= parms$min_nevent ){
                        result <- c(-Inf,0)
                    } else{
                        result <- tryCatch({
                            left_val <- get_node_val(formulay1, formulay2, y[1:i,], test_stat= parms$test_stat, stable=parms$stable)
                            right_val <- get_node_val(formulay1, formulay2, y[(i+1):n,], test_stat= parms$test_stat, stable=parms$stable)
                            get_split_utility(root_val, left_val, right_val, parms$test_stat)
                        }, error = function(e){ c(-Inf, sign(1))})#, warning = function(w){c(0,sign(0))})
                    }
                    goodness[i] = result[1]; direction[i] = result[2]
                }
            }
            #goodness <- pmax(0,goodness)
            goodness <- goodness+parms$split_add
        } else {
            # Categorical X variable
            n <- nrow(y)
            ux <- sort(unique(x))
            nx <- length(ux)

            goodness <- rep(-Inf, nx-1)
            direction <- double(nx-1)

            for (i in 1:(nx-1)){
                next_start <- min(which(x == ux[i+1]))
                nevents_l <- sum(y$event[1:(next_start-1)])
                nevents_r <- sum(y$event[next_start:n])
                if (nevents_l <= parms$min_nevent | nevents_r <= parms$min_nevent ){
                    result <- c(-Inf,0)
                } else{
                    result <- tryCatch({
                        left_val <- get_node_val(formulay1, formulay2, y[1:(next_start-1),], test_stat= parms$test_stat, stable=parms$stable)
                        right_val <- get_node_val(formulay1, formulay2, y[next_start:n,], test_stat= parms$test_stat, stable=parms$stable)
                        get_split_utility(root_val, left_val, right_val, parms$test_stat)
                    }, error = function(e){ c(-Inf, sign(1))})#, warning = function(w){c(0,sign(0))})
                }
                goodness[i] = result[1]
            }
            names(goodness) <- ux[1:(nx-1)]
            goodness <- goodness+parms$split_add
            #goodness <- pmax(0,goodness)
            direction <- ux
        }
    }
    list(goodness=goodness, direction=direction)
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
        # sometimes coxph fails, but works again after shuffling the data.
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
    }, error = function(e){Inf}) # Inf, thus the parent node will not split into this child node

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


jlctree <- function
(survival, ## Surv(start,end,event)~ cox ph vars
 classmb, ## Tree split vars
 fixed, ## biomarker ~ fixed effects vars
 random, ## ~ random effects vars
 subject, ## variable name that indicates a subject
 data,
 maxng=6, ## prune trees to have at most maxng nodes
 greedy=TRUE,## if TRUE, set cp=-Inf in rpart.control.
 parms=list(min_nevent=5)
){

    require('lme4')
    survlist <- list(eval=surve, split=survs, init=survi)
    if(greedy){control <- rpart.control(cp=-Inf)} else {control <- rpart.control()}
    if(is.null(parms$LTRC)){ parms$LTRC=1}
    if(is.null(parms$test_stat)){parms$test_stat='lrt'}
    if(is.null(parms$stop_thre)){parms$stop_thre=3.84}
    if(is.null(parms$stable)){parms$stable=TRUE}
    if(is.null(parms$min_nevent)){parms$min_nevent=5}

    surv_vars <- labels(terms(survival)); if(length(surv_vars)==0){surv_vars='1'}
    surv_char <- deparse(survival[[2]]) #setdiff(all.vars(survival),surv_vars)
    y_var <- deparse(fixed[[2]]) #all.vars(fixed)[1]

    rpart_formula <- paste0('cbind(',paste0(c(surv_char,y_var,surv_vars),collapse=','),')',
                            paste0(as.character(classmb),collapse=''))
    rpart_y <- model.frame(formula(rpart_formula),data=data)

    # ---------- Fit jlctree, and prune to have at most maxng nodes.
    tree <- rpart(rpart_y, control=control, method=survlist, parms=parms)

    if(!is.null(maxng)){
        if (max(tree$cptable[,'nsplit']) > maxng-1){
            minsplit <- max(tree$cptable[,'nsplit'][tree$cptable[,'nsplit']<=maxng-1])
            prunecp <- min( tree$cptable[,'CP'][tree$cptable[,'nsplit']<=minsplit])
            tree <- prune(tree, prunecp)
        }
    }

    # ---------- Fit lmm model.
    data$node <- as.factor(tree$where)
    rand_vars <- labels(terms(random)); if (length(rand_vars)==0){rand_vars='1'}
    if(length(unique(tree$where))==1){
        lmer_formula <- as.formula(paste0(Reduce(paste0, deparse(fixed)), '+', 
                                          #'(', paste0(rand_vars,collapse='+') ,'|node)+',
                                          paste0('(1|',subject,')')))
    } else {
        lmer_formula <- as.formula(paste0(Reduce(paste0, deparse(fixed)), '+', 
                                          '(', paste0(rand_vars,collapse='+') ,'|node)+',
                                          paste0('(1|',subject,')')))
    }
    lmm_model <- lmer(lmer_formula,data=data)

    # ---------- Fit three versions of cox model:
    # diffh: different baseline hazards, diffs: different cox slopes.

    if(length(unique(tree$where))==1){
        coxph_formula_diffh_diffs <- as.formula(paste0(deparse(survival[[2]]),'~',
                                                       paste0(surv_vars,collapse='+')))
        coxph_formula_diffh <- coxph_formula_diffh_diffs 
        coxph_formula_diffs <- coxph_formula_diffh_diffs 
    } else {
        coxph_formula_diffh_diffs <- as.formula(paste0(deparse(survival[[2]]),'~',
                                                       paste0(paste0(surv_vars,'*node'),collapse='+'),'+', 'strata(node)'))
        coxph_formula_diffh <- as.formula(paste0(deparse(survival[[2]]),'~',
                                                 paste0(paste0(surv_vars),collapse='+'),'+', 'strata(node)'))
        coxph_formula_diffs <- as.formula(paste0(deparse(survival[[2]]),'~',
                                                 paste0(paste0(surv_vars,'*node'),collapse='+')))
    }

    coxph_model_diffh_diffs <- tryCatch(coxph(coxph_formula_diffh_diffs, data, model=TRUE), error=function(e) NULL)
    coxph_model_diffh <- tryCatch(coxph(coxph_formula_diffh, data, model=TRUE), error=function(e) NULL)
    coxph_model_diffs <- tryCatch(coxph(coxph_formula_diffs, data, model=TRUE), error=function(e) NULL)


    return(list(tree=tree, lmm_model=lmm_model, 
           coxph_model_diffh_diffs=coxph_model_diffh_diffs,coxph_model_diffh=coxph_model_diffh,coxph_model_diffs=coxph_model_diffs))

}


