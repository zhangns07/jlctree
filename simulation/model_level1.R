# model1: time varying, section 5 in biostats supp
# 1. y ~ (1|id) + noise
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
beta_y <- rep(0,5)
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
opt2$betay <- 0
opt2$model <- 'timevary'
basefilename <- paste0(paste0(names(opt2),"_",opt2),collapse="_")

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
    # biomarker
    ranef <- rnorm(Nsub, sd=sd_ranef) 

    all_Z <- switch(FLAGS$type,'I'=c(0,1), 'II'=c(0,1,0,1))
    LTRC_data <- ldply(array(seq_len(Nsub)),function(i){
                           num_i <- num_pseudo[i]
                           tmp_time <- pmin(c(0,times[i,],Inf),time_Y[i])[1:(num_i+1)]
                           ret <- cbind(sub=i, y = rep(ranef[i],num_i), 
                                        X1=rep(X1[i],num_i),X2=all_Z[1:num_i],X3=rep(X3[i],num_i),
                                        time_L=tmp_time[1:num_i],time_Y=tmp_time[-1],
                                        delta=c(rep(0,num_i-1),delta[i]))})

    Nsub_pseudo <- nrow(LTRC_data)
    LTRC_data$y <- LTRC_data$y  + rnorm(Nsub_pseudo, sd=sd_e)
    LTRC_data$X4 <- round(runif(Nsub_pseudo),1)
    LTRC_data$X5 <- sample(c(1:5),Nsub_pseudo,replace=TRUE)

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
        filename <- paste0('simret_model1/',basefilename,'_s_',stop_thre,'_k_',kse,'.csv')
    write.table(cbind(nsplit_RET, nsplit_prune_RET[,(kse+1)]), file=filename, sep=',',col.names=TRUE, quote=FALSE)
}

filename <- paste0('simret_model1/',basefilename,'_root.csv')
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
    cond_ind_tree <- rpart(model.frame
                           (cbind(time_L,time_Y,delta,y,X1,X2,X3,X4,X5)~X1+X2+X3+X4+X5,
                            data = data), control=rpart.control(minsplit=5),
                           method=survlist, 
                           parms=list(LTRC=1, test_stat=FLAGS$test, stop_thre=stop_thre, min_nevent=6))

data <- data.table(data)
y <- data[,list(time_L,time_Y,delta,y,X1,X2,X3,X4,X5)]
x <- data[,X4]
xorder <- order(x);y<-y[xorder,];x<-x[xorder]
wt <- rep(1,nrow(data))
parms=list(LTRC=1, test_stat='lrt', stop_thre=0, min_nevent=6)
continuous <- FALSE

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
