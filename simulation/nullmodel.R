# Move k loop inside sim
library(optparse)
library(data.table)
source('../0.init.R')
option_list <- list(
                    make_option(c("-N", "--Nsub"), type="numeric", default=50, help=""),
                    make_option(c("-c", "--censor"), type="numeric", default=0, help=""),
                    make_option(c("-d", "--dist"), type="character", default='exponential', help=""),
                    make_option(c("-s", "--stop_thre"), type="numeric", default=0, help=""),
                    make_option(c("-l", "--left"), type="logical", default=TRUE, help="Wether left truncated"),
                    make_option(c("-t", "--test"), type="character", default='lrt', help="Test statistics, rsq, lrt or wald.")
                    )

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

FLAGS <- opt

Nsub <- FLAGS$Nsub
censor_rate <- FLAGS$censor
dist <- FLAGS$dist
stop_thre <- FLAGS$stop_thre

gen_survival <- function(n,
 dist=c('exponential','weibull','lognormal')[1], 
 parms=list(rate=0.25, shape=0.9, scale=2, meanlog=1.4, sdlog=0.4)){
    if (dist=='exponential'){
        time_T <- rexp(n, rate=parms$rate)
    } else if (dist=='weibull'){
        time_T <- rweibull(n, shape=parms$shape, scale=parms$scale)
    } else if (dist=='lognormal'){
        time_T <- rlnorm(n, meanlog=parms$meanlog, sdlog=parms$sdlog)
    }
    return (time_T)
}
sd_ranef <- 2
sd_e <- 1
beta <- rep(0, 5)
Nsim <- 200

survs <- survs_v3
survlist <- list(eval=surve, split=survs, init=survi)

basefilename <- paste0('Nsub_',Nsub, '_dist_',dist, '_censor_',censor_rate,
                       '_beta0_ranef_',sd_ranef,'_noise_',sd_e)


nsplit_RET <- rep(0, Nsim)
which_RET <- rep(0, Nsim)
#nsplit_prune_RET <- rep(0,Nsim)
nsplit_prune_RET <- matrix(0,nrow=Nsim,ncol=4)
which_prune_RET <- matrix(0,nrow=Nsim,ncol=4)

for (sim in c(1:Nsim)){
    set.seed(sim)

    # generate time
    if (FLAGS$left){
        time_L <- runif(2*Nsub, min=0, max=2)
    } else {
        time_L <- rep(0, 2*Nsub)
    }
    time_T <- gen_survival(2*Nsub, dist)


    # only keep those where L < T
    time_tokeep <- time_L < time_T
    time_L <- time_L[time_tokeep][1:Nsub]
    time_T <- time_T[time_tokeep][1:Nsub]

    # censor time 
    if (censor_rate==0){
        time_C <- time_T
    } else{
        time_C <- time_L + rexp(Nsub, rate=censor_rate)
    }

    # (L,Y,delta) triplet
    delta <- as.numeric(time_C >= time_T)
    time_Y <- pmin(time_T, time_C)

    # Xs
    X1 <- runif(Nsub, min=1, max=2)
    X2 <- runif(Nsub, min=1, max=2)
    X3 <- round(runif(Nsub),1)
    X4 <- as.numeric(runif(Nsub)>0.5)
    X5 <- as.numeric(runif(Nsub)>0.5)
    X <- cbind(X1, X2, X3, X4,X5)

    # biomarker
    y <-  c(X %*% beta + rnorm(Nsub, sd=sd_ranef) + rnorm(Nsub, sd=sd_e))
    data <- data.frame(cbind(time_L, time_Y, delta, X, y))

    if (FLAGS$left){
        cond_ind_tree <- rpart(model.frame
                               (cbind(time_L, time_Y, delta, y, X1,X2,X3,X4,X5) ~ X1+X2+X3+X4+X5,
                                data = data), control=rpart.control(minsplit=5),
                               method=survlist, 
                               parms=list(LTRC=1, test_stat=FLAGS$test, stop_thre=stop_thre, min_nevent=10))
    } else {
        cond_ind_tree <- rpart(model.frame
                               (cbind(time_Y, delta, y, X1,X2,X3,X4,X5) ~ X1+X2+X3+X4+X5,
                                data = data), control=rpart.control(minsplit=5),
                               method=survlist, 
                               parms=list(LTRC=0, test_stat=FLAGS$test, stop_thre=stop_thre, min_nevent=10))

    }

    nsplit <- max(cond_ind_tree$cptable[,'nsplit'])
    nsplit_RET[sim] <- nsplit
    if (nsplit==0){ whichsplit <- 0 } else { whichsplit <- as.numeric(gsub('X','',cond_ind_tree$frame[1,'var'])) }
    which_RET[sim] <- whichsplit

    bo<-0
    while(bo!=10){
        xfit <- tryCatch({xpred.rpart(cond_ind_tree, xval=5, cp=cond_ind_tree$cptable[,'CP'])},error=function(e){NA})
        if (any(is.na(xfit))){
            bo<-bo+1
        } else break
    }

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

        if (nsplit_prune==0){ whichsplit_prune <- 0 } else { 
            whichsplit_prune <- as.numeric(gsub('X','',cond_ind_tree_prune$frame[1,'var'])) }
        which_prune_RET[sim] <- whichsplit_prune
    }

   cat('sim: ',sim,'\n')

}

#for (kse in c(0:3)){
#    if (FLAGS$left){
#        filename <- paste0('simret/',basefilename,'_s_',stop_thre,'_k_',kse,'_test_',FLAGS$test,'.csv')
#    } else {
#        filename <- paste0('simret/',basefilename,'_left_FALSE_s_',stop_thre,'_k_',kse,'_test_',FLAGS$test,'.csv')
#    }
#    write.table(cbind(nsplit_RET, nsplit_prune_RET[,(kse+1)]), file=filename, sep=',',col.names=TRUE, quote=FALSE)
#}

data <- data.table(data)
y <- data[,list(time_L, time_Y, delta, y, X1,X2,X3,X4,X5)]
x <- data$X5
xorder <- order(x);y <- y[xorder,]; x <- x[xorder]
parms<-list(min_nevent=1, stop_thre=FLAGS$stop_thre, test_stat=FLAGS$test)

y = data.frame(y)
colnames(y)[1:4] <- c('start','end','event','biomarker')
formulay1 <- Surv(start, end, event) ~ . - biomarker
formulay2 <- Surv(start, end, event) ~ .
nevents <- sum(y[,'event'])
root_val <- get_node_val(formulay1, formulay2, y, test_stat= parms$test_stat)
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
                cat(left_val, right_val,'\n')
                get_split_utility(root_val, left_val, right_val, parms$test_stat)
            }, error = function(e){ c(0, sign(1))})#, warning = function(w){c(0,sign(0))})
        }
        goodness[i] = result[1]; direction[i] = result[2]
    }
}

