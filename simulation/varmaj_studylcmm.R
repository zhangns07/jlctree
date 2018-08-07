# The main model.
library(optparse)
library(data.table)
library(plyr)
library(lcmm)
source('0.init.R')
source('0.gen_survival.R')
gen_data <- function(FLAGS, PARMS, seed){

    Nsub <- FLAGS$Nsub
    censor_rate <- FLAGS$censor
    dist <- FLAGS$dist
    parms <- PARMS$parms; slopes <- PARMS$slopes; lam_D <- PARMS$lam_D
    sd_ranef <- 0.2
    sd_e <- 0.1

    set.seed(seed)
    if(FLAGS$numpred>= 5){
        # X1 - X5
        Xg1 <- round(runif(2*Nsub),2)
        Xg2 <- round(runif(2*Nsub),2)

        Xo1 <- as.numeric(runif(2*Nsub)>0.5)
        Xo2 <- round(runif(2*Nsub),1)
        Xo3 <- sample(c(1:5),2*Nsub,replace=TRUE)

        Xg <- cbind(Xg1,Xg2)
        Xo <- cbind(Xo1,Xo2,Xo3)
    } 
    if (FLAGS$numpred>= 10){
        Xg3 <- round(runif(2*Nsub),2)
        Xg4 <- round(runif(2*Nsub),2)

        Xo4 <- as.numeric(runif(2*Nsub)>0.5)
        Xo5 <- round(runif(2*Nsub),1)
        Xo6 <- sample(c(1:5),2*Nsub,replace=TRUE)

        Xg <- cbind(Xg,Xg3,Xg4)
        Xo <- cbind(Xo,Xo4,Xo5,Xo6)
        slopes <- cbind(slopes,slopes)
    } 
    if (FLAGS$numpred>=20){
        Xg5 <- round(runif(2*Nsub),2)
        Xg6 <- round(runif(2*Nsub),2)
        Xg7 <- round(runif(2*Nsub),2)
        Xg8 <- round(runif(2*Nsub),2)

        Xo7<- as.numeric(runif(2*Nsub)>0.5)
        Xo8 <- round(runif(2*Nsub),1)
        Xo9 <- sample(c(1:5),2*Nsub,replace=TRUE)
        Xo10 <- as.numeric(runif(2*Nsub)>0.5)
        Xo11 <- round(runif(2*Nsub),1)
        Xo12 <- sample(c(1:5),2*Nsub,replace=TRUE)

        Xg <- cbind(Xg,Xg5,Xg6,Xg7,Xg8)
        Xo <- cbind(Xo,Xo7,Xo8,Xo9,Xo10,Xo11,Xo12)
        slopes <- cbind(slopes,slopes)
    }

    X <- cbind(Xg,Xo)
    # change slopes when there are more predictors
    slopes <- slopes/ (FLAGS$numpred/5)
    g <- get_latent_class(Xg,FLAGS$struct, FLAGS$member, seed=seed, FLAGS$majprob)

    if(FLAGS$dist == 'lognormal'){ 
        ebx <- rep(1 , 2*Nsub)
        tmp_parms <- parms[g,]
    } else{
        ebx <- exp(rowSums(slopes[g,] * Xo))
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
                           tmp_time_L <- rep(time_L[i],num_i)
                           tmp_time_Y <- rep(time_Y[i],num_i)
                           tmp_delta <- rep(delta[i],num_i)

                           ret <- cbind(ID=i, 
                                        do.call("rbind", replicate(num_i, Xg[i,], simplify = FALSE)),
                                        do.call("rbind", replicate(num_i, Xo[i,], simplify = FALSE)),
                                        time_L=tmp_time_L, time_Y=tmp_time_Y, delta=tmp_delta)})

    ranef <- rnorm(Nsub, sd=sd_ranef)
    ranefs <- ranef[LTRC_data$ID]
    fixef <- c(0,1,1,2)
    pseudo_g <- g[LTRC_data$ID]
    y <- fixef[pseudo_g] + ranefs + rnorm(nrow(LTRC_data), sd=sd_e) 
    data <- cbind(LTRC_data,y)

    return(list(data=data,pseudo_g=pseudo_g))
}


get_latent_class <- function(X,struct,member,seed=0,majprob=NULL){

    if(is.null(majprob)){
        Wmul <- 1
    } else {
        Wlist <- switch(struct,
                        tree=list('0.25' = 0, '0.5'=1, '0.7'=2, '0.85'=4.5, '0.99'=50), 
                        linear=list('0.25' = 0, '0.5'=1.5, '0.7'=3.8, '0.85'=9, '0.99'=150))
        Wmul <- Wlist[[as.character(majprob)]]
    }

    set.seed(seed)
    W <- matrix(rnorm(ncol(X) * 4),byrow=TRUE,nrow=4)
    W <- W/ apply(W,1,function(x){sqrt(sum(x^2))})

    if (member == 'partition' ){
        g <- apply(X,1,function(x){
                       weights <- W %*% (x-0.5)
                       which.max(weights)})
    } else if (member == 'multinomial'){
        W <- W*Wmul
        g <- apply(X,1,function(x){
                       weights <- exp(W %*% (x-0.5))
                       weights <- weights/sum(weights)
                       cumweights <- cumsum(weights)
                       which(runif(1)  < cumweights)[1] })
    }
    return (g)
}


option_list <- list(make_option(c("-N", "--Nsub"), type="numeric", default=100, help=""),
                    make_option(c("-c", "--censor"), type="numeric", default=0, help=""),
                    make_option(c("-d", "--dist"), type="character", default='exponential', help=""),
                    make_option(c("-p", "--numpred"), type="numeric", default=5, help="number of predictors"),
                    make_option(c("-o", "--outdir"), type="character", default='.', help="Output directory."),
                    make_option(c("-l", "--minsim"), type="numeric", default=1, help="min sim"),
                    make_option(c("-u", "--maxsim"), type="numeric", default=20, help="max sim"),
                    make_option(c("-m", "--majprob"), type="numeric", default=0.85, help="Maximum probablity for majority family")
                    )

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);
FLAGS <- opt


PARMS <- get_parms(FLAGS$dist); parms <- PARMS$parms; slopes <- PARMS$slopes; lam_D <- PARMS$lam_D; 
Nsim <- FLAGS$maxsim-FLAGS$minsim+1

opt2 <- FLAGS; opt2$help <- NULL; opt2$outdir<- NULL; 
RETbasefilename <- paste0(paste0(names(opt2),"_",opt2),collapse="_")

opt3 <- FLAGS; opt3$help <- NULL; opt3$outdir<- NULL; opt3$minsim<-NULL; opt3$maxsim<-NULL
Rbasefilename <- paste0(paste0(names(opt3),"_",opt3),collapse="_")

FLAGS$continuous <- FALSE
FLAGS$inter <- FALSE
FLAGS$test <- 'lrt'
FLAGS$member <- 'multinomial'
FLAGS$struct <- 'linear'
FLAGS$alg <- 'jlcmm'


RET <- matrix(0,ncol=7,nrow=Nsim)
colnames(RET) <- c('sim','t1','t2','t3','t4','t5','t6')
RET_iter <- 1

for (sim in c(FLAGS$minsim:FLAGS$maxsim)){
    set.seed(sim)

    DATA <- gen_data(FLAGS, PARMS,seed=sim)
    data <- DATA$data; pseudo_g <- DATA$pseudo_g

    runtimes <- rep(0, 6)
    BICs <- rep(Inf, 6)

    if(FLAGS$numpred==5){
        tik <- Sys.time()
        m1<-Jointlcmm(fixed=y~Xg1+Xg2+Xo1+Xo2+Xo3,
                      subject='ID',
                      survival = Surv(time_L, time_Y, delta) ~ Xo1 + Xo2 + Xo3,
                      hazard="Weibull",hazardtype="Specific",ng=1,data=data)
        tok <- Sys.time()
        runtimes[1] <- round(difftime(tok,tik,units='secs'),4)
        classmb <- ~Xg1+Xg2

        for(ng in c(2:6)){
            tik <- Sys.time()
            mod <- Jointlcmm(fixed=y~Xg1+Xg2+Xo1+Xo2+Xo3,
                             classmb=classmb,subject='ID',mixture=~1,
                             survival=Surv(time_L,time_Y,delta)~mixture(Xo1+Xo2+Xo3),
                             hazard="Weibull",hazardtype="Specific",ng=ng,data=data,B=m1)

            if (mod$conv %in% c(4,5)){
                mod <- Jointlcmm(fixed=y~Xg1+Xg2+Xo1+Xo2+Xo3,
                                 classmb=classmb,subject='ID',mixture=~1,
                                 survival=Surv(time_L,time_Y,delta)~mixture(Xo1+Xo2+Xo3),
                                 hazard="Weibull",hazardtype="Specific",ng=ng,data=data)
            }

            if (mod$conv %in% c(4,5)){
                mod <- Jointlcmm(fixed=y~~Xg1+Xg2+Xo1+Xo2+Xo3,
                                 classmb=classmb,subject='ID',mixture=~1,
                                 survival=Surv(time_L,time_Y,delta)~mixture(Xo1+Xo2+Xo3),
                                 hazard="Weibull",hazardtype="Specific",ng=ng,data=data,B=random(m1))
            }

            tok <- Sys.time()
            runtimes[ng] <- round(difftime(tok,tik,units='secs'),4)

            assign(paste0('m',ng), mod)
            if (!mod$conv %in% c(4,5)){ BICs[ng] <- mod$BIC } 
        }

    } else if (FLAGS$numpred==10){

        tik <- Sys.time()
        m1 <- Jointlcmm(fixed=y~Xg1+Xg2+Xg3+Xg4+Xo1+Xo2+Xo3+Xo4+Xo5+Xo6,
                        subject='ID',
                        survival=Surv(time_L,time_Y,delta)~Xo1+Xo2+Xo3+Xo4+Xo5+Xo6,
                        hazard="Weibull",hazardtype="Specific",ng=1,data=data)
        tok <- Sys.time()
        runtimes[1] <- round(difftime(tok,tik,units='secs'),4)
        classmb <- ~Xg1+Xg2+Xg3+Xg4

        for(ng in c(2:6)){
            tik <- Sys.time()
            mod <- Jointlcmm(fixed=y~Xg1+Xg2+Xg3+Xg4+Xo1+Xo2+Xo3+Xo4+Xo5+Xo6,
                             classmb=classmb,subject='ID',mixture=~1,
                             survival=Surv(time_L,time_Y,delta)~mixture(Xo1+Xo2+Xo3+Xo4+Xo5+Xo6),
                             hazard="Weibull",hazardtype="Specific",ng=ng,data=data,B=m1)

            if (mod$conv %in% c(4,5)){
                mod <- Jointlcmm(fixed=y~Xg1+Xg2+Xg3+Xg4+Xo1+Xo2+Xo3+Xo4+Xo5+Xo6,
                                 classmb=classmb,subject='ID',mixture=~1,
                                 survival=Surv(time_L,time_Y,delta)~mixture(Xo1+Xo2+Xo3+Xo4+Xo5+Xo6),
                                 hazard="Weibull",hazardtype="Specific",ng=ng,data=data)
            }

            if (mod$conv %in% c(4,5)){
                mod <- Jointlcmm(fixed=y~Xg1+Xg2+Xg3+Xg4+Xo1+Xo2+Xo3+Xo4+Xo5+Xo6,
                                 classmb=classmb,subject='ID',mixture=~1,
                                 survival=Surv(time_L,time_Y,delta)~mixture(Xo1+Xo2+Xo3+Xo4+Xo5+Xo6),
                                 hazard="Weibull",hazardtype="Specific",ng=ng,data=data,B=random(m1))
            }

            tok <- Sys.time()
            runtimes[ng] <- round(difftime(tok,tik,units='secs'),4)

            assign(paste0('m',ng), mod)
            if (!mod$conv %in% c(4,5)){ BICs[ng] <- mod$BIC } 
        }
    } else if (FLAGS$numpred==20){
        tik <- Sys.time()
        m1 <- Jointlcmm(fixed= y~Xg1+Xg2+Xg3+Xg4+Xg5+Xg6+Xg7+Xg8+Xo1+Xo2+Xo3+Xo4+Xo5+Xo6+Xo7+Xo8+Xo9+Xo10+Xo11+Xo12,
                        subject='ID',
                        survival = Surv(time_L,time_Y,delta)~Xo1+Xo2+Xo3+Xo4+Xo5+Xo6+Xo7+Xo8+Xo9+Xo10+Xo11+Xo12,
                        hazard="Weibull",hazardtype="Specific",ng=1,data=data)
        tok <- Sys.time()
        runtimes[1] <- round(difftime(tok,tik,units='secs'),4)
        classmb <- ~Xg1+Xg2+Xg3+Xg4+Xg5+Xg6

        for(ng in c(2:6)){
            tik <- Sys.time()
            mod <- Jointlcmm(fixed= y~Xg1+Xg2+Xg3+Xg4+Xg5+Xg6+Xg7+Xg8+Xo1+Xo2+Xo3+Xo4+Xo5+Xo6+Xo7+Xo8+Xo9+Xo10+Xo11+Xo12,
                             classmb=classmb,subject='ID',mixture=~1,
                             survival=Surv(time_L,time_Y,delta)~mixture(Xo1+Xo2+Xo3+Xo4+Xo5+Xo6+Xo7+Xo8+Xo9+Xo10+Xo11+Xo12),
                             hazard="Weibull",hazardtype="Specific",ng=ng,data=data,B=m1)

            if (mod$conv %in% c(4,5)){
                mod <- Jointlcmm(fixed= y~Xg1+Xg2+Xg3+Xg4+Xg5+Xg6+Xg7+Xg8+Xo1+Xo2+Xo3+Xo4+Xo5+Xo6+Xo7+Xo8+Xo9+Xo10+Xo11+Xo12,
                                 classmb=classmb,subject='ID',mixture=~1,
                                 survival = Surv(time_L,time_Y,delta)~mixture(Xo1+Xo2+Xo3+Xo4+Xo5+Xo6+Xo7+Xo8+Xo9+Xo10+Xo11+Xo12),
                                 hazard="Weibull",hazardtype="Specific",ng=ng,data=data)
            }

            if (mod$conv %in% c(4,5)){
                mod <- Jointlcmm(fixed= y~Xg1+Xg2+Xg3+Xg4+Xg5+Xg6+Xg7+Xg8+Xo1+Xo2+Xo3+Xo4+Xo5+Xo6+Xo7+Xo8+Xo9+Xo10+Xo11+Xo12,
                                 classmb=classmb,subject='ID',mixture=~1,
                                 survival = Surv(time_L,time_Y,delta)~mixture(Xo1+Xo2+Xo3+Xo4+Xo5+Xo6+Xo7+Xo8+Xo9+Xo10+Xo11+Xo12),
                                 hazard="Weibull",hazardtype="Specific",ng=ng,data=data,B=random(m1))
            }

            tok <- Sys.time()
            runtimes[ng] <- round(difftime(tok,tik,units='secs'),4)

            assign(paste0('m',ng), mod)
            if (!mod$conv %in% c(4,5)){ BICs[ng] <- mod$BIC } 
        }
    }


    best_ng <- which.min(BICs)
    mod <- get(paste0('m',best_ng))
    Rfilename <- paste0(FLAGS$outdir,'/',Rbasefilename,'_ng_',best_ng,'_sim_',sim,'.RData')
    save(list=paste0('m',best_ng), file=Rfilename)

    RET[RET_iter,] <- c(sim, runtimes)
    RET_iter <- RET_iter+1
    cat('sim: ',sim,'\n')

    filename <- paste0(FLAGS$outdir,'/',RETbasefilename,'.csv')
    write.table(RET, file=filename, sep=',',col.names=TRUE, quote=FALSE)
}




