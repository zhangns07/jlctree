# The main model.
# Options:
# - data generating: Nsub, censor, dist
# - bookkeeping: outdir, minsim, maxsim
# - latent class structure: structure=c('tree','linear','nonlinear')
# - latent class membership: membership=c('partition','multinomial')
# - alg=c('jlctree','jclmm')
# - jlctree specific: stop_thre, test
library(optparse)
library(data.table)
library(plyr)
library(lcmm)
source('../0.init.R')
source('../0.gen_survival.R')
option_list <- list(make_option(c("-N", "--Nsub"), type="numeric", default=200, help=""),
                    make_option(c("-c", "--censor"), type="numeric", default=0, help=""),
                    make_option(c("-d", "--dist"), type="character", default='lognormal', help=""),
                    make_option(c("-o", "--outdir"), type="character", default='../..', help="Output directory."),
                    make_option(c("-l", "--minsim"), type="numeric", default=1, help="min sim"),
                    make_option(c("-u", "--maxsim"), type="numeric", default=200, help="max sim"),
                    make_option(c("-a", "--struct"), type="character", default='tree',
                                help="latent class generating structure, tree (4 classes), XOR (2 classes), linear, or nonlinear."),
                    make_option(c("-b", "--member"), type="character", default='partition',
                                help="latent class membership, partition or multinomial"),
                    make_option(c("-g", "--alg"), type="character", default='jlctree', 
                                help="algorithm: jlctree or jlcmm"),
                    make_option(c("-s", "--stop_thre"), type="numeric", default=NULL, help=""),
                    make_option(c("-t", "--test"), type="character", default=NULL, help="Test statistics, rsq, lrt or wald."),
                    make_option(c("-i", "--inter"), type="logical", default=NULL, help="Whether to use interaction term in classmb"))


opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

FLAGS <- opt
if (FLAGS$alg == 'jlcmm'){
    if (is.null(FLAGS$inter)){
        stop('Missing argument: inter.')
    }
}

Nsub <- FLAGS$Nsub
censor_rate <- FLAGS$censor
dist <- FLAGS$dist
PARMS <- get_parms(dist); parms <- PARMS$parms; slopes <- PARMS$slopes; lam_D <- PARMS$lam_D; 

sd_ranef <- 0.2
sd_e <- 0.1
Nsim <- FLAGS$maxsim-FLAGS$minsim+1
beta_y <- rep(0,5)

opt2 <- FLAGS; opt2$help <- NULL; opt2$outdir<- NULL; 
RETbasefilename <- paste0(paste0(names(opt2),"_",opt2),collapse="_")

opt3 <- FLAGS; opt3$help <- NULL; opt3$outdir<- NULL; opt3$minsim<-NULL; opt3$maxsim<-NULL
Rbasefilename <- paste0(paste0(names(opt3),"_",opt3),collapse="_")

filename <- paste0(FLAGS$outdir,'/simret_main/',RETbasefilename,'.csv')
RET <- read.table(filename, sep=',')
RET_iter <- 1

for (sim in c(FLAGS$minsim:FLAGS$maxsim)){
    set.seed(sim)

    # X1 - X5
    if(FLAGS$struct == 'linear'){
        X1 <- sample(c(1:3),2*Nsub,replace=TRUE)
        X2 <- sample(c(1:3),2*Nsub,replace=TRUE)
    } else if (FLAGS$struct == 'nonlinear'){
	stop("No need for nonlinear.")
    } else {
        X1 <- as.numeric(runif(2*Nsub)>0.5)
        X2 <- as.numeric(runif(2*Nsub)>0.5)
    }
    X3 <- as.numeric(runif(2*Nsub)>0.5)
    X4 <- round(runif(2*Nsub),1)
    X5 <- sample(c(1:5),2*Nsub,replace=TRUE)
    X <- cbind(X1,X2,X3,X4,X5)

    g <- get_latent_class(X1,X2,FLAGS$struct, FLAGS$member, seed=sim)
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
			   tmp_time_L <- rep(time_L[i],num_i)
			   tmp_time_Y <- rep(time_Y[i],num_i)
			   tmp_delta <- rep(delta[i],num_i)

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

    if (FLAGS$alg == 'jlcmm'){
	    currRET <- RET[RET_iter,]
		    best_ng <- currRET['bestng']

		    Rfilename <- paste0(FLAGS$outdir,'/simret_main_RData/',Rbasefilename,'_ng_',best_ng,'_sim_',sim,'.RData')
		    load(Rfilename)
		    mod <- get(paste0('m',best_ng))
		    tmpeval <- eval_lcmm_pred(data,dist,slopes,parms,mod,g=pseudo_g,inter=FLAGS$inter)

		    RET[RET_iter,'ISE'] <- tmpeval['ISE']
		    RET[RET_iter,'MSE_b'] <- tmpeval['MSE_b']
		    RET[RET_iter,'MSE_y'] <- tmpeval['MSE_y']
		    RET[RET_iter,'purity'] <- tmpeval['purity']
		    RET_iter <- RET_iter+1
    } 

    cat('sim: ',sim,'\n')

    filename <- paste0(FLAGS$outdir,'/simret_main/',RETbasefilename,'.csv')
    write.table(RET, file=filename, sep=',',col.names=TRUE, quote=FALSE)
}

