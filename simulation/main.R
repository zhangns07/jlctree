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
source('0.init.R')
source('0.gen_survival.R')
option_list <- list(make_option(c("-N", "--Nsub"), type="numeric", default=100, help=""),
                    make_option(c("-c", "--censor"), type="numeric", default=0, help=""),
                    make_option(c("-d", "--dist"), type="character", default='exponential', help=""),
                    make_option(c("-o", "--outdir"), type="character", default='.', help="Output directory."),
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
                    make_option(c("-i", "--inter"), type="logical", default=NULL, help="Whether to use interaction term in classmb"),
                    make_option(c("-x", "--continuous"), type="logical", default=FALSE, help="Whether the predictors X1, X2 are continuous")
                    )


opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

FLAGS <- opt
if (FLAGS$alg == 'jlctree'){
#    FLAGS$stop_thre=3.84; FLAGS$test='lrt'
    if (is.null(FLAGS$stop_thre) | is.null(FLAGS$test)){
        stop('Missing argument: stop_thre and test.')
    }
}

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

if (FLAGS$alg == 'jlcmm'){
    RET <- matrix(0,ncol=12,nrow=Nsim)
    colnames(RET) <- c('sim','runtime','bestng','B2','B3','B4','B5','B6','ISE','MSE_b','MSE_y','purity')
} else if (FLAGS$alg == 'jlctree'){
    if (FLAGS$stop_thre==-1){ RET <- matrix(0,ncol=9,nrow=Nsim)
    } else { RET <- matrix(0,ncol=9,nrow=Nsim*5) }
    colnames(RET) <- c('sim','k','runtime','nsplit','nnode','ISE','MSE_b','MSE_y','purity')
}
RET_iter <- 1

for (sim in c(FLAGS$minsim:FLAGS$maxsim)){
    set.seed(sim)

    # X1 - X5
    if (FLAGS$continuous){
        if(FLAGS$struct == 'linear'){
            X1 <- round(runif(2*Nsub, min=1,max=3),1)
            X2 <- round(runif(2*Nsub, min=1,max=3),1)
        } else {
            X1 <- round(runif(2*Nsub),1)
            X2 <- round(runif(2*Nsub),1)
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

    if (FLAGS$alg == 'jlcmm'){
        m1 <- Jointlcmm(fixed=y~X1+X2+X3+X4+X5,
                        subject='ID',
                        survival = Surv(time_L,time_Y,delta)~X3+X4+X5,
                        hazard="Weibull",hazardtype="Specific",ng=1,data=data)

        if(FLAGS$inter){ classmb <- ~X1*X2+X3+X4+X5 } else { classmb <- ~X1+X2+X3+X4+X5}

        BICs <- rep(Inf, 6)
        initB <- rep(1,6) #1:B=m1; 2:B=NULL; 3:B=random(m1); 4:none worked
        tik <- Sys.time()
        for(ng in c(2:6)){
            mod <- Jointlcmm(fixed=y~X1+X2+X3+X4+X5,mixture=~1,
                             classmb=classmb,subject='ID',
                             survival = Surv(time_L,time_Y,delta)~mixture(X3+X4+X5),
                             hazard="Weibull",hazardtype="Specific",ng=ng,data=data,B=m1)

            if (mod$conv %in% c(4,5)){
                mod <- Jointlcmm(fixed=y~X1+X2+X3+X4+X5,mixture=~1,
                                 classmb=classmb,subject='ID',
                                 survival = Surv(time_L,time_Y,delta)~mixture(X3+X4+X5),
                                 hazard="Weibull",hazardtype="Specific",ng=ng,data=data)
                initB[ng] <- 2
            }

            if (mod$conv %in% c(4,5)){
                mod <- Jointlcmm(fixed=y~X1+X2+X3+X4+X5,mixture=~1,
                                 classmb=classmb,subject='ID',
                                 survival = Surv(time_L,time_Y,delta)~mixture(X3+X4+X5),
                                 hazard="Weibull",hazardtype="Specific",ng=ng,data=data,B=random(m1))
                initB[ng] <- 3
            }

            if (!mod$conv %in% c(4,5)){
                BICs[ng] <- mod$BIC
            } else {
                initB[ng] <- 4
            }

            assign(paste0('m',ng), mod)
            Rfilename <- paste0(FLAGS$outdir,'/simret_main_RData/',Rbasefilename,'_ng_',ng,'_sim_',sim,'.RData')
            save(list=paste0('m',ng), file=Rfilename)
        }
        tok <- Sys.time()

        best_ng <- which.min(BICs)
        mod <- get(paste0('m',best_ng))
        runtime <- round(difftime(tok,tik,units='secs'),4)
        RET[RET_iter,] <- c(sim, runtime, best_ng, initB[2:6], 
                            eval_lcmm_pred(data,dist,slopes,parms,mod,g=pseudo_g,inter=FLAGS$inter))
        RET_iter <- RET_iter+1
    } else if (FLAGS$alg == 'jlctree'){

        if (FLAGS$stop_thre==-1){
            RET[RET_iter,] <- c(sim,k=-Inf, runtime=0, nsplit=-1, nnode=1,
                                eval_tree_pred(data,dist, slopes, parms, rep(1,nrow(data)),g=pseudo_g))
            RET_iter <- RET_iter+1
        } else {
            survs <- survs_v3
            survlist <- list(eval=surve, split=survs, init=survi)

            tik <- Sys.time()
            cond_ind_tree <- rpart(model.frame
                                   (cbind(time_L,time_Y,delta,y,X3,X4,X5)~X1+X2+X3+X4+X5, data = data), 
                                   control=rpart.control(minsplit=5),
                                   method=survlist, 
                                   parms=list(LTRC=1, test_stat=FLAGS$test, 
                                              stop_thre=FLAGS$stop_thre, min_nevent=4,stable=TRUE))

            xfit <- xpred.rpart(cond_ind_tree, xval=5)
            xerror <- apply(xfit,2,mean)
            xstd <- apply(xfit,2,sd)/sqrt(nrow(data))
            cptable <- cbind(cond_ind_tree$cptable, xerror,xstd)

            cventry <- which.min(cptable[, "xerror"])
            xerrorcv <- cptable[cventry, "xerror"]
            tok <- Sys.time()
            runtime <- round(difftime(tok,tik,units='secs'),4)

            nsplit <- max(cond_ind_tree$cptable[,'nsplit'])
            nnode <- sum(grepl('leaf',cond_ind_tree$frame$var))
            RET[RET_iter,] <- c(sim,k=-Inf, runtime, nsplit, nnode,
                                eval_tree_pred(data,dist, slopes, parms, cond_ind_tree$where,g=pseudo_g))
            RET_iter <- RET_iter+1

            for (kse in c(0:3)){
                sexerrorcv <- xerrorcv + kse*cptable[cventry, "xstd"] 
                cpcvse <- cptable[which.max(cptable[, "xerror"] <= sexerrorcv), "CP"]
                cond_ind_tree_prune <- prune(cond_ind_tree, cp=cpcvse)
                nsplit_prune <- max(cond_ind_tree_prune$cptable[,'nsplit'])
                nnode_prune  <- sum(grepl('leaf',cond_ind_tree_prune$frame$var))

                if (nsplit_prune == nsplit){
                    RET[RET_iter,] <- RET[RET_iter-1,]
                    RET[RET_iter,2] <- kse
                } else {
                    RET[RET_iter,] <- c(sim,k=kse, runtime, nsplit_prune, nnode_prune,
                                        eval_tree_pred(data,dist, slopes, parms, cond_ind_tree_prune$where,g=pseudo_g))
                }
                RET_iter <- RET_iter+1

            }
            Rfilename <- paste0(FLAGS$outdir,'/simret_main_RData/',Rbasefilename,'_sim_',sim,'.RData')
            save(cond_ind_tree, cptable, file=Rfilename)
        }
    }

    cat('sim: ',sim,'\n')

    filename <- paste0(FLAGS$outdir,'/simret_main/',RETbasefilename,'.csv')
    write.table(RET, file=filename, sep=',',col.names=TRUE, quote=FALSE)
}




