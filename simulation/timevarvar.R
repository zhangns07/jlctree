# The main model.
# All X1 - X5 will change.
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
                    make_option(c("-g", "--alg"), type="character", default='jlctree', 
                                help="algorithm: jlctree or jlcmm"),
                    make_option(c("-s", "--stop_thre"), type="numeric", default=NULL, help=""),
                    make_option(c("-t", "--test"), type="character", default=NULL, help="Test statistics, rsq, lrt or wald."),
                    make_option(c("-i", "--inter"), type="logical", default=NULL, help="Whether to use interaction term in classmb"),
                    make_option(c("-x", "--continuous"), type="logical", default=TRUE, help="Whether the predictors X1, X2 are continuous"),
                    make_option(c("-m", "--majprob"), type="numeric", default=NULL, help="Maximum probablity for majority family")
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

PARMS <- get_parms(FLAGS$dist); parms <- PARMS$parms; slopes <- PARMS$slopes; lam_D <- PARMS$lam_D; 
Nsim <- FLAGS$maxsim-FLAGS$minsim+1

opt2 <- FLAGS; opt2$help <- NULL; opt2$outdir<- NULL; 
RETbasefilename <- paste0(paste0(names(opt2),"_",opt2),collapse="_")

opt3 <- FLAGS; opt3$help <- NULL; opt3$outdir<- NULL; opt3$minsim<-NULL; opt3$maxsim<-NULL
Rbasefilename <- paste0(paste0(names(opt3),"_",opt3),collapse="_")

if(is.null(FLAGS$continuous)){FLAGS$continuous <- FALSE}
if (FLAGS$alg == 'jlcmm'){
    RET <- matrix(0,ncol=8,nrow=Nsim)
    colnames(RET) <- c('sim','runtime','bestng','B2','B3','B4','B5','B6')
} else if (FLAGS$alg == 'jlctree'){
    if (FLAGS$stop_thre==-1){ RET <- matrix(0,ncol=5,nrow=Nsim)
    } else { RET <- matrix(0,ncol=5,nrow=Nsim*5) }
    colnames(RET) <- c('sim','k','runtime','nsplit','nnode')
}
RET_iter <- 1

for (sim in c(FLAGS$minsim:FLAGS$maxsim)){
    set.seed(sim)

    DATA <- gen_data_timevar(FLAGS, PARMS,seed=sim,survvar=TRUE)
    data <- DATA$data; pseudo_g <- DATA$pseudo_g

    DATA_TEST <- gen_data_timevar(FLAGS,PARMS,seed=sim+623, survvar=TRUE)
    data_test <- DATA_TEST$data; pseudo_g_test <- DATA_TEST$pseudo_g

    if (FLAGS$alg == 'jlcmm'){
        m1 <- Jointlcmm(fixed=y~X1+X2+X3+X4+X5,
                        subject='ID',
                        survival = Surv(time_L,time_Y,delta)~X3+X4+X5,
                        hazard="Weibull",hazardtype="Specific",ng=1,data=data)

        if(FLAGS$inter){ classmb <- ~X1*X2+X3+X4+X5 } else { classmb <- ~X1+X2+X3+X4+X5 }

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
        }
        tok <- Sys.time()

        best_ng <- which.min(BICs)
        Rfilename <- paste0(FLAGS$outdir,'/simret_timevar_RData/',Rbasefilename,'_ng_',best_ng,'_sim_',sim,'.RData')
        save(list=paste0('m',best_ng), file=Rfilename)

        mod <- get(paste0('m',best_ng))
        runtime <- round(difftime(tok,tik,units='secs'),4)
        RET[RET_iter,] <- c(sim, runtime, best_ng, initB[2:6])
        #eval_lcmm_pred_inout_timevar(data,data_test, FLAGS$dist, PARMS$slopes, PARMS$parms, mod, pseudo_g, pseudo_g_test)

        RET_iter <- RET_iter+1
    } else if (FLAGS$alg == 'jlctree'){

        if (FLAGS$stop_thre==-1){
            RET[RET_iter,] <- c(sim,k=-Inf, runtime=0, nsplit=-1, nnode=1)
            #eval_tree_pred_inout_timevar(data,data_test,FLAGS$dist, PARMS$slopes, PARMS$parms, rep(1,nrow(data)), rep(1,nrow(data_test)), pseudo_g, pseudo_g_test)
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
            RET[RET_iter,] <- c(sim,k=-Inf, runtime, nsplit, nnode)
            #idx <- cond_ind_tree$where; idx_test <- predict_class(cond_ind_tree, data_test); eval_tree_pred_inout_timevar(data,data_test,FLAGS$dist, PARMS$slopes, PARMS$parms, idx, idx_test, pseudo_g, pseudo_g_test)
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
                    RET[RET_iter,] <- c(sim,k=kse, runtime, nsplit_prune, nnode_prune)
                }

                #idx <- cond_ind_tree_prune$where; idx_test <- predict_class(cond_ind_tree_prune, data_test); eval_tree_pred_inout_timevar(data,data_test,FLAGS$dist, PARMS$slopes, PARMS$parms, idx, idx_test, pseudo_g, pseudo_g_test)

                RET_iter <- RET_iter+1

            }
            Rfilename <- paste0(FLAGS$outdir,'/simret_timevar_RData/',Rbasefilename,'_sim_',sim,'.RData')
            save(cond_ind_tree, cptable, file=Rfilename)
        }
    }

    cat('sim: ',sim,'\n')

    filename <- paste0(FLAGS$outdir,'/simret_timevar/',RETbasefilename,'.csv')
    write.table(RET, file=filename, sep=',',col.names=TRUE, quote=FALSE)
}


