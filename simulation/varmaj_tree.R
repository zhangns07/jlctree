# The main model.
library(optparse)
library(data.table)
library(plyr)
library(lcmm)
source('0.init.R')
source('0.gen_survival.R')
source('../jlctree.R')
source('../util.R')
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
                    make_option(c("-x", "--continuous"), type="logical", default=TRUE, help="Whether the predictors X1, X2 are continuous"),
                    make_option(c("-m", "--majprob"), type="numeric", default=0.85, help="Maximum probablity for majority family")
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

PARMS <- get_parms(FLAGS$dist); parms <- PARMS$parms; slopes <- PARMS$slopes; lam_D <- PARMS$lam_D; 
Nsim <- FLAGS$maxsim-FLAGS$minsim+1

opt2 <- FLAGS; opt2$help <- NULL; opt2$outdir<- NULL; 
RETbasefilename <- paste0(paste0(names(opt2),"_",opt2),collapse="_")

opt3 <- FLAGS; opt3$help <- NULL; opt3$outdir<- NULL; opt3$minsim<-NULL; opt3$maxsim<-NULL
Rbasefilename <- paste0(paste0(names(opt3),"_",opt3),collapse="_")

if(is.null(FLAGS$continuous)){FLAGS$continuous <- FALSE}
RET <- matrix(0,ncol=18,nrow=Nsim) 
colnames(RET) <- c('sim','k','runtime','nsplit','nnode',
                   "ISE1","ISE2","ISE3", 
                   "MSE_b1","MSE_b2","MSE_b3", 
                   "MSE_y","purity",
                   "ISE_test1","ISE_test2","ISE_test3",
                   "MSE_y_test","purity_test")

RET_iter <- 1

for (sim in c(FLAGS$minsim:FLAGS$maxsim)){
    set.seed(sim)

    DATA <- gen_data(FLAGS, PARMS,seed=sim)
    data <- DATA$data; pseudo_g <- DATA$pseudo_g

    DATA_TEST <- gen_data(FLAGS,PARMS,seed=sim+623)
    data_test <- DATA_TEST$data; pseudo_g_test <- DATA_TEST$pseudo_g

    if (FLAGS$stop_thre==-1){
        tik <- Sys.time()
        tree <- jlctree(survival=Surv(time_L,time_Y,delta)~X3+X4+X5,
                        fixed=y~X1+X2+X3+X4+X5,
                        random=~1,
                        classmb=~X1+X2+X3+X4+X5,
                        subject='ID',maxng=6,data=data,
                        parms=list(test_stat=FLAGS$test, 
                                   stop_thre=Inf,min_nevent=3))
        tok <- Sys.time()
        runtime <- round(difftime(tok,tik,units='secs'),4)
        ret <- c(sim,k=-Inf, runtime=runtime, nsplit=-1, nnode=1)

        idx <- rep(1, nrow(data)); idx_test <- rep(1,nrow(data_test))
        EVALS <- eval_tree_pred_inout_jlctree(data,data_test,FLAGS$dist, PARMS$slopes, PARMS$parms, idx, idx_test, pseudo_g, pseudo_g_test,tree,timevar=FALSE)
        ret <- c(ret,EVALS)

        RET[RET_iter,] <- ret
        RET_iter <- RET_iter+1
    } else {
        tik <- Sys.time()
        tree <- jlctree(survival=Surv(time_L,time_Y,delta)~X3+X4+X5,
                fixed=y~X1+X2+X3+X4+X5,
                random=~1,
                classmb=~X1+X2+X3+X4+X5,
                subject='ID',maxng=6,data=data,
                parms=list(test_stat=FLAGS$test, 
                           stop_thre=FLAGS$stop_thre,min_nevent=3))
        tok <- Sys.time()
        runtime <- round(difftime(tok,tik,units='secs'),4)
        #Rfilename <- paste0(FLAGS$outdir,'/simret_varmaj_RData/',Rbasefilename,'_sim_',sim,'.RData')
        #save(tree, file=Rfilename)

        nsplit <- max(tree$tree$cptable[,'nsplit'])
        nnode <- sum(grepl('leaf',tree$tree$frame$var))
        ret <-  c(sim=sim,k=-Inf, runtime=runtime, nsplit=nsplit, nnode=nnode)

        idx <- tree$tree$where
        idx_test <- predict_class(tree$tree, data_test)
        EVALS <- eval_tree_pred_inout_jlctree(data,data_test,FLAGS$dist, PARMS$slopes, PARMS$parms, idx, idx_test, pseudo_g, pseudo_g_test,tree,timevar=FALSE)
        ret <- c(ret,EVALS)

        RET[RET_iter,] <- ret
        RET_iter <- RET_iter+1
    }

    cat('sim: ',sim,'\n')

    filename <- paste0(FLAGS$outdir,'/simret_varmaj/',RETbasefilename,'.csv')
    write.table(RET, file=filename, sep=',',col.names=TRUE, quote=FALSE)
}


