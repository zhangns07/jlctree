# The main model.
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

PARMS <- get_parms(FLAGS$dist)
parms <- PARMS$parms; 
Nsim <- FLAGS$maxsim-FLAGS$minsim+1

opt2 <- FLAGS; opt2$help <- NULL; opt2$outdir<- NULL; 
RETbasefilename <- paste0(paste0(names(opt2),"_",opt2),collapse="_")

opt3 <- FLAGS; opt3$help <- NULL; opt3$outdir<- NULL; opt3$minsim<-NULL; opt3$maxsim<-NULL
Rbasefilename <- paste0(paste0(names(opt3),"_",opt3),collapse="_")

if (FLAGS$alg == 'jlcmm'){
    RET <- matrix(0,ncol=17,nrow=Nsim)
    colnames(RET) <- c('sim','runtime','bestng','B2','B3','B4','B5','B6',
                       'ISE','MSE_b','MSE_y','purity',
                       'ISE_test_max','MSE_y_test_max','ISE_test_avg','MSE_y_test_avg','purity_test')
} else if (FLAGS$alg == 'jlctree'){
    if (FLAGS$stop_thre==-1){ RET<- matrix(0,ncol=11,nrow=Nsim)
    } else { RET <- matrix(0,ncol=12,nrow=Nsim*5) }
    colnames(RET) <- c('sim','k','runtime','nsplit','nnode',
                       'ISE','MSE_b','MSE_y','purity',
                       'ISE_test','MSE_y_test','purity_test')
}
RET_iter <- 1

filename <- paste0(FLAGS$outdir,'/simret_main/',RETbasefilename,'.csv')
INFO <- read.table(filename, sep=',')

for (sim in c(FLAGS$minsim:FLAGS$maxsim)){
    DATA <- gen_data(FLAGS, PARMS,seed=sim)
    data <- DATA$data; pseudo_g <- DATA$pseudo_g

    DATA_TEST <- gen_data(FLAGS,PARMS,seed=sim+623)
    data_test <- DATA_TEST$data; pseudo_g_test <- DATA_TEST$pseudo_g

    if (FLAGS$alg == 'jlcmm'){
        currINFO <- INFO[RET_iter,1:8]
        best_ng <- currINFO['bestng']

        Rfilename <- paste0(FLAGS$outdir,'/simret_main_RData/',Rbasefilename,'_ng_',best_ng,'_sim_',sim,'.RData')
        load(Rfilename)
        mod <- get(paste0('m',best_ng))

        EVALS <- eval_lcmm_pred_inout(data,data_test, FLAGS$dist, PARMS$slopes,
                                      PARMS$parms, mod, pseudo_g, pseudo_g_test)
        RET[RET_iter,] <- c(currINFO, EVALS)
        RET_iter <- RET_iter+1

    } else if (FLAGS$alg == 'jlctree'){
        if (FLAGS$stop_thre==-1){
            idx <- rep(1,nrow(data))
            idx_test <- rep(1,nrow(data_test))
            EVALS <- eval_tree_pred_inout(data,data_test,FLAGS$dist, PARMS$slopes, PARMS$parms,
                                          idx, idx_test, pseudo_g, pseudo_g_test)
            RET[RET_iter,] <- c(sim,k=-Inf, runtime=0, nsplit=-1, nnode=1,EVALS)
            RET_iter <- RET_iter+1

        } else {
            currINFO <- INFO[RET_iter, c(1:5)]
            Rfilename <- paste0(FLAGS$outdir,'/simret_main_RData/',Rbasefilename,'_sim_',sim,'.RData')
            load(Rfilename)

            survs <- survs_v3
            survlist <- list(eval=surve, split=survs, init=survi)

            idx <- cond_ind_tree$where
            idx_test <- predict_class(cond_ind_tree, data_test)
            EVALS <- eval_tree_pred_inout(data,data_test,FLAGS$dist, PARMS$slopes, PARMS$parms,
                                          idx, idx_test, pseudo_g, pseudo_g_test)
            RET[RET_iter,] <- c(currINFO, EVALS)
            RET_iter <- RET_iter+1

            cventry <- which.min(cptable[, "xerror"])
            xerrorcv <- cptable[cventry, "xerror"]

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
                    idx <- cond_ind_tree_prune$where
                    idx_test <- predict_class(cond_ind_tree_prune, data_test)
                    EVALS <- eval_tree_pred_inout(data,data_test,FLAGS$dist, PARMS$slopes, PARMS$parms,
                                                  idx, idx_test, pseudo_g, pseudo_g_test)
                    RET[RET_iter,] <- c(sim,k=kse, currINFO['runtime'], nsplit_prune, nnode_prune, EVALS)
                }
                RET_iter <- RET_iter+1

            }
        }
    }

    cat('sim: ',sim,'\n')


    filename <- paste0(FLAGS$outdir,'/simret_main_eval/',RETbasefilename,'.csv')
    write.table(RET, file=filename, sep=',',col.names=TRUE, quote=FALSE)
}




