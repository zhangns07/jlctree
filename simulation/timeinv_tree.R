# The main model.
library(optparse)
library(data.table)
library(plyr)
library(jlctree)
source('0.init.R')
source('../util.R')
option_list <- list(make_option(c("-N", "--Nsub"), type="numeric", default=500, help=""),
                    make_option(c("-c", "--censor"), type="numeric", default=0.2, help=""),
                    make_option(c("-d", "--dist"), type="character", default='weibulli', help=""),
                    make_option(c("-o", "--outdir"), type="character", default='.', help="Output directory."),
                    make_option(c("-l", "--minsim"), type="numeric", default=1, help="min sim"),
                    make_option(c("-u", "--maxsim"), type="numeric", default=200, help="max sim"),
                    make_option(c("-a", "--struct"), type="character", default='tree',
                                help="latent class generating structure, tree (4 classes), XOR (2 classes), linear, or nonlinear."),
                    make_option(c("-s", "--stop_thre"), type="numeric", default=3.84, help="when set to Inf, do not split"),
                    make_option(c("-t", "--test"), type="character", default='lrt', help="Test statistics, rsq, lrt or wald."),
                    make_option(c("-m", "--majprob"), type="numeric", default=0.85, help="Maximum probablity for majority family"),
                    make_option(c("-x", "--continuous"), type="logical", default=TRUE, help="Whether the predictors X1, X2 are continuous"),
                    make_option(c("-e", "--greedy"), type="logical", default=NULL, help="Whether use greedy jlctree")
                    )


opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

FLAGS <- opt
FLAGS$alg <- 'jlctree' # needed when generating data
if (is.null(FLAGS$stop_thre) | is.null(FLAGS$test)){
    stop('Missing argument: stop_thre and test.')
}
if(is.null(FLAGS$greedy)){FLAGS$greedy<- TRUE}

PARMS <- get_parms(FLAGS$dist); parms <- PARMS$parms; slopes <- PARMS$slopes; lam_D <- PARMS$lam_D; 
Nsim <- FLAGS$maxsim-FLAGS$minsim+1

opt2 <- FLAGS; opt2$help <- NULL; opt2$outdir<- NULL; 
RETbasefilename <- paste0(paste0(names(opt2),"_",opt2),collapse="_")


if(is.null(FLAGS$continuous)){FLAGS$continuous <- FALSE}
RET <- matrix(0,ncol=23,nrow=Nsim) 
colnames(RET) <- c('sim','k','runtime','nsplit','nnode',
                   "ISE1","ISE2","ISE3", 
                   "MSE_b1","MSE_b2","MSE_b3", 
                   "MSE_y","purity",
                   "ISE_test1","ISE_test2","ISE_test3",
                   "MSE_y_test","purity_test",
                   "CIcover1", "CIcover2", "CIcover3","MSE_beta","CIcover_beta")

RET_iter <- 1
for (sim in c(FLAGS$minsim:FLAGS$maxsim)){
    set.seed(sim)

    DATA <- gen_data(FLAGS, PARMS,seed=sim)
    data <- DATA$data; pseudo_g <- DATA$pseudo_g

    DATA_TEST <- gen_data(FLAGS,PARMS,seed=sim+623)
    data_test <- DATA_TEST$data; pseudo_g_test <- DATA_TEST$pseudo_g


    if (FLAGS$stop_thre==-1){FLAGS$stop_thre = Inf}
    tik <- Sys.time()
    tree <- jlctree(survival=Surv(time_L,time_Y,delta)~X3+X4+X5,
                    fixed=y~X1+X2+X3+X4+X5,
                    random=~1,
                    classmb=~X1+X2+X3+X4+X5,
                    subject='ID',data=data,
                    control=list(cp=-Inf,maxdepth=3),
                    parms=list(teststat=FLAGS$test, 
                               stop.thre=FLAGS$stop_thre, min.nevents=3,maxng=6))
    tok <- Sys.time()
    runtime <- round(difftime(tok,tik,units='secs'),4)

    if(FLAGS$stop_thre==Inf){
        idx <- rep(1, nrow(data)); idx_test <- rep(1,nrow(data_test))
        ret <- c(sim,k=-Inf, runtime=runtime, nsplit=-1, nnode=1)
        data$node <- 1
        lmmmodel <- lmer(y ~  (1|ID),data=data)
        tree$lmmmodel <- lmmmodel
    } else {
        # fit a lmm model with node-specific intercept
        data$node <- tree$tree$where
        lmmmodel <- lmer(y ~ 0+factor(node) + (1|ID),data=data)
        tree$lmmmodel <- lmmmodel
        nsplit <- max(tree$tree$cptable[,'nsplit'])
        nnode <- sum(grepl('leaf',tree$tree$frame$var))
        ret <-  c(sim=sim,k=-Inf, runtime=runtime, nsplit=nsplit, nnode=nnode)
        idx <- tree$tree$where
        idx_test <- predict_class(tree$tree, data_test)
    }

    EVALS <- eval_tree_pred_inout_jlctree(data,data_test,FLAGS$dist, PARMS$slopes, PARMS$parms, 
                                          idx, idx_test, pseudo_g, pseudo_g_test,tree,timevar=FALSE)
    if(FLAGS$stop_thre==Inf){
        # Because there is only one node, all three cox models are the same. 
        # In particular, diffh_diffs cannot be fit. Thus set its value.
        EVALS['ISE1']  <- EVALS['ISE2']; EVALS['ISE_test1']  <- EVALS['ISE_test2']; 
        EVALS['MSE_b1']  <-  EVALS['MSE_b2']; EVALS['CIcover1'] <- EVALS['CIcover2']
    }

    ret <- c(ret,EVALS)

    RET[RET_iter,] <- ret
    RET_iter <- RET_iter+1

    cat('sim: ',sim,'\n')
    filename <- paste0(FLAGS$outdir, '/', RETbasefilename,'_timeinv.csv')
    write.table(RET, file=filename, sep=',',col.names=TRUE, quote=FALSE)
}



