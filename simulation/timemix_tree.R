# The main model.
library(optparse)
library(data.table)
library(plyr)
library(jlctree)
source('util.R')
option_list <- list(make_option(c("-N", "--Nsub"), type="numeric", default=500, help=""),
                    make_option(c("-c", "--censor"), type="numeric", default=0.2, help=""),
                    make_option(c("-d", "--dist"), type="character", default='exponential', help=""),
                    make_option(c("-o", "--outdir"), type="character", default='.', help="Output directory."),
                    make_option(c("-l", "--minsim"), type="numeric", default=1, help="min sim"),
                    make_option(c("-u", "--maxsim"), type="numeric", default=200, help="max sim"),
                    make_option(c("-a", "--struct"), type="character", default='tree',
                                help="latent class generating structure, tree (4 classes), XOR (2 classes), linear, or nonlinear."),
                    make_option(c("-g", "--alg"), type="character", default='jlctree', 
                                help="algorithm: jlctree or jlcmm"),
                    make_option(c("-s", "--stop_thre"), type="numeric", default=3.84, help="if set to Inf, no split"),
                    make_option(c("-t", "--test"), type="character", default='lrt', help="Test statistics, rsq, lrt or wald."),
                    make_option(c("-m", "--majprob"), type="numeric", default=0.85, help="Maximum probablity for majority family"),
                    make_option(c("-x", "--continuous"), type="logical", default=TRUE, help="Whether the predictors X1, X2 are continuous"),
                    make_option(c("-e", "--greedy"), type="logical", default=TRUE, help="Whether use greedy jlctree")
                    )


opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

FLAGS <- opt
FLAGS$alg <- 'jlctree'
if (is.null(FLAGS$stop_thre) | is.null(FLAGS$test)){
    stop('Missing argument: stop_thre and test.')
}

if(is.null(FLAGS$greedy)){FLAGS$greedy <- TRUE}
PARMS <- get_parms(FLAGS$dist); parms <- PARMS$parms; slopes <- PARMS$slopes; lam_D <- PARMS$lam_D; 
Nsim <- FLAGS$maxsim-FLAGS$minsim+1

opt2 <- FLAGS; opt2$help <- NULL; opt2$outdir<- NULL; 
RETbasefilename <- paste0(paste0(names(opt2),"_",opt2),collapse="_")


RET <- matrix(0,ncol=78,nrow=Nsim) 
RET_iter <- 1

for (sim in c(FLAGS$minsim:FLAGS$maxsim)){
    set.seed(sim)

    DATA <- gen_data_timemix(FLAGS, PARMS,seed=sim)
    data <- DATA$data; pseudo_g <- DATA$pseudo_g

    DATA_TEST <- gen_data_timemix(FLAGS,PARMS,seed=sim+623)
    data_test <- DATA_TEST$data; pseudo_g_test <- DATA_TEST$pseudo_g

    data <- data.table(data); data_test <- data.table(data_test)
    data[,X1_one:=rep(X1[1],.N),by=ID]; data_test[,X1_one:=rep(X1[1],.N),by=ID]
    data[,X2_one:=rep(X2[1],.N),by=ID]; data_test[,X2_one:=rep(X2[1],.N),by=ID]
    data[,X3_one:=rep(X3[1],.N),by=ID]; data_test[,X3_one:=rep(X3[1],.N),by=ID]
    data[,X4_one:=rep(X4[1],.N),by=ID]; data_test[,X4_one:=rep(X4[1],.N),by=ID]
    data[,X5_one:=rep(X5[1],.N),by=ID]; data_test[,X5_one:=rep(X5[1],.N),by=ID]
    data <- as.data.frame(data); data_test <- as.data.frame(data_test)

    treelist <- list()
    runtime <- rep(0,4)

    tik <- Sys.time()
    treelist[[1]] <- tryCatch({jlctree(survival=Surv(time_L,time_Y,delta)~X3_one+X4_one+X5_one,
        fixed=y~X1+X2+X3+X4+X5, random=~1,
        classmb=~X1_one+X2_one+X3_one+X4_one+X5_one,
        subject='ID',data=data,
        control=rpart.control(cp=-Inf, maxdepth=3),
        parms=list(test.stat=FLAGS$test, 
                   stop.thre=FLAGS$stop_thre,
                   min.nevents=3,
                   maxng=6))}, error=function(e){NULL})
    tok <- Sys.time(); runtime[1] <- round(difftime(tok,tik,units='secs'),4)

    tik <- Sys.time()
    treelist[[2]] <- tryCatch({jlctree(survival=Surv(time_L,time_Y,delta)~X3_one+X4_one+X5_one,
        fixed=y~X1+X2+X3+X4+X5, random=~1,
        classmb=~X1+X2+X3+X4+X5,
        subject='ID',data=data,
        control=rpart.control(cp=-Inf, maxdepth=3),
        parms=list(test.stat=FLAGS$test, 
                   stop.thre=FLAGS$stop_thre,
                   min.nevents=3, maxng=6))}, error=function(e){NULL})
    tok <- Sys.time(); runtime[2] <- round(difftime(tok,tik,units='secs'),4)

    tik <- Sys.time()
    treelist[[3]] <- tryCatch({jlctree(survival=Surv(time_L,time_Y,delta)~X3+X4+X5,
        fixed=y~X1+X2+X3+X4+X5, random=~1,
        classmb=~X1_one+X2_one+X3_one+X4_one+X5_one,
        subject='ID',data=data,
        control=rpart.control(cp=-Inf, maxdepth=3),
        parms=list(test.stat=FLAGS$test, 
                   stop.thre=FLAGS$stop_thre,
                   min.nevents=3, maxng=6))}, error=function(e){NULL})
    tok <- Sys.time(); runtime[3] <- round(difftime(tok,tik,units='secs'),4)


    tik <- Sys.time()
    treelist[[4]] <- tryCatch({jlctree(survival=Surv(time_L,time_Y,delta)~X3+X4+X5,
        fixed=y~X1+X2+X3+X4+X5, random=~1,
        classmb=~X1+X2+X3+X4+X5,
        subject='ID',data=data,
        control=rpart.control(cp=-Inf, maxdepth=3),
        parms=list(test.stat=FLAGS$test, 
                   stop.thre=FLAGS$stop_thre,
                   min.nevents=3, maxng=6))}, error=function(e){NULL})
    tok <- Sys.time(); runtime[4] <- round(difftime(tok,tik,units='secs'),4)

    modelname <- c('var1_var1','var_var1','var1_var','var_var')
    nsplit <- unlist(lapply(treelist,function(x){max(x$tree$cptable[,'nsplit'])}))
    nnode <- unlist(lapply(treelist,function(x){sum(grepl('leaf',x$tree$frame$var))}))

    names(runtime) <- paste0('runtime_',modelname)
    names(nsplit) <- paste0('nsplit_',modelname)
    names(nnode) <- paste0('nnode_',modelname)

    ret <- c(sim=sim, k=-Inf,runtime, nsplit, nnode)
    EVALSname  <- c("ISE1","ISE2","ISE3","MSE_b1","MSE_b2","MSE_b3","MSE_y","purity",
                    "ISE_test1","ISE_test2","ISE_test3","MSE_y_test","purity_test",
                    "CIcover1","CIcover2","CIcover3")

    for(i in c(1:4)){
        if (!is.null(treelist[[i]])){
            if(FLAGS$stop_thre == Inf){
                idx <- rep(1, nrow(data)); idx_test <- rep(1,nrow(data_test))
            } else {
                idx <- treelist[[i]]$tree$where
                idx_test <- predict_class(treelist[[i]]$tree, data_test); 
            }
            EVALS <- eval_tree_pred_inout_jlctree(data,data_test,FLAGS$dist, PARMS$slopes, PARMS$parms, 
                                                  idx, idx_test, pseudo_g, pseudo_g_test,treelist[[i]],timevar=TRUE)

            if(FLAGS$stop_thre==Inf){
                # Because there is only one node, all three cox models are the same. 
                # In particular, diffh_diffs cannot be fit. Thus set its value.
                EVALS['ISE1']  <- EVALS['ISE2']; EVALS['ISE_test1']  <- EVALS['ISE_test2']; 
                EVALS['MSE_b1']  <-  EVALS['MSE_b2']; EVALS['CIcover1'] <- EVALS['CIcover2']
            }
        } else {
            EVALS  <- rep(NA, length(EVALSname))
        }
        names(EVALS) <- paste0(EVALSname,'_',modelname[i])
        ret <- c(ret,EVALS)
    }

    RET[RET_iter,] <- ret
    RET_iter <- RET_iter+1
    cat('sim: ',sim,'\n')

    colnames(RET) <- names(ret)
    filename <- paste0(FLAGS$outdir,'/', RETbasefilename,'_timemix.csv')
    write.table(RET, file=filename, sep=',',col.names=TRUE, quote=FALSE)
}

