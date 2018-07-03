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
                    make_option(c("-g", "--alg"), type="character", default='jlctree', 
                                help="algorithm: jlctree or jlcmm"),
                    make_option(c("-s", "--stop_thre"), type="numeric", default=NULL, help=""),
                    make_option(c("-t", "--test"), type="character", default=NULL, help="Test statistics, rsq, lrt or wald."),
                    make_option(c("-x", "--continuous"), type="logical", default=TRUE, help="Whether the predictors X1, X2 are continuous"),
                    make_option(c("-m", "--majprob"), type="numeric", default=0.85, help="Maximum probablity for majority family"),
                    make_option(c("-w", "--survvar"), type="logical", default=TRUE, help="Whether X3-X5 is time varying"),
                    make_option(c("-e", "--greedy"), type="logical", default=NULL, help="Whether use greedy jlctree")
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

if(is.null(FLAGS$greedy)){FLAGS$greedy <- TRUE}
PARMS <- get_parms(FLAGS$dist); parms <- PARMS$parms; slopes <- PARMS$slopes; lam_D <- PARMS$lam_D; 
Nsim <- FLAGS$maxsim-FLAGS$minsim+1

opt2 <- FLAGS; opt2$help <- NULL; opt2$outdir<- NULL; opt2$survvar <- NULL
RETbasefilename <- paste0(paste0(names(opt2),"_",opt2),collapse="_")

opt3 <- FLAGS; opt3$help <- NULL; opt3$outdir<- NULL; opt3$minsim<-NULL; opt3$maxsim<-NULL; opt3$survvar <- NULL
Rbasefilename <- paste0(paste0(names(opt3),"_",opt3),collapse="_")

if(is.null(FLAGS$continuous)){FLAGS$continuous <- FALSE}
if(is.null(FLAGS$survvar)){FLAGS$survvar <- FALSE}

if (FLAGS$survvar){
    RDatadir <- 'simret_survvar_RData/' 
    RETdir <- 'simret_survvar/'
} else {
    RDatadir <- 'simret_timevar_RData/' 
    RETdir <- 'simret_timevar/'
} 

RET <- matrix(0,ncol=66,nrow=Nsim) 
filename <- paste0(FLAGS$outdir,'/',RETdir, RETbasefilename,'.csv')
if(file.exists(filename)){
    RET <- read.table(filename,sep=",",header=TRUE)
    RET <- as.matrix(RET)
    lastsim <- which(RET[,2] ==0)[1]
    if(is.na(lastsim)){
        stop('All simululations are done.')
    }
} else {
    lastsim <- FLAGS$minsim
}


RET_iter <- lastsim
for (sim in c(lastsim:FLAGS$maxsim)){
    set.seed(sim)

    DATA <- gen_data_timevar(FLAGS, PARMS,seed=sim, FLAGS$survvar)
    data <- DATA$data; pseudo_g <- DATA$pseudo_g

    DATA_TEST <- gen_data_timevar(FLAGS,PARMS,seed=sim+623, FLAGS$survvar)
    data_test <- DATA_TEST$data; pseudo_g_test <- DATA_TEST$pseudo_g

    data <- data.table(data); data_test <- data.table(data_test)
    data[,X1_one:=rep(X1[1],.N),by=ID]; data_test[,X1_one:=rep(X1[1],.N),by=ID]
    data[,X2_one:=rep(X2[1],.N),by=ID]; data_test[,X2_one:=rep(X2[1],.N),by=ID]
    data[,X3_one:=rep(X3[1],.N),by=ID]; data_test[,X3_one:=rep(X3[1],.N),by=ID]
    data[,X4_one:=rep(X4[1],.N),by=ID]; data_test[,X4_one:=rep(X4[1],.N),by=ID]
    data[,X5_one:=rep(X5[1],.N),by=ID]; data_test[,X5_one:=rep(X5[1],.N),by=ID]
    data <- as.data.frame(data); data_test <- as.data.frame(data_test)


    if (FLAGS$stop_thre==-1){
        treelist <- list()
        runtime <- rep(0,4)
        tik <- Sys.time()
        treelist[[1]] <- jlctree(survival=Surv(time_L,time_Y,delta)~X3_one+X4_one+X5_one,
                         fixed=y~X1+X2+X3+X4+X5, random=~1,
                         classmb=~X1_one+X2_one+X3_one+X4_one+X5_one,
                         subject='ID',maxng=6,data=data,greedy=greedy,
                         parms=list(test_stat=FLAGS$test, 
                                    stop_thre=Inf,min_nevent=3))
        tok <- Sys.time(); runtime[1] <- round(difftime(tok,tik,units='secs'),4)

        tik <- Sys.time()
        treelist[[2]] <- jlctree(survival=Surv(time_L,time_Y,delta)~X3_one+X4_one+X5_one,
                fixed=y~X1+X2+X3+X4+X5, random=~1,
                classmb=~X1+X2+X3+X4+X5,
                subject='ID',maxng=6,data=data,greedy=greedy,
                parms=list(test_stat=FLAGS$test, 
                           stop_thre=Inf,min_nevent=3))
        tok <- Sys.time(); runtime[2] <- round(difftime(tok,tik,units='secs'),4)

        tik <- Sys.time()
        treelist[[3]] <- jlctree(survival=Surv(time_L,time_Y,delta)~X3+X4+X5,
                fixed=y~X1+X2+X3+X4+X5, random=~1,
                classmb=~X1_one+X2_one+X3_one+X4_one+X5_one,
                subject='ID',maxng=6,data=data,greedy=greedy,
                parms=list(test_stat=FLAGS$test, 
                           stop_thre=Inf,min_nevent=3))
        tok <- Sys.time(); runtime[3] <- round(difftime(tok,tik,units='secs'),4)


        tik <- Sys.time()
        treelist[[4]] <- jlctree(survival=Surv(time_L,time_Y,delta)~X3+X4+X5,
                fixed=y~X1+X2+X3+X4+X5, random=~1,
                classmb=~X1+X2+X3+X4+X5,
                subject='ID',maxng=6,data=data,greedy=greedy,
                parms=list(test_stat=FLAGS$test, 
                           stop_thre=Inf,min_nevent=3))
        tok <- Sys.time(); runtime[4] <- round(difftime(tok,tik,units='secs'),4)

        modelname <- c('var1_var1','var_var1','var1_var','var_var')
        nsplit <- unlist(lapply(treelist,function(x){max(x$tree$cptable[,'nsplit'])}))
        nnode <- unlist(lapply(treelist,function(x){sum(grepl('leaf',x$tree$frame$var))}))

        names(runtime) <- paste0('runtime_',modelname)
        names(nsplit) <- paste0('nsplit_',modelname)
        names(nnode) <- paste0('nnode_',modelname)
        ret <- c(sim=sim, k=-Inf,runtime, nsplit, nnode)

        idx <- rep(1, nrow(data)); idx_test <- rep(1,nrow(data_test))
        for(i in c(1:4)){
            EVALS <- eval_tree_pred_inout_jlctree(data,data_test,FLAGS$dist, PARMS$slopes, PARMS$parms, 
                                                  idx, idx_test, pseudo_g, pseudo_g_test,treelist[[i]],timevar=TRUE)
            names(EVALS) <- paste0(names(EVALS),'_',modelname[i])
            ret <- c(ret,EVALS)
        }
        RET[RET_iter,] <- ret
        RET_iter <- RET_iter+1

    } else {
        treelist <- list()
        runtime <- rep(0,4)

        tik <- Sys.time()
        treelist[[1]] <- jlctree(survival=Surv(time_L,time_Y,delta)~X3_one+X4_one+X5_one,
                         fixed=y~X1+X2+X3+X4+X5, random=~1,
                         classmb=~X1_one+X2_one+X3_one+X4_one+X5_one,
                         subject='ID',maxng=6,data=data,greedy=greedy,
                         parms=list(test_stat=FLAGS$test, 
                                    stop_thre=FLAGS$stop_thre,min_nevent=3))
        tok <- Sys.time(); runtime[1] <- round(difftime(tok,tik,units='secs'),4)

        tik <- Sys.time()
        treelist[[2]] <- jlctree(survival=Surv(time_L,time_Y,delta)~X3_one+X4_one+X5_one,
                fixed=y~X1+X2+X3+X4+X5, random=~1,
                classmb=~X1+X2+X3+X4+X5,
                subject='ID',maxng=6,data=data,greedy=greedy,
                parms=list(test_stat=FLAGS$test, 
                           stop_thre=FLAGS$stop_thre,min_nevent=3))
        tok <- Sys.time(); runtime[2] <- round(difftime(tok,tik,units='secs'),4)

        tik <- Sys.time()
        treelist[[3]] <- jlctree(survival=Surv(time_L,time_Y,delta)~X3+X4+X5,
                fixed=y~X1+X2+X3+X4+X5, random=~1,
                classmb=~X1_one+X2_one+X3_one+X4_one+X5_one,
                subject='ID',maxng=6,data=data,greedy=greedy,
                parms=list(test_stat=FLAGS$test, 
                           stop_thre=FLAGS$stop_thre,min_nevent=3))
        tok <- Sys.time(); runtime[3] <- round(difftime(tok,tik,units='secs'),4)


        tik <- Sys.time()
        treelist[[4]] <- jlctree(survival=Surv(time_L,time_Y,delta)~X3+X4+X5,
                fixed=y~X1+X2+X3+X4+X5, random=~1,
                classmb=~X1+X2+X3+X4+X5,
                subject='ID',maxng=6,data=data,greedy=greedy,
                parms=list(test_stat=FLAGS$test, 
                           stop_thre=FLAGS$stop_thre,min_nevent=3))
        tok <- Sys.time(); runtime[4] <- round(difftime(tok,tik,units='secs'),4)

        modelname <- c('var1_var1','var_var1','var1_var','var_var')
        nsplit <- unlist(lapply(treelist,function(x){max(x$tree$cptable[,'nsplit'])}))
        nnode <- unlist(lapply(treelist,function(x){sum(grepl('leaf',x$tree$frame$var))}))

        names(runtime) <- paste0('runtime_',modelname)
        names(nsplit) <- paste0('nsplit_',modelname)
        names(nnode) <- paste0('nnode_',modelname)

        ret <- c(sim=sim, k=-Inf,runtime, nsplit, nnode)
        for(i in c(1:4)){
            idx <- treelist[[i]]$tree$where
            idx_test <- predict_class(treelist[[i]]$tree, data_test); 
            EVALS <- eval_tree_pred_inout_jlctree(data,data_test,FLAGS$dist, PARMS$slopes, PARMS$parms, 
                                                  idx, idx_test, pseudo_g, pseudo_g_test,treelist[[i]],timevar=TRUE)
            names(EVALS) <- paste0(names(EVALS),'_',modelname[i])
            ret <- c(ret,EVALS)
        }

        RET[RET_iter,] <- ret
        RET_iter <- RET_iter+1

        #Rfilename <- paste0(FLAGS$outdir,'/',RDatadir,Rbasefilename,'_sim_',sim,'.RData')
        #save(cond_ind_tree, cptable, file=Rfilename)
    }

    cat('sim: ',sim,'\n')

    colnames(RET) <- names(ret)
    filename <- paste0(FLAGS$outdir,'/',RETdir, RETbasefilename,'.csv')
    write.table(RET, file=filename, sep=',',col.names=TRUE, quote=FALSE)
}


