# The main model.
library(optparse)
library(data.table)
library(plyr)
library(lcmm)
source('0.init.R')
option_list <- list(make_option(c("-N", "--Nsub"), type="numeric", default=100, help=""),
                    make_option(c("-c", "--censor"), type="numeric", default=0, help=""),
                    make_option(c("-d", "--dist"), type="character", default='exponential', help=""),
                    make_option(c("-o", "--outdir"), type="character", default='.', help="Output directory."),
                    make_option(c("-l", "--minsim"), type="numeric", default=1, help="min sim"),
                    make_option(c("-u", "--maxsim"), type="numeric", default=200, help="max sim"),
                    make_option(c("-a", "--struct"), type="character", default='tree',
                                help="latent class generating structure, tree (4 classes), XOR (2 classes), linear, or nonlinear."),
                    make_option(c("-x", "--continuous"), type="logical", default=TRUE, help="Whether the predictors X1, X2 are continuous"),
                    make_option(c("-m", "--majprob"), type="numeric", default=0.85, help="Maximum probablity for majority family")
                    )


opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

FLAGS <- opt; FLAGS$alg <- 'jlcmm'
PARMS <- get_parms(FLAGS$dist); parms <- PARMS$parms; slopes <- PARMS$slopes; lam_D <- PARMS$lam_D; 
Nsim <- FLAGS$maxsim-FLAGS$minsim+1

opt2 <- FLAGS; opt2$help <- NULL; opt2$outdir<- NULL; 
RETbasefilename <- paste0(paste0(names(opt2),"_",opt2),collapse="_")


RET <- matrix(0,ncol=17,nrow=Nsim)
colnames(RET) <- c('sim','runtime','bestng','B2','B3','B4','B5','B6',
                   'ISE','MSE_b','MSE_y','purity',
                   'ISE_test_max','MSE_y_test_max','ISE_test_avg','MSE_y_test_avg','purity_test')

RET_iter <- 1

for (sim in c(FLAGS$minsim:FLAGS$maxsim)){
    set.seed(sim)

    DATA <- gen_data(FLAGS, PARMS,seed=sim)
    data <- DATA$data; pseudo_g <- DATA$pseudo_g

    DATA_TEST <- gen_data(FLAGS,PARMS,seed=sim+623)
    data_test <- DATA_TEST$data; pseudo_g_test <- DATA_TEST$pseudo_g

    m1 <- Jointlcmm(fixed=y~X1+X2+X3+X4+X5,
                    subject='ID',
                    survival = Surv(time_L,time_Y,delta)~X3+X4+X5,
                    hazard="Weibull",hazardtype="Specific",ng=1,data=data)

    classmb <- ~X1*X2+X3+X4+X5

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
    mod <- get(paste0('m',best_ng))
    runtime <- round(difftime(tok,tik,units='secs'),4)
    EVALS <- eval_lcmm_pred_inout(data,data_test, FLAGS$dist, PARMS$slopes,
                                  PARMS$parms, mod, pseudo_g, pseudo_g_test)
    ret <- c(sim, runtime, best_ng, initB[2:6], EVALS)
    RET[RET_iter,] <- ret 


    RET_iter <- RET_iter+1
    cat('sim: ',sim,'\n')

    filename <- paste0(FLAGS$outdir,RETbasefilename,'_timeinv.csv')
    write.table(RET, file=filename, sep=',',col.names=TRUE, quote=FALSE)
}

