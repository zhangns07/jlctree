# model4: joint lcmm model
# with y being longitudinal
# with ng chosing by BIC
library(optparse)
library(data.table)
library(plyr)
library(lcmm)
source('0.init.R')
source('0.gen_survival.R')
option_list <- list(
                    make_option(c("-N", "--Nsub"), type="numeric", default=100, help=""),
                    make_option(c("-c", "--censor"), type="numeric", default=0, help=""),
                    make_option(c("-d", "--dist"), type="character", default='exponential', help=""),
                    make_option(c("-o", "--outdir"), type="character", default='.', help="Output directory."),
                    make_option(c("-l", "--minsim"), type="numeric", default=1, help="min sim"),
                    make_option(c("-u", "--maxsim"), type="numeric", default=200, help="max sim"))


opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

FLAGS <- opt
Nsub <- FLAGS$Nsub
censor_rate <- FLAGS$censor
dist <- FLAGS$dist


sd_ranef <- 0.2
sd_e <- 0.1
Nsim <- FLAGS$maxsim-FLAGS$minsim+1
beta_y <- rep(0,5)
parms <- switch(dist,
                exponential=list(lambda=0.1),
                weibulld=list(alp=0.9, beta=1),
                weibulli=list(alp=3, beta=2))


# By average ratios
slopes <- switch(dist,
                 exponential=matrix(c(0,0,0,0.56,0.56,0.09,0.92,0.92,0.15,1.46,1.46,0.24),ncol=3,byrow=TRUE),
                 weibulld=matrix(c(-1.17,-1.17,-0.19,-0.66,-0.66,-0.11,-0.55,-0.55,-0.09,0,0,0),ncol=3,byrow=TRUE),
                 weibulli=matrix(c(-3.22,-3.22,-0.54,-2.26,-2.26,-0.38, -1.53,-1.53,-0.26, 0,0,0),ncol=3,byrow=TRUE))

lam_D <- switch(dist,
                exponential=list('0.2'=c(0.026,0.048,0.086,0.156),
                                 '0.5'=c(0.114,0.232,0.396,0.74)),
                weibulld=list('0.2'=c(0.024,0.072,0.084,0.238),
                              '0.5'=c(0.13,0.318,0.384,0.956)),
                weibulli=list('0.2'=c(0.024,0.04,0.076,0.19),
                              '0.5'=c(0.084,0.168,0.224,0.64)))

opt2 <- FLAGS; opt2$help <- NULL; opt2$outdir<- NULL; opt2$model <- 'model4'; opt2$alg<-'lcmm'
RETbasefilename <- paste0(paste0(names(opt2),"_",opt2),collapse="_")

opt3 <- FLAGS; opt3$help <- NULL; opt3$outdir<- NULL; opt3$model <- 'model4'; opt3$alg<-'lcmm'; opt3$minsim<-NULL; opt3$maxsim<-NULL
Rbasefilename <- paste0(paste0(names(opt3),"_",opt3),collapse="_")


RET <- matrix(0,ncol=12,nrow=Nsim)
colnames(RET) <- c('sim','time','bestng','B2','B3','B4','B5','B6','ISE','MSE_b','MSE_y','purity')
RET_iter <- 1

W <- matrix(c( -1,-1, -1,1,1,-1,1,1),ncol=2,byrow=TRUE)

for (sim in c(FLAGS$minsim:FLAGS$maxsim)){
    set.seed(sim)

    # X1 - X5
    X1 <- as.numeric(runif(2*Nsub)>0.5)
    X2 <- as.numeric(runif(2*Nsub)>0.5)
    X3 <- as.numeric(runif(2*Nsub)>0.5)
    X4 <- round(runif(2*Nsub),1)
    X5 <- sample(c(1:5),2*Nsub,replace=TRUE)
    X <- cbind(X1,X2,X3,X4,X5)

    trueg <- 2*(X1)+X2+1
    g <- apply(cbind(X1,X2),1,function(x){
                   weights <- exp(2*W %*% (x-0.5))
                   weights <- weights/sum(weights)
                   cumweights <- cumsum(weights)
                   which(runif(1)  < cumweights)[1] })

    ebx <- exp(rowSums(slopes[g,] * cbind(X3,X4,X5)))
    time_T <- gen_model3_survival(ebx, dist, parms)
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
                           ret <- cbind(ID=i, X1=rep(X[i,1],num_i),
                                        X2=rep(X[i,2],num_i), X3=rep(X[i,3],num_i),
                                        X4=rep(X[i,4],num_i), X5=rep(X[i,5],num_i),
                                        time_L=rep(time_L[i],num_i),time_Y=rep(time_Y[i],num_i),
                                        delta=rep(delta[i],num_i)) })


    ranef <- rnorm(Nsub, sd=sd_ranef)
    ranefs <- ranef[LTRC_data$ID]
    fixef <- c(0,1,1,2)
    pseudo_g <- g[LTRC_data$ID]
    y <- fixef[pseudo_g] + ranefs + rnorm(nrow(LTRC_data), sd=sd_e) 
    data <- cbind(LTRC_data,y)

    m1 <- Jointlcmm(fixed=y~X1+X2+X3+X4+X5,
                    subject='ID',
                    survival = Surv(time_L,time_Y,delta)~X3+X4+X5,
                    hazard="Weibull",hazardtype="Specific",ng=1,data=data)

    BICs <- rep(Inf, 6)
    initB <- rep(1,6)#1:B=m1; 2:B=NULL; 3:B=random(m1); 4:none worked
    tik <- Sys.time()
    for(ng in c(2:6)){
        mod <- Jointlcmm(fixed=y~X1+X2+X3+X4+X5,mixture=~1,
                        classmb=~X1+X2+X3+X4+X5,subject='ID',
                        survival = Surv(time_L,time_Y,delta)~mixture(X3+X4+X5),
                        hazard="Weibull",hazardtype="Specific",ng=ng,data=data,B=m1)

        if (mod$conv %in% c(4,5)){
            mod <- Jointlcmm(fixed=y~X1+X2+X3+X4+X5,mixture=~1,
                            classmb=~X1+X2+X3+X4+X5,subject='ID',
                            survival = Surv(time_L,time_Y,delta)~mixture(X3+X4+X5),
                            hazard="Weibull",hazardtype="Specific",ng=ng,data=data)#,maxiter=30)
            initB[ng] <- 2
        }

        if (mod$conv %in% c(4,5)){
            mod <- Jointlcmm(fixed=y~X1+X2+X3+X4+X5,mixture=~1,
                            classmb=~X1+X2+X3+X4+X5,subject='ID',
                            survival = Surv(time_L,time_Y,delta)~mixture(X3+X4+X5),
                            hazard="Weibull",hazardtype="Specific",ng=ng,data=data,B=random(m1))#,maxiter=30)
            initB[ng] <- 3
        }

        if (!mod$conv %in% c(4,5)){
            BICs[ng] <- mod$BIC
        } else {
            initB[ng] <- 4
        }

        assign(paste0('m',ng), mod)
        Rfilename <- paste0(FLAGS$outdir,'/simret_model4_RData/',Rbasefilename,'_ng_',ng,'_sim_',sim,'.RData')
        save(list=paste0('m',ng), file=Rfilename)
    }
    tok <- Sys.time()

    best_ng <- which.min(BICs)
    mod <- get(paste0('m',best_ng))
    runtime <- round(difftime(tok,tik,units='secs'),4)
    RET[RET_iter,] <- c(sim, runtime, best_ng, initB[2:6], eval_lcmm_pred(data,dist,slopes,parms,mod,g=pseudo_g))
    RET_iter <- RET_iter+1

    cat('sim: ',sim,'\n')

    filename <- paste0(FLAGS$outdir,'/simret_model4/',RETbasefilename,'.csv')
    write.table(RET, file=filename, sep=',',col.names=TRUE, quote=FALSE)
}
