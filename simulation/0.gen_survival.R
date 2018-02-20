gen_vary_survival <- function
(n,
 x, times, type,
 dist=c('exponential','weibull','gompertz')[1],
 parms
 ){
    if ((type=='I' & ncol(times) !=1) |
        (type=='II' & ncol(times) !=3)){
        stop("Wrong times matrix.")
    }
    if (type=='I'){
        gen_vary_survival_I(n,x,times,dist,parms)
    } else if (type=='II'){
        gen_vary_survival_II(n,x,times,dist,parms)
    } else {
        stop('Wrong type. Should be I or II.')
    }
}

gen_vary_survival_I <- function
(n,
 x, t0,
 dist=c('exponential','weibull','gompertz')[1], 
 parms
 ){
    if (length(x)!=n | length(t0)!=n){
        stop("Wrong input length.")
    }
    lam <- parms$lambda; nu <- alp <- parms$shape
    b <- parms$beta; bz <- parms$betaz
    logu <- log(runif(n))
    lebx <- lam*exp(b*x)

    if (dist=='exponential'){
        time_T <- 
            ifelse(-logu <  lebx*t0,
                   -logu/lebx,
                   (-logu-lebx*t0+ lebx*exp(bz)*t0)/(lebx*exp(bz)))
    } else if (dist=='weibull'){
        time_T <- 
            ifelse(-logu <  lebx*(t0^nu),
                   (-logu/lebx)^(1/nu),
                   ((-logu-lebx*(t0^nu)+lebx*exp(bz)*(t0^nu))/(lebx*exp(bz)))^(1/nu))
    } else if (dist == 'gompertz'){
        time_T <- 
            ifelse(-logu < lebx*(exp(alp*t0)-1)/alp,
                   log(1-alp*logu/lebx)/alp,
                   log(-alp*logu/(lebx*exp(bz))-(exp(alp*t0)-1-exp(bz+alp*t0))/exp(bz))/alp)
    }
    return(time_T)
}

gen_vary_survival_II <- function
(n,
 x, times,
 dist=c('exponential','weibull','gompertz')[1], 
 parms
 ){
    if (length(x)!=n | nrow(times)!=n){
        stop("Wrong input length.")
    }
    lam <- parms$lambda; nu <- alp <- parms$shape
    b <- parms$beta; bz <- parms$betaz
    logu <- log(runif(n))
    lebx <- lam*exp(b*x)
    t1 <- times[,1]; t2 <- times[,2]; t3 <- times[,3];

    if (dist=='exponential'){
        R1 <- lebx*t1
        T1 <- -logu/lebx
        R2 <- R1 + lebx*exp(bz)*(t2-t1)
        T2 <- (-logu -R1+lebx*exp(bz)*t1)/(lebx*exp(bz))
        R3 <- R2 + lebx*(t3-t2)
        T3 <- (-logu -R2+lebx*t2)/lebx
        T4 <- (-logu-R3+lebx*exp(bz)*t3)/(lebx*exp(bz))

    } else if (dist=='weibull'){
        R1 <- lebx*(t1^nu)
        T1 <- (-logu/lebx)^(1/nu)
        R2 <- R1 + lebx*exp(bz)*(t2^nu-t1^nu)
        T2 <- ((-logu -R1+lebx*exp(bz)*(t1^nu))/(lebx*exp(bz)))^(1/nu)
        R3 <- R2 + lebx*(t3^nu-t2^nu)
        T3 <- ((-logu -R2+lebx*t2^nu)/lebx)^(1/nu)
        T4 <- ((-logu-R3+lebx*exp(bz)*t3^nu)/(lebx*exp(bz)))^(1/nu)

    } else if (dist == 'gompertz'){
        R1 <- lebx*(exp(alp*t1)-1)/alp
        T1 <- log(1-alp*logu/lebx)/alp
        R2 <- R1 + lebx*(exp(bz+alp*t2)-exp(bz+alp*t1))/alp
        T2 <- log(-alp*logu/(lebx*exp(bz))-(exp(alp*t1)-1-exp(bz+alp*t1))/(exp(bz)))/alp
        R3 <- R2 + lebx*(exp(alp*t3)-exp(alp*t2))/alp
        T3 <- log(1-alp*logu/lebx-exp(alp*t1)-exp(bz+alp*t2)+exp(bz+alp*t1)+exp(alp*t2))/alp
        T4 <- log(-alp*logu/(lebx*exp(bz))-
                  (exp(alp*t1)-1+exp(bz+alp*t2)-exp(bz+alp*t1)+exp(alp*t3)-exp(alp*t2)-exp(bz+alp*t3))/(exp(bz)))/alp

    }
    time_T <- ifelse(-logu <  R1, T1,
                     ifelse(-logu < R2, T2,
                            ifelse(-logu < R3, T3, T4)))

    return(time_T)
}


gen_model3_survival <- function
(ebx, #exp(beta^T x) 
 dist=c('exponential','weibulld','weibulli')[1],#,'bathtub')[1],
 parms
 ){
    n <- length(ebx)
    logu <- log(runif(n))

    if (dist=='exponential'){
        time_T <- -logu / (ebx*parms$lambda)
    } else if (dist=='weibulld' | dist=='weibulli'){
        time_T <- parms$beta * (-logu / ebx)^(1/parms$alp)
    } 
}
