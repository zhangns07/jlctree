library(optparse)
library(data.table)
library(plyr)
library(lcmm)
library(ggplot2)
source('0.init.R')
source('0.gen_survival.R')

X1 <- runif(20000,min=0,max=1)
X2 <- runif(20000,min=0,max=1)
gg <- get_latent_class(X1,X2,'tree', 'multinomial', seed=1)


data <- cbind(X1,X2,gg)
g <- ggplot()
g <- g + geom_point(aes(x=X1,y=X2,color=factor(gg)))
print(g)

#----------
# Tree
pdf('struct_tree.pdf')
plot(x=NULL,y=NULL,xlim=c(0,1),ylim=c(0,1),xlab="X1",ylab="X2",main=NULL,xaxs="i", yaxs="i")
abline(v=0.5); abline(h=0.5)
text(x=0.25,y=0.75,'g=2',cex=2)
text(x=0.25,y=0.25,'g=1',cex=2)
text(x=0.75,y=0.25,'g=3',cex=2)
text(x=0.75,y=0.75,'g=4',cex=2)
dev.off()

#----------
# linear
#W <- matrix(c(0.8,-0.6,0.9,0.5,-0.8,0.6,0.5,0.9),byrow=TRUE,ncol=2) # random sample from sphere
#pdf('struct_linear.pdf')
#plot(x=NULL,y=NULL,xlim=c(1,3),ylim=c(1,3),xlab="X1",ylab="X2",main=NULL,xaxs="i", yaxs="i")
#y <- c(0,2); x <- 1.2/1.6*(y-2)+2; lines(x,y)
#y <- c(2,4); x <- c(2,3); lines(x,y)
#x <- c(0,2); y <- -13*(x-2)/3+2; lines(x,y)
#x <- c(2,5); y <- -(x-2)/11+2; lines(x,y)
#text(x=1.3,y=2,'g=3',cex=2)
#text(x=2.2,y=1.4,'g=1',cex=2)
#text(x=2.6,y=2.2,'g=2',cex=2)
#text(x=2.1,y=2.6,'g=4',cex=2)
#dev.off()


W <- matrix(c(0.8,-0.6,0.9,0.5,-0.8,0.6,0.5,0.9),byrow=TRUE,ncol=2) # random sample from sphere
pdf('struct_linear.pdf')
plot(x=NULL,y=NULL,xlim=c(0,1),ylim=c(0,1),xlab="X1",ylab="X2",main=NULL,xaxs="i", yaxs="i")
y <- c(0,2); x <- 1.2/1.6*(y-2)+2; lines((x-1)/2,(y-1)/2)
y <- c(2,4); x <- c(2,3); lines((x-1)/2,(y-1)/2)
x <- c(0,2); y <- -13*(x-2)/3+2; lines((x-1)/2,(y-1)/2)
x <- c(2,5); y <- -(x-2)/11+2; lines((x-1)/2,(y-1)/2)
text(x=0.2,y=0.5,'g=3',cex=2)
text(x=0.6,y=0.2,'g=1',cex=2)
text(x=0.8,y=0.6,'g=2',cex=2)
text(x=0.55,y=0.8,'g=4',cex=2)
dev.off()

#----------
# nonlinear
pdf('struct_nonlinear.pdf')
plot(x=NULL,y=NULL,xlim=c(0,1),ylim=c(0,1),xlab="X1",ylab="X2",main=NULL,xaxs="i", yaxs="i")
radius <- 0.75; 
theta <- seq(0, 0.5 * pi, length = 200); lines(x = radius * cos(theta), y = radius * sin(theta))
theta <- seq(1.5*pi, 2 * pi, length = 200); lines(x = radius * cos(theta), y = radius * sin(theta)+1)
text(x=0.5,y=0.2,'g=1',cex=2)
text(x=0.5,y=0.8,'g=2',cex=2)
text(x=0.8,y=0.5,'g=3',cex=2)
text(x=0.2,y=0.5,'g=4',cex=2)
dev.off()

#----------
# asym
pdf('struct_asym.pdf')
plot(x=NULL,y=NULL,xlim=c(0,1),ylim=c(0,1),xlab="X1",ylab="X2",main=NULL,xaxs="i", yaxs="i")
abline(v=0.75); segments(0,0.33,0.75,0.33); segments(0,0.67,0.75,0.67);
text(x=0.87,y=0.5,'g=1',cex=2)
text(x=0.4,y=0.16,'g=2',cex=2)
text(x=0.4,y=0.49,'g=3',cex=2)
text(x=0.4,y=0.82,'g=4',cex=2)
dev.off()

#----------
# Quadratic
pdf('struct_quad.pdf')
plot(x=NULL,y=NULL,xlim=c(0,1),ylim=c(0,1),xlab="X1",ylab="X2",main=NULL,xaxs="i", yaxs="i")
abline(v=0.25); abline(v=0.5); abline(v=0.75); 
text(x=0.125,y=0.5,'g=1',cex=2)
text(x=0.375,y=0.5,'g=4',cex=2)
text(x=0.625,y=0.5,'g=2',cex=2)
text(x=0.875,y=0.5,'g=3',cex=2)
dev.off()




#----------
# Random stuff
load('../../explore/simulation/RData/main/Nsub_500_censor_0_dist_exponential_struct_tree_member_partition_alg_jlcmm_inter_TRUE_continuous_TRUE_ng_6_sim_1.RData')
load('../../explore/simulation/RData/main/Nsub_500_censor_0_dist_exponential_struct_asym_member_partition_alg_jlcmm_inter_TRUE_continuous_TRUE_ng_6_sim_1.RData')
load('../../explore/simulation/RData/main/Nsub_500_censor_0_dist_exponential_struct_quad_member_partition_alg_jlcmm_inter_TRUE_continuous_TRUE_ng_6_sim_1.RData')

X1 <- runif(10000,min=0,max=1)
X2 <- runif(10000,min=0,max=1)
X3 <- as.numeric(runif(10000)>0.5)
X4 <- round(runif(10000),1)
X5 <- sample(c(1:5),10000,replace=TRUE)
X <- cbind(X1,X2,X3,X4,X5)

pred_gg <- apply(predict_class(m6,X),1,which.max)
gg <- get_latent_class(X1,X2,'quad', 'partition', seed=1)
table(gg,pred_gg)

g <- ggplot()
g <- g + geom_point(aes(x=X1,y=X2,color=factor(pred_gg)))
print(g)


obj <- m6
nclasses <- ncol(obj$pprob)-2
coefs <- obj$best
coefend <- min(which(grepl('Weibull',names(coefs))))-1
coefs_multilogit <- matrix(coefs[1:coefend],nrow=nclasses-1)
tmpX <- cbind(1,newdata[,c('X1','X2','X3','X4','X5')])
if (ncol(coefs_multilogit)==7){ tmpX <- cbind(tmpX,newdata[,'X1']*newdata[,'X2']) }

linearval <- as.matrix(tmpX) %*% t(coefs_multilogit)
test_class <- t(apply(linearval, 1,function(x){ exps=exp(c(x,0)); exps/sum(exps)}))


data <- data.frame(cbind(X1,X2, y = as.numeric(linearval[,5]>0)))
g <- ggplot(data)
g <- g + geom_point(aes(x=X1,y=X2,color=factor(y)))
print(g)



X1 <- runif(10000,min=0,max=1)
X2 <- runif(10000,min=0,max=1)
gg <- get_latent_class(X1,X2,'asym', 'partition', seed=1)
data <- data.frame(cbind(X1,X2,gg))
ret <- rpart(gg ~ X1+X2, data)
pdf('asym_tree.pdf')
plot(ret,uniform=TRUE)
text(ret)
dev.off()
