library(optparse)
library(data.table)
library(plyr)
library(lcmm)
source('0.init.R')
source('0.gen_survival.R')
X1 <- runif(10000,min=1,max=3)
X2 <- runif(10000,min=1,max=3)
gg <- get_latent_class(X1,X2,'linear', 'partition', seed=sim)


data <- cbind(X1,X2,gg)
g <- ggplot()
g <- g + geom_point(aes(x=X1,y=X2,color=factor(gg)))
print(g)
W <- matrix(c(0.8,-0.6,0.9,0.5,-0.8,0.6,0.5,0.9),byrow=TRUE,ncol=2) # random sample from sphere



pdf('struct_linear.pdf')
plot(x=NULL,y=NULL,xlim=c(1,3),ylim=c(1,3),xlab="X1",ylab="X2",main=NULL,xaxs="i", yaxs="i")
y <- c(0,2); x <- 1.2/1.6*(y-2)+2; lines(x,y)
y <- c(2,4); x <- c(2,3); lines(x,y)
x <- c(0,2); y <- -13*(x-2)/3+2; lines(x,y)
x <- c(2,5); y <- -(x-2)/11+2; lines(x,y)
text(x=1.3,y=2,'g=3',cex=2)
text(x=2.2,y=1.4,'g=1',cex=2)
text(x=2.6,y=2.2,'g=2',cex=2)
text(x=2.1,y=2.6,'g=4',cex=2)
dev.off()

