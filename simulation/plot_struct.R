library(optparse)
library(data.table)
library(plyr)
library(ggplot2)
library(latex2exp)

source('util.R')

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
par(mar=c(4,4.7,4,4))
plot(x=NULL,y=NULL,xlim=c(0,1),ylim=c(0,1),
#     xlab="X1",ylab="X2",
     xlab=expression(X[1]), ylab=expression(X[2]),
     main=NULL,xaxs="i", yaxs="i",
     ,cex.axis=2,cex.lab=2)

abline(v=0.5); abline(h=0.5)
text(x=0.25,y=0.75,'g=2',cex=3)
text(x=0.25,y=0.25,'g=1',cex=3)
text(x=0.75,y=0.25,'g=3',cex=3)
text(x=0.75,y=0.75,'g=4',cex=3)
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
par(mar=c(4,4.7,4,4))
plot(x=NULL,y=NULL,xlim=c(0,1),ylim=c(0,1),
#     xlab="X1",ylab="X2",
     xlab=expression(X[1]), ylab=expression(X[2]),
     main=NULL, xaxs="i", yaxs="i", cex.axis=2, cex.lab=2)
y <- c(0,2); x <- 1.2/1.6*(y-2)+2; lines((x-1)/2,(y-1)/2)
y <- c(2,4); x <- c(2,3); lines((x-1)/2,(y-1)/2)
x <- c(0,2); y <- -13*(x-2)/3+2; lines((x-1)/2,(y-1)/2)
x <- c(2,5); y <- -(x-2)/11+2; lines((x-1)/2,(y-1)/2)
text(x=0.2,y=0.5,'g=3',cex=3)
text(x=0.6,y=0.2,'g=1',cex=3)
text(x=0.8,y=0.6,'g=2',cex=3)
text(x=0.55,y=0.8,'g=4',cex=3)
dev.off()

#----------
# nonlinear
pdf('struct_nonlinear.pdf')
par(mar=c(4,4.7,4,4))
plot(x=NULL,y=NULL,xlim=c(0,1),ylim=c(0,1),
#     xlab="X1",ylab="X2",
     xlab=expression(X[1]), ylab=expression(X[2]),
     main=NULL,xaxs="i", yaxs="i"
     ,cex.axis=2,cex.lab=2)
radius <- 0.75; 
theta <- seq(0, 0.5 * pi, length = 200); lines(x = radius * cos(theta), y = radius * sin(theta))
theta <- seq(1.5*pi, 2 * pi, length = 200); lines(x = radius * cos(theta), y = radius * sin(theta)+1)
text(x=0.5,y=0.2,'g=1',cex=3)
text(x=0.5,y=0.8,'g=2',cex=3)
text(x=0.8,y=0.5,'g=3',cex=3)
text(x=0.2,y=0.5,'g=4',cex=3)
dev.off()

#----------
# asym
pdf('struct_asym.pdf')
par(mar=c(4,4.7,4,4))
plot(x=NULL,y=NULL,xlim=c(0,1),ylim=c(0,1),
#     xlab="X1",ylab="X2",
     xlab=expression(X[1]), ylab=expression(X[2]),
     main=NULL,xaxs="i", yaxs="i",
     cex.axis=2,cex.lab=2)
abline(v=0.75); segments(0,0.33,0.75,0.33); segments(0,0.67,0.75,0.67);
text(x=0.87,y=0.5,'g=1',cex=3)
text(x=0.4,y=0.16,'g=2',cex=3)
text(x=0.4,y=0.49,'g=4',cex=3)
text(x=0.4,y=0.82,'g=3',cex=3)
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




