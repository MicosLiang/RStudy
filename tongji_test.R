'<-'(
  chi_square,
  function(xv, xe=NULL, typ=0)
  {
    if(is.null(xe))
    {
      xe <- mean(xv)
    }
    else
    {
      if(typ==1)
      {
        xv <- xv / sum(xv)
      }
      else if(typ==2)
      {
        xe <- xe * sum(xv)
      } 
    }
    return(sum((xv-xe)^2/xe))
  }
)

'<-'(
  x2_live,
  function(x2,df)
  {
    k <- df/2
    return(0.5^k/my_gamma(k) * x2 ^ (k-1) * exp(1)^(-x2/2))
  }
)

'<-'(
  ingamma,
  function(z,s,isdown=F)
  {
    b <- 1/s*z^s*exp(1)^-z
    fs <- 1
    M <- 1
    for(i in 1:1000)
    {
      tmp <- M
      fs <- fs*(s+i)
      tmp <- tmp + z^i/fs
      if(!(tmp < Inf) && tmp>=0)
      {
        M <- tmp 
      }
      else
      {
        break
      }
    }
    #it goes well that it can cal p>0.005
    if(isdown)
    {
      return(my_gamma(s)-b*M)
    }
    return(b*M)
  }
)

'<-'(
  my_gamma,
  function(z)
  {
    #it seems well when dof<100
    return(sqrt(2*pi/z)*((1/exp(1))*(z+1/((12*z)-(1/10*z))))^z)
  }
)

'<-'(
  chi_p,
  function(x,k)
  {
    return(ingamma(x/2,k/2,F)/my_gamma(k/2))
  }
)

'<-'(
  chi_test,
  function(x1, x2=NULL)
  {
    len1 <- length(x1)
    if(!is.null(x2))
    {
      len2 <- length(x2)
      if(len1==len2)
      {
        square <- chi_square(x1,x2)
      }
      else
      {
        warning('x1 and x2 must have same length!')
        return(c(1,0,1))
      }
    }
    else
    {
      square <- chi_square(x1)
    }
    
    ans <- numeric(3)
    ans[1] <- square
    ans[2] <- len1-1
    ans[3] <- chi_p(square,len1-1)
    names(ans) <- c('X^2','dof','p-value')
    return(ans)
  }
)

'<-'(
  odist,
  function(v1,v2)
  {
    return(sqrt(sum((v1-v2)^2)))
  }
)

'<-'(
  class_in,
  function(yang, zu)
  {
    return(mean(unlist(apply(zu,1,odist,yang))))
  }
)

'<-'(
  class_out,
  function(yang, zu)
  {
    return(min(unlist(apply(zu,1,odist,yang))))
  }
)

'<-'(
  class_right,
  function(yang, zu1, zu2)
  {
    ai <- class_in(yang,zu1)
    bi <- class_out(yang,zu2)
    return((bi-ai)/max(ai,bi))
  }
)

'<-'(
  sil_coe,
  function(zu1, zu2)
  {
    nr <- dim(zu1)[1]
    tmp <- 0
    for(i in 1:nr)
    {
      tmp <- tmp + class_right(zu1[i,],zu1[-i,],zu2)
    }
    ans <- c(tmp/nr)
    names(ans) <- 'Silhouette Coefficient '
    return(ans)
  }
)

'<-'(
  gradient,
  function(x,y,eta=0.0001,iter=100)
  {
    #x <- (x - mean(x)) / sd(x)
    d <- dim(x)[2]
    if(length(d)==0)
    {
      d <- 1
    }
    w <- rnorm(d+1)
    x <- cbind(1,x)
    for(i in 1:iter)
    {
      yp <- x%*%w
      ye <- y - yp
      j <- eta*colSums(apply(x,2,function(xi){return(xi*ye)}))
      w <- w + j
    }
    return(w)
  }
)

'<-'(
  my_lm,
  function(x,y,d=1,draw=F)
  {
    ox <- x
    x <- sapply(0:d,function(v){return(x^v)})
    e <- try(w <- solve((t(x)%*%x))%*%(t(x)%*%y))
    if(typeof(e)=='character')
    {
      warning('not good dim')
      return(x,y,d-1,draw)
    }
    if(draw)
    {
      plot(ox,y,col=2,xlab='x',ylab='y',main=paste('lm: dim=',d))
      tx <- seq(min(ox)-10,max(ox)*1.2,0.1)
      x <- sapply(0:d,function(v){return(tx^v)})
      lines(tx, as.vector(x%*%w), col=3, type='l')
    }
    return(w)
  }
)

'<-'(
  my_lm2,
  function(x,y)
  {
    x <- cbind(1,x)
    w <- solve((t(x)%*%x))%*%(t(x)%*%y)
    return(w)
  }
)

'<-'(
  lm_R,
  function(x,y,w)
  {
    d <- length(w)-1
    x <- sapply(0:d,function(v){return(x^v)})
    yp <- x%*%w
    sr <- sum((yp-mean(y))^2)
    se <- sum((y-yp)^2)
    ans <- c(sr/(sr+se))
    names(ans) <- 'R'
    return(ans)
  }
)

'<-'(
  my_t,
  function(x,n)
  {
    t <- my_gamma((n+1)/2)/my_gamma(n/2)*(1+(x^2)/n)^(-1*((n+1)/2))
    return(t)
  }
)

'<-'(
  lm_t,
  function(x,y,w)
  {
    x <- cbind(1,x)
    yp <- x%*%w
    n <- dim(x)[1]
    p <- length(w)
    c <- solve(t(x)%*%x)
    tj <- numeric(p)
    s <- n-p-1
    t_sd <- sqrt((1/s)*sum((y-yp)^2))
    for(i in 1:p)
    {
      tj[i] <- -1 * w[i] / sqrt(c[i,i]*t_sd)
    }
    tj <- abs(tj)
    return(1-pt(tj,s))
  }
)

#
alpha = 0.05

#first
print('(1 ------------------------------')
print(dat1)
ans1 <- chi_test(dat1[,1],dat1[,2])
print('use Chi square test: ')
print(ans1)
if(ans1[3]>alpha)
{
  print('Because the p-value of data is lager than alpha, so we think that H0 is true')
}

#second
print('(2 ------------------------------')
print(dat2)
plot(1:10,dat2[,1],'b',col=2,ylim=c(0,max(dat2)*1.2),xlab='tag',ylab='weight',main='Effect of feed on mice')
lines(1:10,dat2[,2],col=3,'b')
legend('topright',pch=c(15,15),legend=c('control','experiment'),col=c(2,3),bty='n')
ans2 <- chi_test(dat2[,1],dat2[,2])
print('use Chi square test: ')
print(ans2)
if(ans2[3]<=alpha)
{
  print('Because the p-value of data is small than alpha, so we think that H0 is true')
}

#third
print('(3 ------------------------------')
print(dat3)
plot(c(0.5,1,2),dat3[1,-1],'b',col=2,ylim=c(0,max(dat3)*1.2),xlab='dose',ylab='weight',main='Dimension reduction Plot')
for(i in 2:10)
{
  lines(c(0.5,1,2),dat3[i,-1],'b',col=2)
}
for(i in 11:20)
{
  lines(c(0.5,1,2),dat3[i,-1],'b',col=3)
}
legend('topright',pch=c(15,15),legend=c('juice','vitman c'),col=c(2,3),bty='n')

#y~ w0 + w1x1 + w2x2 + w3x1x2
#use multiple linear regression, consider fees_way = juice = 0, vitamen = 1
pts_1 <- c(c(rep(0,30),rep(1,30)),c(c(rep(0.5,10),rep(1,10),rep(2,10)),c(rep(0.5,10),rep(1,10),rep(2,10))))
pts_1 <- matrix(pts_1, 60, 2)
pts_1 <- cbind(pts_1, pts_1[,1] * pts_1[,2])
pts_1 <- cbind(pts_1, c(dat3[1:10,2],dat3[11:20,2],dat3[1:10,3],dat3[11:20,3],dat3[1:10,4],dat3[11:20,4]))
#ans3_4 <- numeric(4)
#for(i in 1:10)
#{
#  ans3_4 <- ans3_4 + gradient(pts_1[,1:3],pts_1[,4])
#}
#ans3_4 <- ans3_4 / 10
ans3_4 <- my_lm2(pts_1[,1:3],pts_1[,4])
model_3 <- paste('model : weight = ', ans3_4[1],'+',ans3_4[2],'* feed_way +', ans3_4[3],'* dose +',ans3_4[4],'* (feed_way * dose)')
print(model_3)
ans3_5 <- lm_t(pts_1[,1:3],pts_1[,4],ans3_4)
names(ans3_5) <- c('no-meaning','feed_way','dose','feed_way*dose')
print('p-value: ')
print(ans3_5)
if(ans3_5[2] < alpha)
{
  print('The feed way have the effect to weight')
}
if(ans3_5[3] < alpha)
{
  print('The dose have the effect to weight')
}
if(ans3_5[4] < alpha)
{
  print('The feed way * dose have the effect to weight')
}

#draw picture
plot(pts_1[1:30,2],pts_1[1:30,4],col=2,xlab='dose/mg',ylab='weight',main='Multiple linear regression')
lines(pts_1[31:60,2],pts_1[31:60,4],col=3,type='p')
ans3_6 <- my_lm(pts_1[1:30,2],pts_1[1:30,4])
ans3_7 <- my_lm(pts_1[31:60,2],pts_1[31:60,4])
ans3_8 <- my_lm(pts_1[,2],pts_1[,4])
#lines(c(0.5,2),c(c(1,0.5)%*%ans3_6,c(1,2)%*%ans3_6),col=2,type='l')
#lines(c(0.5,2),c(c(1,0.5)%*%ans3_7,c(1,2)%*%ans3_7),col=3,type='l')
lines(c(0.5,2),c(c(1,0.5)%*%ans3_8,c(1,2)%*%ans3_8),col=4,type='l')
lines(c(0.5,2),c(c(1,0,0.5,0)%*%ans3_4,c(1,0,2,0)%*%ans3_4),col=5,type='l')
lines(c(0.5,2),c(c(1,1,0.5,0.5)%*%ans3_4,c(1,1,2,2)%*%ans3_4),col=6,type='l')
legend('bottomright',pch=c(15,15,15,15),legend=c('juice','vitamen c','all_base','feed_way=jucie','feed_way=vitamen c'),col=c(2,3,4,5,6),bty='n')

#fourth
print('(4 ------------------------------')
ans4 <- my_lm(dat4[,1], dat4[,2], d = 2, draw = T)
model_4 <- paste('model : weight =',ans4[1],'+',ans4[2],'* height +',ans4[3],'* height^2')
print(model_4)
ans4_1 <- lm_R(dat4[,1], dat4[,2], ans4)
print(ans4_1)
if(ans4_1>0.9)
{
  print('The ploy fit result is very well')
}