#liang
#2022-5-20

forest <- sapply(1:100,function(i){return(runif(100,10+0.5*i,10+2*i))})
trees <- sample(1:100,10)
heights <- round(runif(10,25,80))
spaces <- matrix(0,100,100)
for(x in heights)
{
  spaces[trees,1:x] <- 1
  forest[trees,1:x] <- forest[trees,1:x] * 0.6
}
spaces[trees,1] <- 2

hsy <- numeric(75)
for(it in 1:75)
{
  cnt <- 0
  for(i in trees)
  {
    cnt <- cnt + 1
    x <- 1
    while(spaces[i,x]==2)
    {
      x <- x + 1
    }
    x <- x - 1
    growstem <- 0.4
    if(runif(1) < growstem)
    {
      u <- i
      while(u <= 99 && (spaces[u,x]==2 || spaces[u,x]==3))
      {
        u <- u + 1
      }
      if(spaces[u,x]==0)
      {
        spaces[u,x] <- 3
      }
      else
      {
        u <- i
        while(u > 0 && (spaces[u,x]==2 || spaces[u,x]==3))
        {
          u <- u - 1
        }
        if(spaces[u,x]==0)
        {
          spaces[u,x] <- 3
        }
      }
    }
    else
    {
      if(x < heights[cnt])
      {
        spaces[i,x+1] <- 2
      }
    }
  }
  
  hsy[it] <- sum(forest[which(spaces==3)])
}
plot(1:75,hsy,'l')

plot(1:75,a1,'l',col=2,xlab='days',ylab='benfit')
lines(1:75,a2,'l',col=3)
legend('topleft',inset=.05, title="APT",c('p=0.6','p=0.4'),lty=c(1,1),col=c(2,3))