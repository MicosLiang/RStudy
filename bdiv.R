rm(list=ls())

'<-'(
  aDiv,
  function(spv)
  {
    return(length(which(spv==T)))
  }
)

'<-'(
  gDiv,
  function(spm)
  {
    return(length(which(sapply(1:dim(spm)[2],function(x){return(any(spm[,x]))}))))
  }
)

'<-'(
  bDiv,
  function(mean_alpha, gam)
  {
    return(1-mean_alpha/gam)
  }
)

'<-'(
  juan, #一个类似卷积的操作
  function(dm,func)
  {
    len <- dim(dm)[1]
    ans <- NULL
    ans_scale <- character()
    ans_id <- character()
    al_cnt <- 1
    for(scale in 1:(len-1))
    {
      square <- (scale+1)^2
      end <- len - scale
      cnt <- 1
      for(posx in end:1)
      {
        for(posy in 1:end)
        {
          edgex <- posx + scale
          edgey <- posy + scale
          ans[al_cnt] <- func(dm[posx:edgex,posy:edgey])
          #ans_scale[al_cnt] <- paste('scale',scale,sep="",collaspe="")
          ans_scale[al_cnt] <- square
          ans_id[al_cnt] <- paste('i',cnt,sep="",collapse = "")
          cnt <- cnt + 1
          al_cnt <- al_cnt + 1
        }
      }
    }
    ans <- cbind(ans_scale,ans_id,ans)
    return(ans)
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
      plot(ox,y,col=3,xlab='x',ylab='y',main=paste('lm: dim=',d))
      tx <- seq(min(ox)-10,max(ox)*1.2,0.1)
      x <- sapply(0:d,function(v){return(tx^v)})
      lines(tx, as.vector(x%*%w), col=2, type='l')
      legend('topright',legend=c('y=Σki*x^i',w))
    }
    return(w)
  }
)

'<-'(
  my_coef,
  function(x,y,do_sc=T)
  {
    if(do_sc)
    {
      #标准化
      x <- scale(x)
      y <- scale(y)
    }
    est <- my_lm(x,y)
    #attr(est, 't') <- lm_t(x,y,est)[-1]
    est <- est[-1]
    return(est)
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

dat <- read.csv("./RStudy/data/bdiv.csv")

spMatrix <- as.matrix(dat[,-(1:6)])

grids <- levels(factor(dat$grid))
plotId <- levels(factor(dat$plot))
plot_len <- 8

datByGrid <- lapply(grids, function(grid){
  pid <- which(dat$grid==grid)
  pdat <- dat[pid,]
  pspMatrix <- spMatrix[pid,]
  adiv <- matrix(NA,plot_len,plot_len,dimnames = list(1:plot_len,LETTERS[1:plot_len]))
  newSpm <- matrix(NA,plot_len,plot_len,dimnames = list(1:plot_len,LETTERS[1:plot_len]))
  biomass <- matrix(NA,plot_len,plot_len,dimnames = list(1:plot_len,LETTERS[1:plot_len]))
  tot.biomass <- matrix(NA,plot_len,plot_len,dimnames = list(1:plot_len,LETTERS[1:plot_len]))
  cnt <- 1
  for(x in plotId)
  {
    v <- which(pdat$plot==x)
    #if(length(v)!=1)
    #{
    #  warning('data error!')
    #}
    adiv[cnt] <- aDiv(pspMatrix[v,])
    newSpm[cnt] <- v
    biomass[cnt] <- pdat[v,]$biomass
    tot.biomass[cnt] <- pdat[v,]$tot.bio
    cnt <- cnt + 1
  }
  adiv <- juan(adiv,mean)
  gdiv <- juan(newSpm,function(spvId){
    return(gDiv(pspMatrix[as.vector(spvId),]))
  })
  bdiv <- 1 - as.numeric(adiv[,3]) / as.numeric(gdiv[,3])
  bio_mean <- juan(biomass,mean)
  bio_sum <- juan(biomass,sum)
  tot_bio_mean <- juan(tot.biomass,mean)
  tot_bio_sum <- juan(tot.biomass,sum)
  scale <- adiv[,1]
  isubsample <- adiv[,2]
  grid <- rep(grid, length(scale))
  adiv <- adiv[,3]
  gdiv <- gdiv[,3]
  bio_mean <- bio_mean[,3]
  bio_sum <- bio_sum[,3]
  tot_bio_mean <- tot_bio_mean[,3]
  tot_bio_sum <- tot_bio_sum[,3]
  return(cbind(grid, scale,isubsample,adiv,gdiv,bdiv,bio_mean,bio_sum,tot_bio_mean,tot_bio_sum))
})

final_ans <- NULL
for(each in datByGrid)
{
  final_ans <- rbind(final_ans, each)
}
#final_ans即所求

final_ans2 <- final_ans[,-3]
final_ans2 <- apply(final_ans2,2,as.numeric)
datByGrid <- lapply(as.numeric(grids), function(grid){
  pid <- which(final_ans2[,1]==grid)
  pdat <- final_ans2[pid,]
  scale <- levels(factor(pdat[,2]))
  ans <- sapply(scale, function(l){
    w <- which(pdat[,2] == l)
    if(length(w)>1)
    {
      return(unlist(apply(pdat[w,-(1:2)], 2, mean)))
    }
    else
    {
      return(pdat[w,-(1:2)])
    }
  })
  ans <- t(ans)
  grid <- rep(grid,dim(ans)[1])
  ans <- cbind(grid,scale,ans)
  return(ans)
})

final_ans2 <- NULL
for(each in datByGrid)
{
  final_ans2 <- rbind(final_ans2, each)
}
final_ans2 <- apply(final_ans2,2,as.numeric)

est <- sapply(levels(factor(final_ans2[,2])),function(sc){
  pid <- which(final_ans2[,2]==sc)
  pdat <- final_ans2[pid,]
  beta_bio <- my_coef(pdat[,5],log(pdat[,9]),T)
  alpha_bio <- my_coef(pdat[,3],log(pdat[,9]),T)
  gamma_bio <- my_coef(pdat[,4],log(pdat[,9]),T)
  print(attr(beta_bio, 't'))
  #gamma_bio <- my_coef(pdat[,4],log(pdat[,9]),F)
  #biomass~gamma的标准化和非标准化系数有极大差异
  return(c(alpha_bio,beta_bio,gamma_bio))
})

est <- rbind(levels(factor(final_ans2[,2])),est)
par(mfrow=c(1,3))
plot(est[1,],est[2,],'b',main='biomass~alphaDiv',xlab='scale(m^2)',ylab='std.Estimate',col=2)
plot(est[1,],est[3,],'b',main='biomass~betaDiv',xlab='scale(m^2)',ylab='std.Estimate',col=4)
plot(est[1,],est[4,],'b',main='biomass~gammaDiv',xlab='scale(m^2)',ylab='std.Estimate',col=6)