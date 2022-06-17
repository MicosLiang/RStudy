'<-'(
  PCA,
  function(A, tag = NULL, center = T, normalize = F, d = 2, draw = T)
  {
    m <- dim(A)[1]
    n <- dim(A)[2]
    if(center)
    {
      A <- apply(A,2,function(x){return(x-mean(x))}) 
    }
    if(normalize)
    {
      A <- apply(A,2,function(x){return(scale(x))})
    }
    B <- 1/(m-1) * t(A) %*% A
    tmp <- eigen(B)
    values <- as.vector(tmp$values)
    disVaule <- values / sum(values)
    vectors <- as.matrix(tmp$vectors)
    vectors <- vectors[,1:d]
    C <- A %*% vectors
    if(d==2 && draw)
    {
      if(is.null(tag))
      {
        plot(C[,1],C[,2],xlab=paste('PC1(',disVaule[1],'%)'),ylab=paste('PC2(',disVaule[2],'%)'))
      }
      else
      {
        types <- levels(factor(tag))
        classes <- lapply(types, function(x){return(which(tag==x))})
        plot(C[classes[[1]],1],
             C[classes[[1]],2],
             col=2,
             xlab=paste('PC1(',round(disVaule[1]*100,2),'%)'),
             ylab=paste('PC2(',round(disVaule[2]*100,2),'%)'),
             xlim=c(min(C[,1]),max(C[,1])),
             ylim=c(min(C[,2]),max(C[,2])),
             pch=16)
        classes <- classes[-1]
        col_cnt <- 2
        for(each in classes)
        {
          col_cnt <- col_cnt + 1
          lines(C[each,1],C[each,2],col=col_cnt,pch=16,type='p')
        }
        legend('topleft',legend = types, col=2:col_cnt, pch=16)
      }
    }
  }
)