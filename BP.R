eta <- 1

'<-'(
  sigmoid,
  function(x)
  {
    return(1/(1+exp(-x)))
  }
)

'<-'(
  prdict,
  function(dat, paras, fitting=FALSE)
  {
    '<-'(
      insdel,
      function(i)
      {
        return(paras$p1[i,] * dat[i])
      }
    )
    ins <- lapply(1:length(dat), insdel)
    ins <- as.matrix(ins)
    ins <- apply(ins,1,unlist)
    hideIns <- unlist(apply(ins, 1, sum))
    '<-'(
      hideCal,
      function(i)
      {
        return(sigmoid(hideIns[i] - paras$p2[i]))
      }
    )
    hideAns <- unlist(lapply(1:length(hideIns), hideCal))
    '<-'(
      hideDel,
      function(i)
      {
        return(hideAns[i] * paras$p3[i,])
      }
    )
    hideOuts <- lapply(1:length(hideAns), hideDel)
    hideOuts <- as.matrix(hideOuts)
    hideOuts <- apply(hideOuts,1,unlist)
    if(is.null(dim(hideOuts)))
    {
      outIns <- sum(hideOuts)
    }
    else
    {
      outIns <- unlist(apply(hideOuts, 1, sum)) 
    }
    '<-'(
      outDel,
      function(i)
      {
        if(length(outIns)==1)
        {
          return(outIns[i]-paras$p4[i])
        }
        return(sigmoid(outIns[i] - paras$p4[i]))
      }
    )
    outPuts <- unlist(lapply(1:length(outIns), outDel))
    if(fitting)
    {
      ans <- list(outPuts, hideAns)
      names(ans) <- c('yp', 'bhs')
      return(ans)
    }
    return(outPuts)
  }
)

'<-'(
  BPCAL,
  function(dat, paras, yr)
  {
    temp <- prdict(dat, paras, TRUE)
    yp <- temp$yp
    bhs <- temp$bhs
    gs <- yp * (1 - yp) * (yr - yp) #vector
    temp_bhs <- bhs * (1 - bhs)
    temp_ehs <- t(t(paras$p3) * gs)
    temp_ehs <- unlist(apply(temp_ehs, 1, sum))
    ehs <- temp_ehs * temp_bhs
    temp_gs <- eta * gs
    whs <- lapply(bhs, '*', temp_gs)
    whs <- as.matrix(whs)
    whs <- t(apply(whs,1,unlist))
    outThres <- -eta * gs
    hideThres <- -eta * ehs
    temp_ehs <- eta * ehs
    vhs <- as.matrix(lapply(dat, '*', temp_ehs))
    vhs <- t(apply(vhs,1,unlist))
    '<-'(
      EkDel,
      function(i)
      {
        return((yp[i]-yr[i])^2)
      }
    )
    Ek <- 0.5 * sum(unlist(lapply(1:length(yp), EkDel)))
    ans <- list(gs, ehs, whs, vhs, outThres, hideThres, Ek)
    names(ans) <- c('gs', 'ehs', 'whs', 'vhs', 'outThres', 'hideThres', 'Ek')
    return(ans)
  }
)

'<-'(
  BPSTD,
  function(dats, paras, rights)
  {
    '<-'(
      oneFit,
      function(i)
      {
        upd <- BPCAL(dats[i,], paras, rights[i,])
        paras[[1]] <<- paras$p1 + upd$vhs
        paras[[2]] <<- paras$p2 + upd$hideThres
        if(length(rights[i,])==1)
        {
          paras[[3]] <<- paras$p3 + t(upd$whs)
        }
        else
        {
          paras[[3]] <<- paras$p3 + upd$whs 
        }
        paras[[4]] <<-  paras$p4 + upd$outThres
        return(upd$Ek)
      }
    )
    Eks <- unlist(lapply(1:length(dats[,1]), oneFit))
    ans <- list(paras, Eks)
    names(ans) <- c('paras', 'Eks')
    return(ans)
  }
)

'<-'(
  BPADD,
  function(dats, paras, rights)
  {
    temp <- as.matrix(apply(dats, 1, prdict))
    yp <- t(yp)
    yr <- dats[,1]
    '<-'(
      gjCal,
      function(i)
      {
        
      }
    )
    #it is not complied
  }
)

'<-'(
  mainFit,
  function(num, dats, hides, outs)
  {
    recorder <- c()
    p1 <- matrix(runif(length(dats[1,])*hides), length(dats[1,]), hides)
    p2 <- runif(hides)
    p3 <- matrix(runif(hides*length(outs[1,])), hides, length(outs[1,]))
    p4 <- runif(length(outs[1,]))
    paras <- list(p1, p2, p3, p4)
    names(paras) <- c('p1','p2','p3','p4')
    '<-'(
      doFit,
      function(i)
      {
        temp <- BPSTD(dats, paras, outs)
        paras <<- temp$paras
        recorder <<- c(recorder, mean(temp$Eks))
      }
    )
    temp.no <- lapply(1:num, doFit)
    plot(1:num, recorder, type = 'o', main='loss-times', xlab='time',ylab='loss')
    return(paras)
  }
)

'<-'(
  testModel, 
  function(paras, dats, rights)
  {
    '<-'(
      doTest,
      function(i)
      {
        temp <- prdict(dats[i,], paras)
        temp <- which.max(temp)
        if(rights[i,temp]==1)
        {
          return(TRUE)
        }
        print(temp)
        print(rights[i,][tr[i,]==1])
        return(FALSE)
      }
    )
    testAns <- unlist(lapply(1:length(dats[,1]), doTest))
    return(length(testAns[testAns==TRUE])/length(dats[,1]))
  }
)
