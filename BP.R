eta <- 0.1

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
        return(dat[i] * paras[[1]][i,])
      }
    )
    ins <- apply(1:length(dat), insdel)
    ins <- as.numeric(as.matrix(ins))
    hideIns <- unlist(apply(ins, 2, sum))
    '<-'(
      hideCal,
      function(i)
      {
        return(sigmoid(hideIns[i] - paras[2][i]))
      }
    )
    hideAns <- unlist(lapply(1:length(hideIns), hideCal))
    '<-'(
      hideDel,
      function(i)
      {
        return(hideAns[i] * paras[[3]][i,])
      }
    )
    hideOuts <- lapply(1:length(hideAns), hideDel)
    hideOuts <- as.numeric(as.matrix(hideOuts))
    outIns <- unlist(apply(hideOuts, 1, sum))
    '<-'(
      outDel,
      function(i)
      {
        return(sigmoid(outIns[i] - paras[4][i]))
      }
    )
    outPuts <- unlist(lapply(1:length(outIns), outDel))
    if(fitting)
    {
      return(list(outPuts, hideAns))
    }
    return(outPuts)
  }
)

'<-'(
  BPCAL,
  function(dat, paras, yr)
  {
    temp <- prdict(dat, paras, TRUE)
    yp <- temp[1]
    bhs <- temp[2]
    gs <- yp * (1 - yp) * (yr - yp) #vector
    temp_bhs <- bhs * (1 - bhs)
    temp_ehs <- unlist(apply(paras[[3]], 1, sum))
    ehs <- unlist(lapply(temp_ehs, '*', temp_bhs))
    temp_gs <- eta * gs
    whs <- lapply(bhs, '*', temp_gs)
    whs <- as.matrix(whs)
    outThres <- -eta * gs
    hideThres <- -eta * ehs
    temp_ehs <- eta * ehs
    vhs <- as.matrix(lapply(dat, '*', temp_ehs))
    '<-'(
      EkDel,
      function(i)
      {
        return((yp[i]-yr[i])^2)
      }
    )
    Ek <- 0.5 * unlist(lapply(1:length(yp), EkDel))
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
      function(dat, paras,  right)
      {
        upd <- BPCAL(dat, paras, right)
        paras[[1]] <<- paras[[1]] + upd$vhs
        paras[2] <<- paras[2] + upd$hideThres
        paras[[3]] <<- paras[3] + upd$whs
        paras[4] <<-  paras[4] + upd$outThres
        return(upd$Ek)
      }
    )
    Eks <- apply(dats, 1, oneFit)
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
    
  }
)
