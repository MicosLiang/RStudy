library(parallel)

'<-'(
  read.bmp24,
  function(fn)
  {
    fn <- file(fn, 'rb')
    temp.no <- readBin(fn, 'raw', 18)
    width <- as.integer(paste(c('0x', rev(as.character(readBin(fn,'raw',4)))),collapse = ''))
    height <- as.integer(paste(c('0x', rev(as.character(readBin(fn,'raw',4)))),collapse = ''))
    temp.no <- readBin(fn, 'raw', 28)
    img <- as.integer(readBin(fn, 'raw', width*height*3))
    close(fn)
    ans = matrix(0, width, height)
    img <- matrix(img, width*height, 3, byrow = TRUE)
    '<-'(
      Gray,
      function(RGB)
      {
        R <- RGB[1]
        G <- RGB[2]
        B <- RGB[3]
        g <- (sum(R^2.2, (1.5*G)^2.2, (0.6*B)^2.2) / 3.7651) ^ 0.4545455
        return(g)
      }
    )
    clist <- makeCluster(detectCores())
    img <- unlist(parApply(clist, img, 1, Gray))
    img <- matrix(img, width, height)
    stopCluster(clist)
    return(img)
  }
)

'<-'(
  waveFilter,
  function(img, cap= matrix(c(1,1,1,1,-8,1,1,1,1),3,3))
  {
    #cap <- matrix(c(1,1,1,1,-8,1,1,1,1),3,3)
    dims <- dim(img)
    width <- dims[1]
    height <- dims[2]
    ans <- matrix(0, width, height)
    ans[2:width,2:height] <- ans[2:width,2:height] + img[1:(width-1),1:(height-1)] * cap[1,1]
    ans[1:width,2:(height)] <- ans[1:(width),2:(height)] + img[1:(width),1:(height-1)] * cap[1,2]
    ans[1:(width-1),1:(height-1)] <- ans[1:(width-1),1:(height-1)] + img[2:(width),2:(height)] * cap[1,3]
    ans[2:(width),1:(height)] <- ans[2:(width),1:(height)] + img[1:(width-1),1:(height)] * cap[2,1]
    ans[2:(width-1),2:(height-1)] <- ans[2:(width-1),2:(height-1)] + img[2:(width-1),2:(height-1)] * cap[2,2]
    ans[1:(width-1),1:(height)] <- ans[1:(width-1),1:(height)] + img[2:width,1:(height)] * cap[2,3]
    ans[2:(width),1:(height-1)] <- ans[2:(width),1:(height-1)] + img[1:(width-1),2:height] * cap[3,1]
    ans[1:(width),1:(height-1)] <- ans[1:(width),1:(height-1)] + img[1:(width),2:height] * cap[3,2]
    ans[1:(width-1),1:(height-1)] <- ans[1:(width-1),1:(height-1)] + img[2:width,2:height] * cap[3,3]
    ans <- ans[2:(width-1),2:(height-1)]
    return(abs(ans))
  }
)

'<-'(
  getPoints,
  function(img, yuzhi=-1)
  {
    dims <- dim(img)
    width <- dims[2]
    height <- dims[1]
    yuzhi <- ifelse(yuzhi==-1,max(img)/2,yuzhi)
    print(yuzhi)
    pts <- which(img > yuzhi)
    pts <- rbind((pts %% height) + 1, floor(pts/width))
    return(pts)
  }
)


'<-'(
  lessNoise,
  function(img)
  {
    img <- waveFilter(img, matrix(c(1,2,1,2,4,2,1,2,1),3,3))
    img <- img / 16
    return(img)
  }
)

#The following is not completed


'<-'(
  outsideFind,
  function(img)
  {
    dims <- dim(img)
    rows <- dims[1]
    cols <- dims[2]
    nowP <- c(1,1)
    moves <- matrix(c(1,0,1,1,0,1,-1,-1,-1,0,-1,1,-1,0,1,-1),2,8)
    nowDc <- 0
    
  }
)

'<-'(
  ployFit,
  function(edge)
  {
    
  }
)
