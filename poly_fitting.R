library(parallel)
source('./RStudy/TSP.R')

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
    pts <- cbind((pts %% height) + 1, floor(pts/width))
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
  getEdge,
  function(pts)
  {
    len <- dim(pts)[1]
    tmp <- matrix(0,max(pts[,1]),max(pts[,2]))
    tmp[((pts[,2]-1)*max(pts[,1]) + pts[,1])] <- 1
    t1 <- sapply(1:max(pts[,1]), function(i){return(which.min(which(tmp[i,]==1)))})
    t2 <- sapply(1:max(pts[,1]), function(i){return(which.max(which(tmp[i,]==1)))})
    t3 <- sapply(1:max(pts[,2]), function(i){return(which.min(which(tmp[,i]==1)))})
    t4 <- sapply(1:max(pts[,2]), function(i){return(which.max(which(tmp[,i]==1)))})
    tz <- c(t1,t2)
    ans1 <- matrix(c(1:length(tz),tz), length(tz), 2)
    tz <- c(t3,t4)
    ans2 <- matrix(c(tz,1:length(tz)), length(tz), 2)
    ans <- rbind(ans1, ans2)
    return(ans)
    #len <- dim(pts)[1]
    #if(len > 100)
    #{
    #  dm <- as.matrix(dist(pts))
    #  vote <- rowSums(dm)
    #  vote <- abs((vote / sum(vote)) - 0.1)
    #  pts <- pts[sample(1:len, 200, prob = vote),]
    #}
    #road <- SaTsp(dist(pts))
    road <- AcoTsp(dist(ans))
    edge <- ans[road,]
    return(edge)
  }
)

'<-'(
  ployFit,
  function(edge)
  {
    
  }
)
