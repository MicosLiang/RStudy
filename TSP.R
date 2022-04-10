#library(Rcpp)
library(parallel)

'<-'(
  AcoTsp,
  function(dm)
  {
    if(class(dm) == 'dist')
    {
      dm <- as.matrix(dm)
    }
    
    nums <- dim(dm)[1]
    alpha <- 1
    beta <- 5
    Q <- 100
    rho <- 0.9 # 1 - 挥发
    pm <- matrix(1, nums, nums) #信息素矩阵
    em <- 1 / dm #距离启发矩阵
    ants <- nums #蚂蚁数量
    nIter <- 100 #迭代次数
    
    distances <- numeric(nIter) #每次迭代最小距离
    minRoad <- integer(nums) #最短路
    minDis <- Inf
    
    if(nums > 50)
    {
      #大于50就并行计算，因为线程启动也要时间
      clist <- parallel::makeCluster(parallel::detectCores(logical = FALSE))
    }
    
    '<-'(
      eachIter,
      function(it)
      {
        ants <- ceiling(nums * (1- (it / nIter))) + 10
        rm <- matrix(0, nums+nums+1, ants) #每只蚂蚁的路径矩阵（每列）,+nums是为了保存路径在距离矩阵中的位置和总距离
        am <- matrix(TRUE, nums, ants) #每只蚂蚁（每列）的允许访问表
        pmTmp <- (em ^ alpha) * (pm ^ beta) #计算转移概率用到的中间变量矩阵
        
        starts <- sample(1:nums, ants, replace = TRUE) #为蚂蚁选择起点
        rm[1,] <- starts
        am[starts + (0:(ants-1) * nums)] <- FALSE
        '<-'(
          eachAnt,
          function(k)
          {
            allows <- am[,k]
            road <- rm[,k]
            distance <- 0
            '<-'(
              findNext,
              function(i)
              {
                if(i==nums-1)
                {
                  return(which(allows))
                }
                n <- road[i] #当前位置
                allow <- which(allows)
                antPm <- pmTmp[n,allow]
                antPm <- antPm / sum(antPm)
                nxt <- sample(allow, 1, replace = FALSE, prob = antPm) #效果类似轮盘赌
                return(nxt)
              }
            )
            '<-'(
              goNext,
              function(i)
              {
                n <- road[i] #当前位置
                if(i<nums)
                {
                  nxt <- findNext(i)
                  distance <<- distance + dm[n,nxt]
                  road[i+1] <<- nxt
                  road[i+nums] <<- (ifelse(n<nxt,n,nxt)-1)*nums + ifelse(n<nxt,nxt,n) #释放信息素的位置，控制在下三角
                  allows[nxt] <<- FALSE
                  return(goNext(i+1))
                }
                nxt <- road[1]
                distance <<- distance + dm[n,nxt]
                road[i+nums] <<- (ifelse(n<nxt,n,nxt)-1)*nums + ifelse(n<nxt,nxt,n)
                return(TRUE)
              }
            )
            try(goNext(1))
            road[nums+nums+1] <- distance
            return(road)
          }
        )
        if(nums > 50)
        {
          rm <- parallel::parSapply(clist, 1:ants, eachAnt)
        }
        else
        {
          rm <- sapply(1:ants, eachAnt)
        }
        wMin <- which.min(rm[nums+nums+1,])
        distances[it] <<- rm[nums+nums+1,wMin]
        if(minDis > rm[nums+nums+1,wMin])
        {
          minRoad <<- rm[1:nums,wMin]
          minDis <<- rm[nums+nums+1,wMin]
        }
        rm[nums+nums+1,] <- Q / rm[nums+nums+1,]
        pm <<- pm * rho #蒸发
        apply(rm, 2, function(x) {pm[x[(nums+1):(nums+nums)]] <<- pm[x[(nums+1):(nums+nums)]] + x[nums+nums+1]})
        pm <<- as.matrix(as.dist(pm))
      }
    )
    try(sapply(1:nIter, eachIter))
    plot(1:nIter, distances, 'l',xlab='iter_times',ylab='iter_distance',main='ACO-TSP')
    
    if(nums>50)
    {
      parallel::stopCluster(clist)
    }
    print(minDis)
    return(minRoad)
  }
)

'<-'(
  SaTsp,
  function(dm, TS = 100, TE = 20,alpha = 0.99)
  {
    if(class(dm) == 'dist')
    {
      #浅浅验证下输入
      dm <- as.matrix(dm)
    }
    
    '<-'(
      interfere,
      function(r)
      {
        #对解进行随机的扰动，随机采用不同的两种方式
        if(runif(1) > 0.5)
        {
          #交换两座城市的顺序
          k <- sample(1:nums, 2)
          r[k] <- r[rev(k)]
        }
        else
        {
          #交换一段城市的顺序
          k <- sort(sample(2:nums, 3))
          nr <- r
          r <- c(r[1:(k[1]-1)],r[k[2]:(k[3]-1)],r[k[1]:(k[2]-1)],r[(k[3]):nums])
        }
        return(r)
      }
    )
    
    '<-'(
      calDistance,
      function(r)
      {
        #计算路程距离，作为评估解好坏的指标
        tr <- integer(nums)
        tr <- r[2:nums]
        tr[nums] <- r[1]
        return(sum(dm[((r-1)*nums + tr)]))
      }
    )
    
    nums <- dim(dm)[1] #城市数量
    Tp <- TS #温度
    road <- sample(1:nums, nums) #随机产生路径作为初始解
    dis <- calDistance(road)
    displot <- numeric(floor(log(TE/TS, alpha))+10)
    
    cout <- 1
    while(Tp > TE)
    {
      cout <- cout + 1
      #你想用这种方式也行
      #cnt <- ceiling(nums*TS/Tp)
      cnt <- 10
      while(cnt > 0)
      {
        #扰动
        nroad <- interfere(road)
        ndis <- calDistance(nroad)
        #评价新解的好坏
        deta <- (ndis - dis)/dis
        #模拟退火的核心公式
        if(deta < 0 || exp(-Tp/deta) > runif(1))
        {
          dis <- ndis
          road <- nroad
          #cnt <- ceiling(nums*Tp/TS)
          cnt <- 10
        }
        #十次没有获得新解就降温
        cnt <- cnt - 1
      }
      Tp <- Tp * alpha #降温
      displot[cout] <- dis
    }
    
    plot(2:cout, displot[2:cout],'l',xlab='iter_times',ylab='iter_distance',main='SA-TSP')
    print(dis)
    return(road)
  }
)


'<-'(
  GATSP,
  function()
  {
    
  }
)

'<-'(
  testAco,
  function(pn=30)
  {
    px <- runif(pn,100,1000)
    py <- runif(pn,100,1000)
    pts <- matrix(c(px,py),pn,2)
    road <- AcoTsp(dist(pts))
    plot(c(pts[road,1], pts[road[1],1]),c(pts[road,2],pts[road[1],2]),'o',xlab='x',ylab='y')
  }
)

'<-'(
  testSa,
  function(pn=30)
  {
    px <- runif(pn,100,1000)
    py <- runif(pn,100,1000)
    pts <- matrix(c(px,py),pn,2)
    road <- SaTsp(dist(pts))
    plot(c(pts[road,1], pts[road[1],1]),c(pts[road,2],pts[road[1],2]),'o')
  }
)

'<-'(
  compare,
  function(pn=30)
  {
    px <- runif(pn,100,1000)
    py <- runif(pn,100,1000)
    pts <- matrix(c(px,py),pn,2)
    ra <- AcoTsp(dist(pts))
    rs <- SaTsp(dist(pts))
    split.screen(c(2,2))
    screen(1)
    plot(pts[,1],pts[,2],xlab='x',ylab='y',main='Original')
    screen(2)
    plot(pts[,1],pts[,2],'o',xlab='x',ylab='y',main='Random')
    screen(3)
    plot(c(pts[ra,1], pts[ra[1],1]),c(pts[ra,2],pts[ra[1],2]),'o',xlab='x',ylab='y',main='ACO')
    screen(4)
    plot(c(pts[rs,1], pts[rs[1],1]),c(pts[rs,2],pts[rs[1],2]),'o',xlab='x',ylab='y',main='SA')
    close.screen(all=TRUE)
  }
)