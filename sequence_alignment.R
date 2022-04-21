#Author Liang
#time 2022-4-1
#
#
#
require('Rcpp')
sourceCpp('./RStudy/sequence_alignment.cpp')
liang_NW(as.matrix(c('c','g')), as.matrix(c('a','t','t','c','g')), print)
temp.no <- liang_NW(dna1, dna2, print)
require('parallel')
'<-'(
  liang.dsaNW,
  function(seq1, seq2, markFunc=NULL)
  {
    '<-'(
      getMark,
      function(b1, b2)
      {
        #默认采用匹配计分
        bb1 <- b1=='-'
        bb2 <- b2=='-'
        return(ifelse(bb1||bb2, ifelse(bb1&&bb2, 0, -5), ifelse(b1==b2, 5, -4)))
      }
    )
    markFunc <- ifelse(is.null(markFunc), getMark, markFunc)
    '<-'(
      getMarks,
      function(b, bs)
      {
        len <- length(bs)
        mark <- numeric(len)
        mark[1] <- markFunc(b, bs[1])
        '<-'(
          eachIter,
          function(i)
          {
            return(mark[i-1] + markFunc(b, bs[i]))
          }
        )
        mark[-1] <- sapply(2:len, eachIter)
        return(mark)
      }
    )
    
    #转化为方便处理的形式
    seq1 <- unlist(seq1)
    seq2 <- unlist(seq2)
    #获取长度，方便后续处理
    len1 <- length(seq1)
    len2 <- length(seq2)
    delSeq <- function(s){return(as.character(unlist(strsplit(s,split=''))))}
    if(len1==1)
    {
      seq1 <- delSeq(seq1)
      len1 <- length(seq1)
    }
    if(len2==1)
    {
      seq2 <- delSeq(seq2)
      len2 <- length(seq2)
    }
    #保证seq1长度最小，方便后面优化
    isChange <- len1>len2
    if(isChange)
    {
      tmp <- seq1
      seq1 <- seq2
      seq2 <- tmp
      tmp <- len1
      len1 <- len2
      len2 <- tmp
    }
    
    #初始化
    marks <- matrix(0, len1+1, len2+1) #计分矩阵
    marks[1,-1] <- getMarks('-', seq2)
    marks[-1,1] <- getMarks('-', seq1)
    roads <- matrix(0, len1+1, len2+1) #匹配的结果，最差情况是完全错位
    roads[1,-1] <- 2
    roads[-1,1] <- 3
    
    #写操作
    '<-'(
      writeMatrix,
      function(sp1, sp2, mark)
      {
        pl1 <- length(sp1)
        pl2 <- length(sp2)
        road <- matrix(0, pl1, pl2)
        '<-'(
          rowIter,
          function(i)
          {
            '<-'(
              colIter,
              function(k)
              {
                tmark <- numeric(3)
                tmark[1] <- mark[i-1,k-1] + markFunc(sp1[i-1], sp2[k-1]) #markLerp
                tmark[2] <- mark[i,k-1] + markFunc('-', sp2[k-1]) #markLeft
                tmark[3] <- mark[i-1,k] + markFunc(sp1[i-1], '-') #markUp
                way <- which.max(tmark)
                road[i-1,k-1] <<- way
                return(tmark[way])
              }
            )
            return(sapply(2:(1+pl2),colIter))
          }
        )
        mark[-1,-1] <- sapply(2:(1+pl1),rowIter)
        tmp <- list(mark, road)
        names(tmp) <- c('mark', 'road')
        return(tmp)
      }
    )
    
    if(len1 > 800 && len2 > 800)
    {
      dx <- floor(len2 / 2)
      dy <- floor(len1 / 2)
      clist <- makeCluster(2)
      pans <- writeMatrix(seq1[1:dy],seq2[1:dx], marks[1:(1+dy),1:(1+dx)])
      marks[1:(1+dx),1:(1+dy)] <- pans$mark
      roads[-1,-1][1:dy,1:dx] <- pans$road

      '<-'(
        divideDo,
        function(w)
        {
          if(w==1)
          {
            pans <- writeMatrix(seq1[1:dy],seq2[(dx+1):len2], marks[1:(1+dy),(dx+1):(len2+1)])
          }
          else
          {
            pans <- writeMatrix(seq1[(dy+1):len1],seq2[1:dx], marks[(dy+1):(len1+1),1:(1+dx)])
          }
          return(pans)
        }
      )
      try(tpans <- parLapply(clist, 1:2, divideDo))
      stopCluster(clist)
      marks[1:(dy+1),(dx+1):(len2+1)] <- tpans[[1]]$mark
      roads[-1,-1][1:dy, (dx+1):len2] <- tpans[[1]]$road
      marks[(dy+1):(len1+1),1:(1+dx)] <- tpans[[2]]$mark
      roads[-1,-1][(dy+1):len1, 1:dx] <- tpans[[2]]$road
      
      pans <- writeMatrix(seq1[(dy+1):len1],seq2[(dx+1):len2], marks[(dy+1):(len1+1),(dx+1):(len2+1)])
      marks[(dy+1):(len1+1),(dx+1):(len2+1)] <- pans$mark
      roads[-1,-1][(dy+1):len1,(dx+1):len2] <- pans$road
    }
    else
    {
      pans <- writeMatrix(seq1, seq2, marks)
      marks <- pans$mark
      roads[-1,-1] <- pans$road 
    }
    
    #回溯
    ans <- matrix('-', len1+len2, 2)
    nowi <- len1 + 1
    nowj <- len2 + 1
    long <- len1+len2

    while(roads[nowi,nowj] > 0)
    {
      sans <- ifelse(roads[nowi,nowj]==1,list(c(seq1[nowi-1],seq2[nowj-1],nowi-1,nowj-1)),ifelse(roads[nowi,nowj]==2, list(c('-',seq2[nowj-1],nowi,nowj-1)),list(c(seq1[nowi-1],'-',nowi-1,nowj))))
      sans <- unlist(sans)
      ans[long, 1] <- sans[1]
      ans[long, 2] <- sans[2]
      nowi <- as.numeric(sans[3])
      nowj <- as.numeric(sans[4])
      long <- long - 1
    }
    #保持和输入顺序一致
    if(isChange)
    {
      ans <- ans[,c(2,1)]
    }
    ans <- ans[-(1:long),]
    return(ans)
  }
)

'<-'(
  liang.dsaAC,
  function(seq1, seq2)
  {
    #转化为方便处理的形式
    seq1 <- unlist(seq1)
    seq2 <- unlist(seq2)
    #获取长度，方便后续处理
    len1 <- length(seq1)
    len2 <- length(seq2)
    delSeq <- function(s){return(as.character(unlist(strsplit(s,split=''))))}
    if(len1==1)
    {
      seq1 <- delSeq(seq1)
      len1 <- length(seq1)
    }
    if(len2==1)
    {
      seq2 <- delSeq(seq2)
      len2 <- length(seq2)
    }
    #保证seq1长度最小，方便后面优化
    isChange <- len1>len2
    if(isChange)
    {
      tmp <- seq1
      seq1 <- seq2
      seq2 <- tmp
      tmp <- len1
      len1 <- len2
      len2 <- tmp
    }
    alen <- len1 + len2 #总长度
    mlen <- (len1+1) * (len2+1) #矩阵元素个数
    
    #参数设置
    #方向权值,1下2右3斜
    dway <- c(1,2,3)
    #防止偏离参数
    H <- 5
    #信息挥发 1-p
    pho <- 0.95
    #随机阈值
    p1 <- 0.6
    p2 <- 0.35
    p3 <- 0.2
    p0 <- p1
    #信息素启发参数，分值启发，偏离启发
    alpha <- -5
    beta <- -3
    cita <- 2
    #迭代次数
    nIter <- 500
    #蚂蚁数量
    ants <- 100
    
    #信息素矩阵
    A <- array(5, dim=c(len1+1,len2+1,3))
    #信息素更新参数
    Q1 <- 0.1
    Q2 <- 0.2
    Q3 <- -50
    #匹配矩阵
    D <- array(0.5, dim=c(len1+1,len2+1,3))
    D[-1,-1,3] <- 1
    D[-(len1+1),-(len2+1),3][t(sapply(seq1, function(b){return(seq2==b)}))] <- 2
    
    gRoad <- numeric(alen)
    gMark <- -Inf
    hsy <- numeric(nIter)
    
    '<-'(
      chooseWay,
      function(i,j,dk)
      {
        if(runif(1) > p0)
        {
          return(3)
        }
        pway <- Tmp[i,j,] + dk^cita
        pway <- pway / sum(pway)
        choose <- sample(1:3, 1, prob = pway)
        return(choose)
      }            
    )
    
    '<-'(
      eachAnt,
      function(who)
      {
        i <- 1
        j <- 1
        long <- 0
        road <- numeric(alen)
        dk <- dway
        mark <- 0
        atmp <- array(0, dim=c(len1+1,len2+1,3)) #记录信息素更新位置
        while(TRUE)
        {
          isLow <- i==(len1+1)
          isRight <- j==(len2+1)
          if(isLow && isRight)
          {
            break
          }
          long <- long + 1
          h <- i - j # 偏离量
          dk[1] <- ifelse(-h > H, dk[1]*0.7, dway[1]) #如果太下，减少向下的权值
          dk[2] <- ifelse(h > H, dk[2]*0.7, dway[2]) #如果太右，减少向右的权值
          way <- ifelse(isRight, 1, ifelse(isLow, 2, chooseWay(i,j,dk)))
          road[long] <- way
          mark <- mark + ifelse(D[i,j,way]!=2, ifelse(D[i,j,way]==1,-1,2), -2)
          #atmp[i,j,way] <- ifelse(D[i,j,way]!=2, ifelse(D[i,j,way]==1,0,2), 0.5)
          atmp[i,j,way] <- 1
          i <- ifelse(way==2,i,i+1)
          j <- ifelse(way==1,j,j+1)
        }
        ans <- list(road[1:long], atmp, mark)
        names(ans) <- c('r','a','mark')
        return(ans)
      }
    )
    for(I in 1:nIter)
    {
      p0 <- ifelse(I < 0.3, p1, ifelse(I < 0.7, p2, p3))
      Tmp <- A^alpha + D^beta #计算转移概率的中间矩阵
      ans <- lapply(1:ants, eachAnt)
      maxMark <- -Inf
      wmax <- 1
      A <- A * pho
      for(i in 1:ants)
      {
        each <- ans[[i]]
        mark <- each$mark
        maxMark <- ifelse(mark>maxMark, mark, maxMark)
        wmax <- ifelse(mark>maxMark, i, wmax)
        A <- A + Q1 * each$a * ifelse(mark>0,mark,1+(mark/Q3))
      }
      #A <- A + (Q2-Q1) * ans[[wmax]]$a #加强
      if(maxMark > gMark)
      {
        gMark <- maxMark
        gRoad <- ans[[wmax]]$r
      }
      hsy[I] <- maxMark
    }
    
    lenRoad <- length(gRoad)
    final <- matrix('-', lenRoad, 2)
    nowi <- 1
    nowj <- 1
    long <- 1
    for(k in 1:lenRoad)
    {
      if(gRoad[k]==1)
      {
        final[long,1] <- seq1[nowi]
        #final[long,2] <- '-'
        nowi <- nowi + 1
      }
      else if(gRoad[k]==2)
      {
        #final[long,1] <- '-'
        final[long,2] <- seq2[nowj]
        nowj <- nowj + 1
      }
      else
      {
        final[long,1] <- seq1[nowi]
        final[long,2] <- seq2[nowj]
        nowi <- nowi + 1
        nowj <- nowj + 1
      }
      long <- long + 1
    }
    if(isChange)
    {
      final <- final[,c(2,1)]
    }
    
    plot(1:nIter, hsy, 'l',xlab='times', ylab='mark', main='times-mark')
    return(final)
  }
)

'<-'(
  liang.dsaSA,
  function(seq1, seq2, TS = 100,TE = 20)
  {
    #转化为方便处理的形式
    seq1 <- unlist(seq1)
    seq2 <- unlist(seq2)
    #获取长度，方便后续处理
    len1 <- length(seq1)
    len2 <- length(seq2)
    delSeq <- function(s){return(as.character(unlist(strsplit(s,split=''))))}
    if(len1==1)
    {
      seq1 <- delSeq(seq1)
      len1 <- length(seq1)
    }
    if(len2==1)
    {
      seq2 <- delSeq(seq2)
      len2 <- length(seq2)
    }
    #保证seq1长度最小，方便后面优化
    isChange <- len1>len2
    if(isChange)
    {
      tmp <- seq1
      seq1 <- seq2
      seq2 <- tmp
      tmp <- len1
      len1 <- len2
      len2 <- tmp
    }
    alen <- len1 + len2 #总长度
    
    '<-'(
      interfere,
      function(ord)
      {
        lord <- length(ord)
        pos0 <- which(ord==0)
        if(length(pos0) && runif(1) > 0.7)
        {
          nord <- numeric(lord-1)
          div <- ifelse(length(pos0)==1, pos0, sample(pos0,1)) #请注意这个sample函数
          if(div==1)
          {
            nord <- ord[2:lord]
          }
          else if(div==lord)
          {
            nord <- ord[1:(lord-1)]
          }
          else
          {
            nord[1:(div-1)] <- ord[1:(div-1)]
            nord[div:(lord-1)] <- ord[(div+1):lord]
          }
        }
        else
        {
          nord <- numeric(lord+1)
          div <- round(runif(1, 0, lord+1))
          nord <- append(ord, 0, div)
        }
        return(nord)
      }
    )
    
    '<-'(
      delOrder,
      function(ord1, ord2)
      {
        lord1 <- length(ord1)
        lord2 <- length(ord2)
        ml <- max(lord1, lord2)
        tm <- matrix(0, 2, ml)
        tm[1,1:lord1] <- ord1
        tm[2,1:lord2] <- ord2
        tm[1,tm[1,]==1] <- seq1
        tm[2,tm[2,]=='1'] <- seq2
        tm[which(tm=='0')] <- '-'
        #不允许空格对空格
        sp <- which(tm[1,]=='-' & tm[2,]=='-')
        if(length(sp))
        {
          tm <- tm[1:2,-sp]
        }
        return(tm)
      }
    )
    
    '<-'(
      calMark,
      function(pm)
      {
        m1 <- length(c(which(pm[1,]=='-'), which(pm[2,]=='-'))) #空格对字母
        m2 <- length(which(pm[1,]==pm[2,])) #字符对字符，而且一致
        m3 <- dim(pm)[2] - m1 - m2 #字符对字符，但不一致
        return(m1 * -2 + m2 * 2 + m3 * -1)
      }
    )
    
    '<-'(
      mInterfere,
      function(ord, times=10)
      {
        for(i in 1:times)
        {
          ord <- interfere(ord)
        }
        return(ord)
      }
    )
    
    nT <- TS
    order1 <- rep(1, len1)
    order2 <- rep(1, len2)
    bestMark <- calMark(delOrder(order1, order2))
    alpha <- 0.99
    cout <- 1
    record <- numeric(floor(log(TE/TS, alpha))+10)
    record[cout] <- bestMark
    hsy <- bestMark
    hsy1 <- order1
    hsy2 <- order2
    while(nT > TE)
    {
      cnt <- 10
      cout <- cout + 1
      while(cnt > 0)
      {
        cnt <- cnt - 1
        to1 <- interfere(order1)
        to2 <- interfere(order2)
        tm <- delOrder(to1, to2)
        tmark <- calMark(tm)
        #和TSP不同的是，我们此处的评价标准变成了分数，分数可以是负数和0，所以deta要相应调整
        deta <- ifelse(bestMark, (bestMark - tmark) / abs(bestMark), 0)
        if(deta <= 0 || exp(-nT/deta) > runif(1))
        {
          bestMark <- tmark
          order1 <- as.numeric(tm[1,]!='-')
          order2 <- as.numeric(tm[2,]!='-')
          cnt <- ifelse(deta==0,cnt,10)
        }
      }
      nT <- nT * alpha
      record[cout] <- bestMark
      if(bestMark >= hsy)
      {
        hsy <- bestMark
        hsy1 <- order1
        hsy2 <- order2
      }
    }
    
    #plot(1:cout, record[1:cout], 'l', xlab='each_Iter', ylab='mark', main='mark-iter_times')
    #print(hsy)
    return(t(delOrder(hsy1, hsy2)))
  }
)

'<-'(
  liang.msaGA,
  function(seqs, gnum = 1000, pv = 0.01)
  {
    snum <- length(seqs)
    maxLen <- ceiling(max(sapply(seqs, function(s){return(length(s))})) * 1.2)
    '<-'(
      randomStart,
      function(id)
      {
        for(i in 1:snum)
        {
          tmp <- rep('-', maxLen)
          tmp[sort(sample(1:maxLen, length(seqs[[i]])))] <- seqs[[i]]
          seqs[[i]] <- tmp
        }
        return(seqs)
      }
    )
    '<-'(
      delAlign,
      function(seqs)
      {
        #删除全空格行
        #生成对齐矩阵
        tl <- sapply(seqs, function(s){return(length(s))})
        ml <- max(tl)
        tm <- matrix('-', ml, snum)
        for(i in 1:snum)
        {
          tm[1:tl[i],i] <- seqs[[i]]
        }
        sp <- which(unlist(apply(tm,1,function(l){return(all(l=='-'))})))
        if(length(sp))
        {
          tm <- tm[-sp,]
        }
        return(tm)
      }
    )
    '<-'(
      SP,
      function(ind)
      {
        tm <- delAlign(ind)
        di <- dim(tm)
        csum <- di[2]
        mark <- 0.001
        for(i in 1:(csum-1))
        {
          b1 <- tm[,(i+1):(csum-1)] == tm[,i] #两个位置一致
          b2 <- tm[,(i+1):(csum-1)] == '-' #被比对位置为空格
          m1 <- length(which(!b1 & b2)) #空格对字母
          m2 <- length(which(b1 & !b2)) #字符对字符，而且一致
          #m3 <- length(which(b1 & b2)) #空格对空格
          m4 <- length(which(!b1 & !b2)) #字符对字符，但不 一致
          mark <- mark + m1 * -2 + m2 * 2 + m4 * 1
        }
        mark <- ifelse(mark==0,1,ifelse(mark<0,-1/mark,mark))
        return(mark)
      }
    )
    '<-'(
      cross,
      function(id, ind1, ind2)
      {
        #仅选择多行横向交叉一种算子
        rd <- round(runif(1, 2, snum))
        p1 <- ind1[rd:snum]
        p2 <- ind2[rd:snum]
        ind1[rd:snum] <- p2
        ind2[rd:snum] <- p1
        #为加快运行速度
        group[[gnum+id]] <<- ind1
        group[[gnum+id+1]] <<- ind2
      }
    )
    '<-'(
      vary,
      function(ind)
      {
        if(runif(1) > pv)
        {
          return(ind)
        }
        vp <- round(runif(1, 1, snum))
        vl <- length(ind[[vp]])
        tmp <- rep('-', vl)
        tmp[sort(sample(1:vl, length(seqs[[vp]])))] <- seqs[[vp]]
        ind[[vp]] <- tmp
        return(ind)
      }
    )
    #产生初始种群
    group <- lapply(1:(2*gnum), function(x){return(NA)})
    group[1:gnum] <- lapply(1:gnum, randomStart)
    imark <- sapply(group[1:gnum], SP)
    cnt <- 0
    record <- numeric(100000)
    best <- NA
    bestMark <- 0
    iter <- 200
    while(iter > 0)
    {
      print(cnt)
      cnt <- cnt + 1
      #交叉，优秀的个体繁衍的机会越大
      for(i in seq(1, gnum, 2))
      {
        tg <- sample(1:gnum, 2, prob = imark)
        cross(i, group[[tg[1]]], group[[tg[2]]])
      }
      #变异
      group <- lapply(group, vary)
      #计算适应度
      imark <- sapply(group, SP)
      #一些记录
      record[cnt] <- max(imark)
      if(record[cnt] > bestMark)
      {
        bestMark <- record[cnt]
        best <- group[[which.max(imark)]]
        iter <- 200
      }
      #轮盘赌选择
      chInd <- sample(1:(2*gnum), gnum, prob = imark)
      #精英选择
      #chInd <- which(order(imark)<=gnum)
      group[1:gnum] <- group[chInd]
      imark <- imark[chInd]

      iter <- iter - 1
    }
    print(bestMark)
    
    plot(1:cnt, record[1:cnt], 'l', col=2, xlab='iter_time',ylab='sp_mark')
    return(delAlign(best))
  }
)

'<-'(
  liang.GASA,
  function(seqs, gnum = 30, pv = 0.01)
  {
    snum <- length(seqs)
    maxLen <- ceiling(max(sapply(seqs, function(s){return(length(s))})) * 1.2)
    '<-'(
      randomStart,
      function(id)
      {
        for(i in 1:snum)
        {
          tmp <- rep('-', maxLen)
          tmp[sort(sample(1:maxLen, length(seqs[[i]])))] <- seqs[[i]]
          seqs[[i]] <- tmp
        }
        return(seqs)
      }
    )
    '<-'(
      delAlign,
      function(seqs)
      {
        #删除全空格行
        #生成对齐矩阵
        tl <- sapply(seqs, function(s){return(length(s))})
        ml <- max(tl)
        tm <- matrix('-', ml, snum)
        for(i in 1:snum)
        {
          tm[1:tl[i],i] <- seqs[[i]]
        }
        sp <- which(unlist(apply(tm,1,function(l){return(all(l=='-'))})))
        if(length(sp))
        {
          tm <- tm[-sp,]
        }
        return(tm)
      }
    )
    '<-'(
      SP,
      function(ind)
      {
        tm <- delAlign(ind)
        di <- dim(tm)
        csum <- di[2]
        mark <- 0.001
        for(i in 1:(csum-1))
        {
          b1 <- tm[,(i+1):(csum-1)] == tm[,i] #两个位置一致
          b2 <- tm[,(i+1):(csum-1)] == '-' #被比对位置为空格
          m1 <- length(which(!b1 & b2)) #空格对字母
          m2 <- length(which(b1 & !b2)) #字符对字符，而且一致
          #m3 <- length(which(b1 & b2)) #空格对空格
          m4 <- length(which(!b1 & !b2)) #字符对字符，但不 一致
          mark <- mark + m1 * -2 + m2 * 1 + m4 * -1
        }
        mark <- ifelse(mark==0,1,ifelse(mark<0,-1/mark,mark))
        return(mark)
      }
    )
    '<-'(
      cross,
      function(id, ind1, ind2)
      {
        #仅选择多行横向交叉一种算子
        rd <- round(runif(1, 2, snum))
        p1 <- ind1[rd:snum]
        p2 <- ind2[rd:snum]
        ind1[rd:snum] <- p2
        ind2[rd:snum] <- p1
        #为加快运行速度
        group[[gnum+id]] <<- ind1
        group[[gnum+id+1]] <<- ind2
      }
    )
    '<-'(
      interfere,
      function(ind)
      {
        rd <- runif(1)
        for(vp in 1:snum)
        {
          tmp <- ind[[vp]]
          tl <- length(tmp)
          posk <- which(tmp=='-')
          if(rd>0.7 && length(posk))
          {
            posk <- ifelse(length(posk)==1, posk, sample(posk,1))
            if(posk==1)
            {
              tmp <- tmp[2:tl]
            }
            else if(posk==tl)
            {
              tmp <- tmp[1:(tl-1)]
            }
            else
            {
              ttmp <- character(tl-1)
              ttmp[1:(posk-1)] <- tmp[1:(posk-1)]
              ttmp[posk:(tl-1)] <- tmp[(posk+1):tl]
              tmp <- ttmp
            }
          }
          else if(rd > 0.4 && length(posk))
          {
            posk <- ifelse(length(posk)==1, posk, sample(posk,1))
            posz <- sample(which(tmp!='-'),1)
            if(posk==1)
            {
              tmp[1:(tl-1)] <- tmp[2:tl]
              tmp[tl] <- '-'
            }
            else if(posk==tl)
            {
              tmp[2:tl] <- tmp[1:(tl-1)]
              tmp[1] <- '-'
            }
            else
            {
              ttmp <- character(tl-1)
              ttmp[1:(posk-1)] <- tmp[1:(posk-1)]
              ttmp[posk:(tl-1)] <- tmp[(posk+1):tl]
              tmp <- append(ttmp, '-', posz-1)
            }
          }
          else
          {
            tmp <- append(tmp, '-', round(runif(1, 0, tl)))
          }
          ind[[vp]] <- tmp
        }
        return(ind)
      }
    )
    '<-'(
      vary,
      function(ind)
      {
        if(runif(1) > pv)
        {
          return(ind)
        }
        vp <- round(runif(1, 1, snum))
        vl <- length(ind[[vp]])
        tmp <- rep('-', vl)
        tmp[sort(sample(1:vl, length(seqs[[vp]])))] <- seqs[[vp]]
        ind[[vp]] <- tmp
        return(ind)
      }
    )
    #产生初始种群
    group <- lapply(1:(2*gnum), function(x){return(NA)})
    group[1:gnum] <- lapply(1:gnum, randomStart)
    imark <- sapply(group[1:gnum], SP)
    cnt <- 0
    record <- numeric(100000)
    best <- NA
    bestMark <- 0
    iter <- 100
    while(iter > 0)
    {
      cnt <- cnt + 1
      #交叉，优秀的个体繁衍的机会越大
      for(i in seq(1, gnum, 2))
      {
        tg <- sample(1:gnum, 2, prob = imark)
        cross(i, group[[tg[1]]], group[[tg[2]]])
      }
      #变异
      group <- lapply(group, vary)
      #计算适应度
      imark <- sapply(group, SP)
      #一些记录
      record[cnt] <- max(imark)
      if(record[cnt] > bestMark)
      {
        bestMark <- record[cnt]
        best <- group[[which.max(imark)]]
        iter <- 100
      }
      #轮盘赌选择进入下一代的个体
      chInd <- sample(1:(2*gnum), gnum, prob = imark)
      group[1:gnum] <- group[chInd]
      imark <- imark[chInd]
      
      iter <- iter - 1
    }
    dcnt <- cnt
    
    #遗传算法局部寻优能力不强，接下来用模拟退火进行局部寻优
    Tp <- 3000
    alpha <- 0.9
    oind <- best
    omark <- bestMark
    while(Tp > 20)
    {
      cnt <- cnt + 1
      tcnt <- 10
      while(tcnt > 0)
      {
        tind <- interfere(oind)
        tmark <- SP(tind)
        deta <- (omark - tmark) / omark
        if(deta <= 0 || exp(-Tp/deta) > runif(1))
        {
          oind <- tind
          omark <- tmark
          if(deta < 0)
          {
            tcnt <- 10
            if(bestMark < tmark)
            {
              bestMark <- tmark
              best <- tind
            }
          }
        }
        tcnt <- tcnt - 1
      }
      record[cnt] <- omark
      Tp <- Tp * alpha
    }
    print(bestMark)
    
    plot((dcnt+1):cnt, record[(dcnt+1):cnt], 'l', col=2, xlab='iter_time',ylab='sp_mark', xlim=c(1,cnt), ylim=c(min(record[1:cnt]),bestMark))
    lines(1:dcnt, record[1:dcnt], 'l', col=3)
    return(delAlign(best))
  }
)


'<-'(
  liang.msaSA,
  function(seqs, TS=3000, TE=100, tim = 10)
  {
    snum <- length(seqs)
    maxLen <- ceiling(max(sapply(seqs, function(s){return(length(s))})) * 1.2)
    '<-'(
      randomStart,
      function(id)
      {
        for(i in 1:snum)
        {
          tmp <- rep('-', maxLen)
          tmp[sort(sample(1:maxLen, length(seqs[[i]])))] <- seqs[[i]]
          seqs[[i]] <- tmp
        }
        return(seqs)
      }
    )
    '<-'(
      delAlign,
      function(seqs)
      {
        #删除全空格行
        #生成对齐矩阵
        tl <- sapply(seqs, function(s){return(length(s))})
        ml <- max(tl)
        tm <- matrix('-', ml, snum)
        for(i in 1:snum)
        {
          tm[1:tl[i],i] <- seqs[[i]]
        }
        sp <- which(unlist(apply(tm,1,function(l){return(all(l=='-'))})))
        if(length(sp))
        {
          tm <- tm[-sp,]
        }
        return(tm)
      }
    )
    '<-'(
      SP,
      function(ind)
      {
        tm <- delAlign(ind)
        di <- dim(tm)
        csum <- di[2]
        mark <- 0.001
        for(i in 1:(csum-1))
        {
          b1 <- tm[,(i+1):(csum-1)] == tm[,i] #两个位置一致
          b2 <- tm[,(i+1):(csum-1)] == '-' #被比对位置为空格
          m1 <- length(which(!b1 & b2)) #空格对字母
          m2 <- length(which(b1 & !b2)) #字符对字符，而且一致
          #m3 <- length(which(b1 & b2)) #空格对空格
          m4 <- length(which(!b1 & !b2)) #字符对字符，但不 一致
          mark <- mark + m1 * -5 + m2 * 5 + m4 * -4
        }
        mark <- ifelse(mark==0,1,ifelse(mark<0,-1/mark,mark))
        return(mark)
      }
    )
    '<-'(
      interfere,
      function(ind)
      {
        rd <- runif(1)
        for(vp in 1:snum)
        {
          tmp <- ind[[vp]]
          tl <- length(tmp)
          posk <- which(tmp=='-')
          if(rd>0.8 && length(posk))
          {
            posk <- ifelse(length(posk)==1, posk, sample(posk,1))
            if(posk==1)
            {
              tmp <- tmp[2:tl]
            }
            else if(posk==tl)
            {
              tmp <- tmp[1:(tl-1)]
            }
            else
            {
              ttmp <- character(tl-1)
              ttmp[1:(posk-1)] <- tmp[1:(posk-1)]
              ttmp[posk:(tl-1)] <- tmp[(posk+1):tl]
              tmp <- ttmp
            }
          }
          else if(rd > 0.5 && length(posk))
          {
            posk <- ifelse(length(posk)==1, posk, sample(posk,1))
            posz <- sample(which(tmp!='-'),1)
            if(posk==1)
            {
              tmp[1:(tl-1)] <- tmp[2:tl]
              tmp[tl] <- '-'
            }
            else if(posk==tl)
            {
              tmp[2:tl] <- tmp[1:(tl-1)]
              tmp[1] <- '-'
            }
            else
            {
              ttmp <- character(tl-1)
              ttmp[1:(posk-1)] <- tmp[1:(posk-1)]
              ttmp[posk:(tl-1)] <- tmp[(posk+1):tl]
              tmp <- append(ttmp, '-', posz-1)
            }
          }
          else
          {
            tmp <- append(tmp, '-', round(runif(1, 0, tl)))
          }
          ind[[vp]] <- tmp
        }
        return(ind)
      }
    )
    Tp <- TS
    alpha <- 0.99
    oind <- randomStart()
    omark <- SP(oind)
    cnt <- 0
    record <- numeric(floor(log(TE/TS, alpha))+10)
    bestMark <- 0
    best <- NA
    while(Tp > TE)
    {
      cnt <- cnt + 1
      tcnt <- tim
      while(tcnt > 0)
      {
        tind <- interfere(oind)
        tmark <- SP(tind)
        deta <- (omark - tmark) / omark
        if(deta <= 0 || exp(-deta/Tp) > runif(1))
        {
          oind <- tind
          omark <- tmark
          if(deta < 0)
          {
            tcnt <- tim
            if(bestMark < tmark)
            {
              bestMark <- tmark
              best <- tind
            }
          }
        }
        tcnt <- tcnt - 1
      }
      record[cnt] <- omark
      Tp <- Tp * alpha
    }
    print(bestMark)
    plot(1:cnt,record[1:cnt],'l',xlab='iter_time',ylab='sp_mark')
    return(delAlign(best))
  }
)

'<-'(
  liang.msaStar,
  function(seqs)
  {
    '<-'(
      SP,
      function(tm)
      {
        di <- dim(tm)
        csum <- di[2]
        mark <- 0
        for(i in 1:(csum-1))
        {
          b1 <- tm[,(i+1):(csum-1)] == tm[,i] #两个位置一致
          b2 <- tm[,(i+1):(csum-1)] == '-' #被比对位置为空格
          m1 <- length(which(!b1 & b2)) #空格对字母
          m2 <- length(which(b1 & !b2)) #字符对字符，而且一致
          #m3 <- length(which(b1 & b2)) #空格对空格
          m4 <- length(which(!b1 & !b2)) #字符对字符，但不 一致
          mark <- mark + m1 * -5 + m2 * 5 + m4 * -4
        }
        mark <- ifelse(mark==0,1,ifelse(mark<0,-1/mark,mark))
        return(mark)
      }
    )
    
    num <- length(seqs)
    marks <- matrix(0, num, num, dimnames = list(names(seqs), names(seqs)))
    seqs <- lapply(seqs, as.matrix)
    for(i in 1:num){
      for(k in 1:num){
        if(i>=k){
          next
        }
        tmp <- liang_NW(seqs[[i]], seqs[[k]], print)
        marks[i,k] <- SP(tmp) 
      }
    }
    marks <- as.matrix(as.dist(marks))
    now <- which.max(rowSums(marks))
    had <- numeric(num)
    ans <- seqs[[now]]
    for(i in 1:(num-1)){
      had[i] <- now
      nxt <- -Inf
      odr <- order(marks[now,])
      for(i in 1:num)
      {
        if(i %in% had)
        {
          next
        }
        nxt <- ifelse(odr[i]>=nxt,i,nxt)
      }
      ans <- liang_NW(seqs[[nxt]], ans, print)
      now <- nxt
    }
    #print(had)
    #print(SP(ans))
    return(ans)
  }
)

'<-'(
  liang.msaStarClimb,
  function(seqs, bear=1000, delSpace=TRUE, hb = NULL)
  {
    '<-'(
      SP,
      function(tm)
      {
        di <- dim(tm)
        csum <- di[2]
        mark <- 0.1
        for(i in 1:(csum-1))
        {
          b1 <- tm[,(i+1):(csum-1)] == tm[,i] #两个位置一致
          b2 <- tm[,(i+1):(csum-1)] == '-' #被比对位置为空格
          m1 <- length(which(!b1 & b2)) #空格对字母
          m2 <- length(which(b1 & !b2)) #字符对字符，而且一致
          #m3 <- length(which(b1 & b2)) #空格对空格
          m4 <- length(which(!b1 & !b2)) #字符对字符，但不 一致
          mark <- mark + m1 * -5 + m2 * 5 + m4 * -4
        }
        mark <- ifelse(mark==0,1,ifelse(mark<0,-1/mark,mark))
        return(mark)
      }
    )
    '<-'(
      Interfere,
      function(tm)
      {
        di <- dim(tm)
        ch <- round(runif(1, 1, di[2]))
        v <- tm[,ch]
        iszero <- which(v=='-')
        if(length(iszero)<=1)
        {
          if(length(iszero))
          {
            p <- iszero
          }
          else
          {
            return(tm)
          }
        }
        else
        {
          p <- sample(iszero, 1)
        }
        if(p==1)
        {
          tmp <- v[2:di[1]]
        }
        else if(p == di[1])
        {
          tmp <- v[1:(di[1]-1)]
        }
        else
        {
          tmp <- character(di[1]-1)
          tmp[1:(p-1)] <- v[1:(p-1)]
          tmp[p:(di[1]-1)] <- v[(p+1):di[1]]
        }
        s <- ifelse(p<5,0,p-5)
        e <- ifelse(p>di[1]-6,di[1]-1,p+5)
        v <- append(tmp, '-', sample(s:e, 1))
        tm[,ch] <- v
        return(tm)
      }
    )
    
    
    yuzhi <- bear
    if(is.null(hb))
    {
      ans <- liang.msaStar(seqs) 
    }
    else
    {
      ans <- hb
    }
    mark <- SP(ans)
    rnum <- dim(ans)[1]
    
    #record <- numeric(1000000)
    #tim <- -1
    while(yuzhi > 0)
    {
      #tim <- tim + 1
      tans <- Interfere(ans)
      tmark <- SP(tans)
      if(mark > tmark)
      {
        yuzhi <- yuzhi - 1 
      }
      else
      {
        yuzhi <- ifelse(tmark>mark, bear, yuzhi)
        ans <- tans
        mark <- tmark
      }
      #record[tim] <- mark
    }
    #print(mark)
    #plot(1:tim, record[1:tim], 'l',xlab='iter', ylab='mark')
    if(delSpace)
    {
      weishu <- dim(ans)[1]
      weishu <- which(sapply(1:weishu, function(x){return(all(ans[x,]=='-'))}))
      if(length(weishu))
      {
        ans <- ans[-weishu,]
      } 
    }
    return(ans)
  }
)

'<-'(
  liang.msaStarClimbGA,
  function(seqs, gnum=100, vp=0.1)
  {
    '<-'(
      SP,
      function(tm)
      {
        di <- dim(tm)
        csum <- di[2]
        mark <- 0
        for(i in 1:(csum-1))
        {
          b1 <- tm[,(i+1):(csum-1)] == tm[,i] #两个位置一致
          b2 <- tm[,(i+1):(csum-1)] == '-' #被比对位置为空格
          m1 <- length(which(!b1 & b2)) #空格对字母
          m2 <- length(which(b1 & !b2)) #字符对字符，而且一致
          #m3 <- length(which(b1 & b2)) #空格对空格
          m4 <- length(which(!b1 & !b2)) #字符对字符，但不 一致
          mark <- mark + m1 * -5 + m2 * 5 + m4 * -4
        }
        mark <- ifelse(mark==0,1,ifelse(mark<0,-1/mark,mark))
        return(mark)
      }
    )
    hb <- liang.msaStar(seqs)
    di <- dim(hb)
    print('序列初步比对完成')
    clist <- makeCluster(detectCores(logical = F))
    sapply(ls(".GlobalEnv"), (function(t){return(clusterExport(clist, t, envir = environment()))}))
    try(group <- parLapply(clist, 1:gnum, function(x){return(liang.msaStarClimb(seqs, 100, FALSE, hb=hb))}))
    print('初始种群建立完成')
    
    '<-'(
      getChild,
      function(ind1, ind2)
      {
        s <- round(runif(1, 1, di[2]))
        e <- round(runif(1, s, di[2]))
        ind <- NA
        if(s==1)
        {
          if(e==di[[2]])
          {
            ind <- ind2
          }
          else
          {
            ind <- cbind(ind2[,1:e],ind1[,(e+1):di[2]]) 
          }
        }
        else if(e==di[2])
        {
          ind <- cbind(ind1[,1:(s-1)],ind2[,s:e])
        }
        else
        {
          ind <- cbind(ind1[,1:(s-1)],ind2[,s:e],ind1[,(e+1):di[2]])
        }
        return(ind)
      }
    )
    
    '<-'(
      vary,
      function(tm)
      {
        di <- dim(tm)
        ch <- round(runif(1, 1, di[2]))
        v <- tm[,ch]
        iszero <- which(v=='-')
        if(length(iszero)<=1)
        {
          if(length(iszero))
          {
            p <- iszero
          }
          else
          {
            return(tm)
          }
        }
        else
        {
          p <- sample(iszero, 1)
        }
        if(p==1)
        {
          tmp <- v[2:di[1]]
        }
        else if(p == di[1])
        {
          tmp <- v[1:(di[1]-1)]
        }
        else
        {
          tmp <- character(di[1]-1)
          tmp[1:(p-1)] <- v[1:(p-1)]
          tmp[p:(di[1]-1)] <- v[(p+1):di[1]]
        }
        s <- ifelse(p<5,0,p-5)
        e <- ifelse(p>di[1]-6,di[1]-1,p+5)
        v <- append(tmp, '-', sample(s:e, 1))
        tm[,ch] <- v
        return(tm)
      }
    )
    
    yuzhi <- 100
    mark <- 0
    best <- NA
    hsy <- numeric(10000)
    tim <- -1
    while(yuzhi > 0)
    {
      tim <- tim + 1
      yuzhi <- yuzhi - 1
      '<-'(
        fanyan,
        function(x)
        {
          c <- sample(1:gnum,2)
          return(getChild(group[[c[1]]], group[[c[2]]]))
        }
      )
      #try(ngroup <- parLapply(clist, 1:gnum, fanyan))
      ngroup <- lapply(1:gnum, fanyan)
      group[(gnum+1):(gnum*2)] <- ngroup
      marks <- sapply(group, SP)
      keep <- order(marks, decreasing = TRUE)[1:gnum]
      group <- group[keep]
      marks <- marks[keep]
      maxMark <- max(marks)
      if(maxMark > mark)
      {
        yuzhi <- 100
        best <- group[[which.max(marks)]]
        mark <- maxMark
      }
      hsy[tim] <- maxMark
      
      '<-'(
        eachVary,
        function(ind)
        {
          if(runif(1)>vp)
          {
            return(ind)
          }
          return(vary(ind))
        }
      )
      for(each in 1:gnum)
      {
        group[[each]] <- eachVary(group[[each]])
      }
      #try(group <- parLapply(clist, group, eachVary))
    }
    stopCluster(clist)
    plot(1:tim, hsy[1:tim],'l',xlab='iter',ylab='mark',main='msa-StarClimbGA')
    print(mark)
    return(best)
  }
)

'<-'(
  liang.msaStarGA,
  function(seqs, gnum = 100)
  {
    start <- liang.msaStar(seqs)
    di <- dim(start)
    '<-'(
      genGen,
      function(ch, tm)
      {
        v <- tm[,ch]
        iszero <- which(v=='-')
        lenZero <- length(iszero)
        if(lenZero)
        {
          p <- iszero[round(runif(1,1,lenZero))]
        }
        else
        {
          return(c(0,0))
        }
        s <- ifelse(p<5,0,p-5)
        e <- ifelse(p>di[1]-6,di[1]-1,p+5)
        r <- round(runif(1,s,e))
        return(c(p,r))
      }
    )
    '<-'(
      genInd,
      function(id, tm)
      {
        return(sapply(1:di[2], genGen, tm))
      }
    )
    '<-'(
      SP,
      function(tm)
      {
        di <- dim(tm)
        csum <- di[2]
        mark <- 0
        for(i in 1:(csum-1))
        {
          b1 <- tm[,(i+1):(csum-1)] == tm[,i] #两个位置一致
          b2 <- tm[,(i+1):(csum-1)] == '-' #被比对位置为空格
          m1 <- length(which(!b1 & b2)) #空格对字母
          m2 <- length(which(b1 & !b2)) #字符对字符，而且一致
          #m3 <- length(which(b1 & b2)) #空格对空格
          m4 <- length(which(!b1 & !b2)) #字符对字符，但不 一致
          mark <- mark + m1 * -5 + m2 * 5 + m4 * -4
        }
        mark <- ifelse(mark==0,1,ifelse(mark<0,-1/mark,mark))
        return(mark)
      }
    )
    '<-'(
      calMark,
      function(gen,tm)
      {
        tm <- tranGen_c(tm, gen)
        return(SP(tm))
      }
    )
    '<-'(
      cross,
      function(ind1, ind2)
      {
        ch <- round(runif(1,1,di[2]-1))
        ind <- matrix(0,2,di[2])
        ind[,1:ch] <- ind1[,1:ch]
        ind[,(ch+1):di[2]] <- ind2[,(ch+1):di[2]]
        return(ind)
      }
    )
    '<-'(
      groupEvo,
      function(group, start)
      {
        marks <- matrix(0,gnum,2)
        marks[,1] <- sapply(group, calMark, start)
        '<-'(
          fanyan,
          function(id)
          {
            pch <- sample(1:gnum, 2)
            return(cross(group[[pch[1]]], group[[pch[2]]]))
          }
        )
        bear <- 100
        nowMax <- which.max(marks[,1])
        maxMark <- marks[nowMax,1]
        best <- group[[nowMax]]
        while(bear > 0)
        {
          ngroup <- lapply(1:gnum,fanyan)
          marks[,2] <- sapply(ngroup, calMark, start)
          group[(gnum+1):(gnum*2)] <- ngroup
          keep <- order(marks, decreasing = TRUE)[1:gnum]
          group <- group[keep]
          marks[,1] <- marks[keep]
          nowMax <- which.max(marks[,1])
          if(marks[nowMax,1] <= maxMark)
          {
            bear <- bear  - 1
          }
          else
          {
            print(maxMark)
            maxMark <- marks[nowMax,1]
            best <- group[[nowMax]]
            bear <- 100
          }
        }
        return(best)
      }
    )
    
    last <- 1
    record <- numeric(1000)
    record[1] <- SP(start)
    tim <- 1
    while(last > 0)
    {
      print(record[tim])
      print(tim)
      tim <- tim + 1
      group <- lapply(1:gnum, genInd, start)
      nxt <- groupEvo(group, start)
      start <- tranGen_c(start, nxt)
      record[tim] <- SP(start)
      last <- last - 1
    }
    plot(1:tim, record[1:tim],'l',xlab='iter',ylab='mark')
    
    return(start)
  }
)

#mt <- sapply(1:10, function(x){return(SP(liang.msaStar(temp)))})

