############################
#Author : liang
#time : 2022-03-07
############################

'<-'(
  nw,
  function(s1, s2)
  {
    #temp.mark <- c(10, -1, -3, -4, -5, -1, 7, -5, -3, -5, -3, -5, 9, 0, -5, -4, -3, 0, 8, -5, -5, -5,-5, -5, -5)
    #temp.name <- c('A', 'G', 'C', 'T', '-')
    #temp.mark <- matrix(data=temp.mark, nrow = 5, ncol = 5,dimnames = list(temp.name, temp.name))
    '<-'(
      changeMark,
      function(ni, nj)
      {
        if(ni == nj)
        {
          return(5)
        }
        else
        {
          return(-4)
        }
      }
    )
    gap <- -5
    
    len1 <- length(s1) + 1
    len2 <- length(s2) + 1
    long <- len1 + len2 - 1
    rn <- c('-', s1)
    cn <- c('-', s2)
    mark <- matrix(data=0, nrow=len1, ncol=len2, dimnames = list(rn ,cn)) #计分矩阵
    road <- matrix(data=0, nrow=len1, ncol=len2, dimnames = list(rn, cn)) #朔源矩阵 1平 2左 3上
    
    nowi <- 0
    nowj <- 0
    namei <- '-'
    namej <- '-'
    '<-'(
      sMark,
      function()
      {
        if(nowi==1)
        {
          if(nowj==1)
          {
            return(0)
          }
          else
          {
            return(gap * (nowj-1))
          }
        }
        else if(nowj==1)
        {
          return(gap * (nowi-1))
        }
        
        lp <- mark[nowi-1, nowj-1] +  changeMark(namei, namej)
        up <- mark[nowi-1, nowj] + gap
        left <- mark[nowi, nowj-1] + gap
        
        return(c(lp, left,  up))
      }
    )    
    
    '<-'(
      writeMark,
      function(n)
      {
        nowi <<- (n %% len1) + 1
        nowj <<- floor(n / len1) + 1
        namei <<- rn[nowi]
        namej <<- cn[nowj]
        #print(c(nowi, nowj, namei, namej))
        
        temp <- sMark()
        #print(temp)
        if(length(temp)==1)
        {
          if(temp==0)
          {
            mark[nowi, nowj] <<- 0
            road[nowi, nowj] <<- 0
          }
          else if(nowi == 1)
          {
            mark[nowi, nowj] <<- temp
            road[nowi, nowj] <<- 2
          }
          else if(nowj == 1)
          {
            mark[nowi, nowj] <<- temp
            road[nowi, nowj] <<- 3
          }
        }
        else
        {
          mark[nowi, nowj] <<- max(temp)
          road[nowi, nowj] <<- which.max(temp)
        }
      }
    )
    
    temp.no <- lapply(0:((len1*len2)-1), writeMark)
    #print(mark)
    
    ans <- matrix(data='-', nrow = 2, ncol = long)
    #nowj <- which.max(mark[len1,])
    #print(road)
    '<-'(
      findBest,
      function(i, j)
      {
        temp.way <- road[i, j]
        #print(temp.way)
        if(temp.way==0)
        {
          return(TRUE)
        }
        else if(temp.way == 1)
        {
          ans[1, long] <<- rn[i]
          ans[2, long] <<- cn[j]
          long <<- long - 1
          return(findBest(i-1, j-1))
        }
        else if(temp.way == 2)
        {
          ans[2, long] <<- cn[j]
          long <<- long - 1
          return(findBest(i, j-1))
        }
        else if(temp.way == 3)
        {
          ans[1, long] <<- rn[i]
          long <<- long - 1
          return(findBest(i-1, j))
        }
        else
        {
          return(FALSE)
        }
      }
    )
    
    temp.right <- findBest(nowi, nowj)
    
    #print(temp.right)
    return(ans[,-(1:long)])
  }
)

'<-'(
  delSeq,
  function(s)
  {
    ans <- as.character(unlist(strsplit(s,split='')))
    return(ans)
  }
)

'<-'(
  read.fasta,
  function(f)
  {
    ans <- list()
    temp.in <- readLines(f)
    long <- length(temp.in)
    temp.ans <- c()
    temp.names <- c()
    '<-'(
      delLine,
      function(w)
      {
        temp.c <- delSeq(temp.in[w])
        if(length(temp.c)==0 || w==long)
        {
          ans <<- c(ans, list(temp.ans))
          temp.ans <<- c()
        }
        else if(length(which(temp.c=='>')))
        {
          temp.names <<- c(temp.names, paste(temp.c[-1], collapse = ''))
        }
        else
        {
          temp.ans <<- c(temp.ans, temp.c)
        }
      }
    )
    lapply(1:long, delLine)
    names(ans) <- temp.names
    return(ans)
  }
)
