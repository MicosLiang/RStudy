

#Author Liang
#e-mail liangjh21@lzu.edu.cn
#
#使用蒙特卡洛树搜索的中国象棋程序
#经测试，由于R语言递归速度不快和代码优化不好等原因，本程序在未使用神经网络作为
#评价函数的情况下，棋力水平仅相当于刚了解规则的初学者
#使用随机走子作为评价函数的情况下，结果稍好但时间代价无法接受
#程序还在写，有一些未修正的BUG，鲁棒性不高

#source('./RStudy/BP.R')
library(parallel)

'<-'(
  boardToVector,
  function(board)
  {
    '<-'(
      sam,
      function(i)
      {
        return(as.numeric(board==i))
      }
    )
    ans <- lapply(1:14, sam)
    ans <- unlist(ans)
    return(ans)
  }
)

'<-'(
  BPjudge,
  function(board, who)
  {
    v <- boardToVector(board)
  }
)

'<-'(
  judgeboard,
  function(board, who)
  {
    red <- board[board < 8]
    dark <- board[board > 7] %% 7
    
    biao <- c(1.5, 4.5, 4, 9, 2, 2)
    '<-'(
      temp.do,
      function(i, wi)
      {
        length(wi[wi == i]) * biao[i]
      }
    )
    red <- sum(sapply(1:6, temp.do, wi = red))
    dark <- sum(sapply(1:6, temp.do, wi = dark))
    if(red == dark)
    {
      return(who)
    }
    else if(red > dark)
    {
      return(1)
    }
    else
    {
      return(2)
    }
  }
)

'<-'(
  getChoices,
  function(board, chess, isRed)
  {
    zi <- chess[1] %% 7
    zi <- ifelse(zi, zi, 7)
    pos <- chess[2:3]
    if(zi == 1)
    {
      ops <- matrix(c(1,2,2,1,pos[2]-1,pos[2]+1), 3, 2)
    }
    else if(zi == 2)
    {
      ops <- matrix(1, 27, 2)
      ops[10:18,1] <- 2
      ops[19:27,1] <- 3
      ops[,2] <- 1:9
    }
    else if(zi == 3)
    {
      ops <- matrix(1, 8, 2)
      ops[5:8,1] <- 3
      temp <- c(pos[2]-1,pos[2]-2,pos[2]+1,pos[2]+2)
      ops[,2] <- temp
    }
    else if(zi == 4)
    {
      ops <- matrix(1, 27, 2)
      ops[10:18,1] <- 2
      ops[19:27,1] <- 3
      ops[,2] <- 1:9
    }
    else if(zi == 5)
    {
      ops <- matrix(1, 4, 2)
      ops[3:4,1] <- 3
      temp <- c(pos[2]-2,pos[2]+2)
      ops[,2] <- temp
    }
    else if(zi == 6)
    {
      ops <- matrix(1, 4, 2)
      ops[3:4,1] <- 3
      temp <- c(pos[2]-1,pos[2]+1)
      ops[,2] <- temp
    }
    else if(zi == 7)
    {
      ops <- matrix(c(1,2,2,1,1,pos[2]-1,pos[2]+1,1), 4, 2)
    }
    else
    {
      return()
    }
    '<-'(
      checkop,
      function(opi)
      {
        op <- ops[opi,]
        after <- getAfter(zi, pos, op, isRed)
        afters[opi,] <<- after
        can <- canMove(board, zi, pos, after, op, isRed)
        return(can)
      }
    )
    
    afters <- ops
    cans <- sapply(1:length(ops[,1]), checkop)
    ops <- ops[cans,]
    afters <- afters[cans,]
    left <- length(cans[cans])
    if(left)
    {
      if(left == 1)
      {
        before <- pos
        luo <- as.matrix(c(before, afters))
      }
      else 
      {
        before <- matrix(pos, left, 2, byrow = TRUE)
        luo <- cbind(before, afters)
        luo <- t(luo)
      }
      return(luo)
    }
    else
    {
      return(NULL)
    }
  }
)

'<-'(
  getPossibles,
  function(board, who)
  {
    '<-'(
      colsDo,
      function(i)
      {
        '<-'(
          eachDo,
          function(j)
          {
            if(who == 1)
            {
              if(board[i,j] <= 7 && board[i, j] > 0)
              {
                return(getChoices(board, c(board[i,j], i, j), isRed))
              }
            }
            else
            {
              if(board[i,j] > 7)
              {
                return(getChoices(board, c(board[i,j], i, j), isRed))
              }
            }
          }
        )
        return(sapply(1:9,eachDo))
      }
    )
    isRed <- who == 1
    allchoices <- sapply(1:10,colsDo)
    allchoices <- unlist(allchoices)
    if(is.null(allchoices))
    {
      print(board)
      return(NULL)
    }
    num <- length(allchoices) / 4
    allchoices <- matrix(allchoices, num, 4, byrow = TRUE)
    return(allchoices)
  }
)

'<-'(
  judgeGameOver,
  function(board, who)
  {
    if(!any(board==7))
    {
      return(2)
    }
    else if(!any(board==14))
    {
      return(1)
    }
    else
    {
      red <- which(board==7)
      dark <- which(board==14)
      red <- ceiling(red / 10)
      dark <- ceiling(dark / 10)
      if(red == dark)
      {
        dui <- board[,red]
        dui <- dui[dui != 0]
        if(length(dui) == 2)
        {
          if(who == 1)
          {
           ans <- 1 
          }
          else
          {
            ans <- 2
          }
          return(ans)
        }
      }
    }
    return(FALSE)
  }
)

'<-'(
  randomGame,
  function(board = newBoard(), step = 1)
  {
    who <- step %% 2
    isover <- judgeGameOver(board, who)
    if(isover)
    {
      return(isover)
    }
    board <- randomGo(board, who)
    return(randomGame(board, step + 1))
  }
)

'<-'(
  randomGo,
  function(board, who)
  {
    possibiles <- getPossibles(board, who)
    if(is.null(possibiles))
    {
      warning('noting to do')
      return(board)
    }
    get <- sample(1:length(possibiles[,1]), 1)
    pos <- possibiles[get,1:2]
    after <- possibiles[get,3:4]
    board[after[1], after[2]] <- board[pos[1], pos[2]]
    board[pos[1], pos[2]] <- 0
    return(board)
  }
)

'<-'(
  MTSC,
  function(board, who, long = 10000, ping = randomGame)
  {
    UC <- 0.5
    
    who <- ifelse(who == 1, 1, 0)
    get <- getPossibles(board, who)
    if(is.null(get))
    {
      warning('nothing to do')
      return(board)
    }
    #tree <- list(numeric(long * 5), list(get), list(board), list(c()), list(c()), 1)
    tree <- list(c(0, 0, 0, 1, who), list(get), list(board), list(c()), list(c()), 1)
    names(tree) <- c('dat', 'ps', 'bs', 'cs', 'had', 'deep')
    tree$dat[4] <- 1
    tree$dat[5] <- who
    '<-'(
      getDat,
      function(i)
      {
        s <- (i-1)*5 + 1
        return(tree$dat[s:(s+4)])
      }
    )
    '<-'(
      selection,
      function(now = 1)
      {
        childs <- unlist(tree$cs[[now]])
        allway <- length(tree$ps[[now]][,1])
        if(length(childs) == allway)
        {
          '<-'(
            UCT,
            function(i)
            {
              cdat <- getDat(i)
              xc <- cdat[2] / cdat[3]
              if(pdat[5] != who)
              {
                xc <- 1 - xc
              }
              ans <- xc + (UC * sqrt(log(pdat[3]) / cdat[3]))
              return(ans)
            }
          )
          pdat <- getDat(now)
          ucts <- sapply(childs, UCT)
          cc <- childs[which.max(ucts)]
          return(selection(cc))
        }
        else
        {
          return(now)
        }
      }
    )
    '<-'(
      expansion,
      function(p)
      {
        pdat <- getDat(p)
        nwho <- as.numeric(!pdat[5])
        tree$deep <<- tree$deep + 1
        tree$dat <<- c(tree$dat, c(p, 0, 0, tree$deep, nwho))
        tree$cs <<- c(tree$cs, list(c()))
        tree$had <<- c(tree$had, list(c()))
        
        tree$cs[[pdat[[4]]]] <<- c(tree$cs[[pdat[[4]]]] ,tree$deep)
      }
    )
    '<-'(
      simluation,
      function()
      {
        ldat <- getDat(tree$deep)
        pdat <- getDat(ldat[1])
        choice <- tree$ps[[pdat[4]]]
        not <- 1:length(choice[,1])
        had <- unlist(tree$had[[pdat[4]]])
        if(!is.null(had))
        {
          not <- not[-had]
        }
        this <- sample(not, 1)
        tree$had[[pdat[4]]] <<- c(had, this)
        choice <- choice[this,]
        
        nboard <- tree$bs[[pdat[4]]]
        nboard[choice[3], choice[4]] <- nboard[choice[1], choice[2]]
        nboard[choice[1], choice[2]] <- 0
        tree$bs[[tree$deep]] <<- nboard
        nwho <- as.numeric(!pdat[5])
        getp <- getPossibles(nboard, nwho)
        tree$ps[[tree$deep]] <<- getp
        
        ifover <- judgeGameOver(nboard, nwho)
        if(ifover)
        {
          ans <- ifover
        }
        else if(is.null(getp))
        {
          ans <- nwho
        }
        else
        {
          ans <- ping(nboard, pdat[5]) 
        }
        ans <- ifelse(ans==1, 1, 0)
        if(ans==who)
        {
          ans <- c(1,1)
        }
        else
        {
          ans <- c(0,1)
        }
        return(ans)
      }
    )
    '<-'(
      backpropagain,
      function(ans, now = tree$deep)
      {
        st <- (now-1)*5
        nxt <- tree$dat[st + 1]
        tree$dat[st + 2] <<- tree$dat[st + 2] + ans[1]
        tree$dat[st + 3] <<- tree$dat[st + 3] + ans[2]
        if(nxt == 0)
        {
          return(TRUE)
        }
        return(backpropagain(ans, nxt))
      }
    )
    '<-'(
      main.do,
      function(u)
      {
        nc <- selection()
        expansion(nc)
        sa <- simluation()
        backpropagain(sa)
      }
    )
    
    temp.no <- sapply(1:long, main.do)
    
    childs <- tree$cs[[1]]
    '<-'(
      choose.do,
      function(i)
      {
        cdat <- getDat(i)
        return(cdat[3])
      }
    )
    best <- childs[which.max(sapply(childs, choose.do))]
    return(tree$bs[[best]])
  }
)

'<-'(
  minmax,
  function(board, who, deep)
  {
    
  }
)

'<-'(
  gamePlay,
  function(num = 100)
  {
    '<-'(
      main.play,
      function(id, board = newBoard(), step = 1, hsty = NULL)
      {
        who <- step %% 2
        isover <- judgeGameOver(board, who)
        if(isover)
        {
          return(isover)
          if(isover == 1)
          {
            #jie <- matrix(c(1,0), step-1, 2, byrow = TRUE)
          }
          else
          {
            #jie <- matrix(c(0,1), step-1, 2, byrow = TRUE)
          }
          #hsty <- cbind(jie, hsty)
          return(hsty)
        }
        if(is.null(hsty))
        {
          #hsty <- boardToVector(board)
        }
        else
        {
          #hsty <- rbind(hsty, boardToVector(board))
        }
        if(who == 1)
        {
          board <- MTSC(board, who, long = 1000, judgeboard)
        }
        else
        {
          board <- randomGo(board, who)
        }
        return(main.play(id ,board, step+1, hsty))
      }
    )
    #clist <- makeCluster(detectCores())
    #tmp.no <- lapply(ls(".GlobalEnv"), (function(t){return(clusterExport(clist, t, envir = environment()))}))
    #ans <- parSapply(clist, 1:num, main.play)
    #stopCluster(clist)
    ans <- sapply(1:num, main.play)
    return(ans)
  }
)

'<-'(
  moveChess,
  function(board, position, op)
  {
    #判断是哪一方
    isRed <- !(board[position[1], position[2]] > 7)
    
    #判断是否有这个子
    if(board[position[1], position[2]])
    {
      board <- chessMove(board, position, op, isRed)
    }
    else
    {
      warning('No chess')
    }
    return(board)
  }
)

#获取移动后的位置，不考虑能不能移动
'<-'(
  getAfter,
  function(chess, position, op, isRed)
  {
    if(chess == 1)
    {
      if(op[1] == 1)
      {
        way <- op[2] * ifelse(isRed, -1, 1)
        return(c((position[1] + way), position[2]))
      }
      else if(op[1] == 2)
      {
        way <- op[2]
        return(c(position[1], way))
      }
    }
    else if(chess == 2)
    {
      if(op[1] == 1)
      {
        way <- op[2] * ifelse(isRed, -1, 1)
        return(c((position[1] + way), position[2]))
      }
      else if(op[1] == 2)
      {
        way <- op[2]
        return(c(position[1], way))
      }
      else if(op[1] == 3)
      {
        way <- op[2] * ifelse(isRed, -1, 1)
        return(c((position[1] - way), position[2]))
      }
    }
    else if(chess == 3)
    {
      way <- (abs(op[2] - position[2]) %% 2) + 1
      way <- way * ifelse(isRed, -1, 1)
      if(op[1] == 1)
      {
        return(c((position[1] + way), op[2]))
      }
      else if(op[1] == 3)
      {
        return(c((position[1] - way), op[2]))
      }
    }
    else if(chess == 4)
    {
      if(op[1] == 1)
      {
        way <- op[2] * ifelse(isRed, -1, 1)
        return(c((position[1] + way), position[2]))
      }
      else if(op[1] == 2)
      {
        way <- op[2]
        return(c(position[1], way))
      }
      else if(op[1] == 3)
      {
        way <- op[2] * ifelse(isRed, -1, 1)
        return(c((position[1] - way), position[2]))
      }
    }
    else if(chess == 5)
    {
      if(op[1] == 1)
      {
        way <- 2 * ifelse(isRed, -1, 1)
        return(c((position[1] + way), op[2]))
      }
      else if(op[1] == 3)
      {
        way <- 2 * ifelse(isRed, -1, 1)
        return(c((position[1] - way), op[2]))
      }
    }
    else if(chess == 6)
    {
      if(op[1] == 1)
      {
        way <- ifelse(isRed, -1, 1)
        return(c((position[1] + way), op[2]))
      }
      else if(op[1] == 3)
      {
        way <- ifelse(isRed, -1, 1)
        return(c((position[1] - way), op[2]))
      }
    }
    else if(chess == 7)
    {
      if(op[1] == 1)
      {
        way <- op[2] * ifelse(isRed, -1, 1)
        return(c((position[1] + way), position[2]))
      }
      else if(op[1] == 2)
      {
        way <- op[2]
        return(c(position[1], way))
      }
      else if(op[1] == 3)
      {
        way <- op[2] * ifelse(isRed, -1, 1)
        return(c((position[1] - way), position[2]))
      }
    }
    return(position)
  }
)

#判断能不能移动
'<-'(
  canMove,
  function(board, chess, position, after, op, isRed)
  {
    '<-'(
      howMuch,
      function(way)
      {
        if(way == 2)
        {
          road <- board[position[1], position[2]:after[2]]
        }
        else
        {
          road <- board[position[1]:after[1], position[2]]
        }
        road <- road[-1]
        road <- rev(road)[-1]
        num <- length(road[road!=0])
        return(num)
      }
    )
    '<-'(
      isout,
      function()
      {
        ans <- after[1] > 0 && after[1] < 11 && after[2] > 0 && after[2] < 10
        return(!ans)
      }
    )
    #
    if(isout())
    {
      return(FALSE)
    }
    
    #判断要走的位置上是否是自己的子
    if(board[after[1], after[2]] != 0)
    {
      if((board[after[1], after[2]] > 7) != isRed)
      {
        return(FALSE)
      }
    }
    
    #各个棋子
    if(chess == 1)
    {
      if(op[1] == 1)
      {
        #兵只能走一格
        if(op[2] != 1)
        {
          return(FALSE)
        }
        
        return(TRUE)
      }
      else if(op[1] == 2)
      {
        #看过没过河
        if(isRed)
        {
          if(position[1] > 5)
          {
            return(FALSE)
          }
        }
        else
        {
          if(position[1] < 6)
          {
            return(FALSE)
          }
        }
        
        #兵只能走一格
        if(abs(op[2] - position[2]) != 1)
        {
          return(FALSE)
        }
        
        return(TRUE)
      }
      else
      {
        #兵不能后退
        return(FALSE)
      }
    }
    else if(chess == 2)
    {
      tai <- howMuch(op[1])
      if(tai)
      {
        if(tai == 1)
        {
          return(board[after[1], after[2]] != 0)
        }
        return(FALSE)
      }
      else
      {
        return(board[after[1], after[2]] == 0)
      }
    }
    else if(chess == 3)
    {
      ge <- position[2] - op[2]
      a_ge <- abs(ge)
      if(a_ge == 1 || a_ge == 2)
      {
        ka <- howMuch(a_ge)
        return(ka == 0)
      }
      return(FALSE)
    }
    else if(chess == 4)
    {
      ge <- howMuch(op[1])
      if(ge == 0)
      {
        return(TRUE)
      }
      return(FALSE)
    }
    else if(chess == 5)
    {
      if(isRed)
      {
        if(after[1] < 6)
        {
          return(FALSE)
        }
      }
      else
      {
        if(after[1] > 5)
        {
          return(FALSE)
        }
      }
      
      if(abs(position[2] - op[2]) == 2)
      {
        zhong <- (position + after) / 2
        return(board[zhong[1],zhong[2]] == 0)
      }
      return(FALSE)
    }
    else if(chess == 6)
    {
      if(after[2] < 4 || after[2] > 6)
      {
        return(FALSE)
      }
      
      if(isRed)
      {
        if(after[1] < 8)
        {
          return(FALSE)
        }
      }
      else
      {
        if(after[1] > 3)
        {
          return(FALSE)
        }
      }
      
      if(abs(op[2] - position[2]) != 1)
      {
        return(FALSE)
      }
      
      return(TRUE)
    }
    else if(chess == 7)
    {
      if(after[2] < 4 || after[2] > 6)
      {
        return(FALSE)
      }
      
      if(isRed)
      {
        if(after[1] < 8)
        {
          return(FALSE)
        }
      }
      else
      {
        if(after[1] > 3)
        {
          return(FALSE)
        }
      }
      
      if(op[1] == 2)
      {
        return(abs(op[2] - position[2]) == 1)
      }
      else
      {
        return(op[2] == 1)
      }
    }
  }
)

'<-'(
  chessMove,
  function(board, position, op, isRed)
  {
    #主要逻辑
    chess <- board[position[1], position[2]] %% 7
    chess <- ifelse(chess, chess, 7)
    after <- getAfter(chess, position, op, isRed)
    if(canMove(board, chess, position, after, op, isRed))
    {
      #移动操作
      board[after[1], after[2]] <- chess + ifelse(isRed, 0, 7)
      board[position[1], position[2]] <- 0
    }
    else
    {
      warning('ill op')
    }
    return(board)
  }
)

'<-'(
  newBoard,
  function()
  {
    ans <- c(11,10,12,13,14,13,12,10,11,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,9,0,8,0,8,0,8,0,8,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,1,0,1,0,2,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,4,3,5,6,7,6,5,3,4)
    ans <- matrix(ans, 10 ,9, byrow = TRUE)
    return(ans)
  }
)