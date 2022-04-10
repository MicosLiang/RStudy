Rcpp::sourceCpp('./RStudy/five_chess_game.cpp')

'<-'(
  randomGame,
  function(id=0, board = matrix(0,15,15), isDark = TRUE)
  {
    spaces <- which(board==0)
    if(length(spaces))
    {
      spaces <- spaces - 1
      #spaces <- sample(spaces, length(spaces))
      #ans <- random_game(board, isDark, spaces)
      len1 <- length(spaces)
      good <- chooseGood(board, isDark, 0)
      spaces <- spaces[-good]
      len2 <- length(spaces)
      spaces <- sample(spaces, len2)
      len3 <- length(good)
      good <- good - 1
      ans <- integer(len1)
      ans[1:len3] <- good
      ans[(len3+1):(len2+len3)] <- spaces
      ans <- random_game(board, isDark, ans)
      return(ans)
    }
    else
    {
      return(ifelse(isDark, 2, 1))
    }
  }
)

'<-'(
  chooseSpaces,
  function(board, isDark)
  {
    spaces <- which(board==0)
    if(length(spaces))
    {
      spaces <- spaces - 1
      len1 <- length(spaces)
      good <- chooseGood(board, isDark, 0)
      spaces <- spaces[-good]
      len2 <- length(spaces)
      spaces <- sample(spaces, len2)
      len3 <- length(good)
      good <- good - 1
      ans <- integer(len1)
      ans[1:len3] <- good
      ans[(len3+1):(len2+len3)] <- spaces
      return(ans)
    }
    else
    {
      return(c())
    }
  }
)

'<-'(
  chooseGood,
  function(board, isDark, deep)
  {
    not <- board==0
    board <- as.vector(board!=0)
    ch <- board
    '<-'(
      rec,
      function(i)
      {
        b <- board[i]
        if(b)
        {
          if(i!=1)
          {
            ch[i-1] <<- TRUE        
          }
          if(i!=225)
          {
            ch[i+1] <<- TRUE
          }
          if(i>16)
          {
            ch[i-15] <<- TRUE
          }
          if(i<209)
          {
            ch[i+15] <<- TRUE
          } 
        }
      }
    )
    sapply(1:225, rec)
    board <- board | ch
    if(deep > 1)
    {
      not <- not & board
      not <- which(not)
      return(not)
    }
    sapply(1:225, rec)
    board <- board | ch
    sapply(1:225, rec)
    board <- board | ch
    not <- not & board
    if(!any(not))
    {
      not[123] <- TRUE
    }
    not <- which(not)
    return(not)
  }
)

'<-'(
  putChess,
  function(board, p, zi)
  {
    board[p] <- zi
    nboard <- board
    return(nboard)
  }
)

'<-'(
  aiChess,
  function(board=matrix(0,15,15), isDark=TRUE, times=10000)
  {
    ans <- MCTS(board, isDark, putChess, chooseGood, chooseSpaces, print, times)
    #board[ans] <- ifelse(isDark, 1, 2);
    pos <- c((ans-1)%%15 + 1, floor((ans-1)/15) + 1)
    return(pos)
  }
)

'<-'(
  boardToData,
  function(board)
  {
    v1 <- as.numeric(board==1)
    v2 <- as.numeric(board==2)
    return(c(v1, v2))
  }
)

'<-'(
  gamePlay,
  function()
  {
    board <- matrix(0,15,15)
    who <- TRUE
    ans <- 0
    step <- 1
    '<-'(
      main.play,
      function(id)
      {
        if(ans)
        {
          return(rep(0, 550))
        }
        pos <- aiChess(board, who)
        board[pos[1], pos[2]] <<- ifelse(who, 1, 2)
        ans <<- judge_game_over(board, pos[1], pos[2])
        print(ans)
        step <<- step + 1
        who <<- !who
        print(board)
        return(boardToData(board))
      }
    )
    dat <- sapply(1:225, main.play)
    print(ans)
    print(board)
    if(ans)
    {
      dat <- t(dat[,1:step])
      yes <- matrix(c(ifelse(ans==1, 1, 0), ifelse(ans==1, 0, 1)), step, 2, byrow = TRUE)
      return(list(dat, yes))
    }
    return(NULL)
  }
)

'<-'(
  randomGen,
  function(long, isDNA = TRUE, isTwo = FALSE)
  {
    db <- c('A','G','C','T')
    rb <- c('A','G','C','U')
    if(isDNA)
    {
      base <- db
    }
    else
    {
      base <- rb
    }
    tmp <- floor(rnorm(long, 14, 10)) %% 4 + 1
    #gen <- sapply(1:long, function(i){return(base[tmp[i]])})
    gen <- character(long)
    sapply(1:4, function(i){gen[which(tmp==i)] <<- base[i]})
    if(isTwo && isDNA)
    {
      gen2 <- compleDNA(gen)
      return(matrix(c(gen, gen2), long, 2))
    }
    return(gen)
  }
)

'<-'(
  compleDNA,
  function(dna)
  {
    base <- c('A','G','C','T')
    len <- length(dna)
    ans <- sapply(dna, function(b){return(base[5 - which(base==b)])})
    return(ans)
  }
)

'<-'(
  revcomDNA,
  function(dna)
  {
    return(compleDNA(rev(dna)))
  }
)

'<-'(
  tranDNA,
  function(dna)
  {
    db <- c('A','G','C','T')
    rb <- c('A','G','C','U')
    len <- length(dna)
    mrna <- sapply(dna, function(b){return(rb[5 - which(db==b)])})
    return(mrna)
  }
)

'<-'(
  tranRNA,
  function(gen)
  {
    #keys() #密码子矩阵，匹配函数
    if(is.matrix(gen))
    {
      #判断是否是DNA，双链
      ans1 <- tranRNA(gen[,1])
      ans2 <- tranRNA(gen[,2])
      return(c(ans1, ans2))
    }
    if(any(gen=='T'))
    {
      #判断是否是DNA单链
      return(tranRNA(tranDNA(gen)))
    }
    #其他情况默认RNA单链
    '<-'(
      readRNA,
      function(rna)
      {
        p <- 1 #指针
        len <- length(rna)
        out <- character(floor(len/3))
        now <- character(3)
        while(p <= len)
        {
          #沿序列搜索
          id <- p %% 3
          now[ifelse(id,id,3)] <- rna[p]
          p <- p + 1
          if(p %% 3 == 0)
          {
            #out[p/3] <- keys(now)
          }
        }
        return(out)
      }
    )
    #考虑三种读框
    tmp <- sapply(1:3, function(i){return(gen[-(1:i)])})
    ans <- lapply(tmp, readRNA)
    return(ans)
  }
)