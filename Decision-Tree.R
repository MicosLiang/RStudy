'<-'(
  Node,
  function(isLeaf, children, values, father, label, id)
  {
    temp <- list(isLeaf, children, father, label, id ,'', values, 0)
    names(temp) <- c('isLeaf', 'children', 'father', 'label', 'id', 'judge', 'values', 'deep')
    return(temp)
  }
)

'<-'(
  Ent,
  function(D)
  {
    if(is.null(dim(D)))
    {
      return(0)
    }
    D <- as.matrix(table(D[,1]))
    al <- sum(D[,1])
    '<-'(
      cal,
      function(x)
      {
        x <- x / al
        return(x * log(x, 2))
      }
    )
    ps <- apply(D, 1, cal)
    #names(ps) <- names(D[,1])
    return(-sum(ps))
  }
)

'<-'(
  Gain,
  function(arf, D)
  {
    av <- levels(factor(D[,arf]))
    rowNumAll <- dim(D)[1]
    '<-'(
      cal_1,
      function(v)
      {
        Dv <- D[which(D[,arf] %in% v),]
        return((dim(Dv)[1]/rowNumAll) * Ent(Dv))
      }
    )
    vs <- lapply(av, cal_1)
    vs <- as.numeric(unlist(vs))
    return(Ent(D)-sum(vs))
  }
)

'<-'(
  TreeGenerate,
  function(Data, mainTree=list(), father=0, give=0, As=list())
  {
    attris <- names(Data[1,-1])
    clas <- factor(Data[,1])
    bigClas <- names(which.max(table(clas)))
    nid <- length(mainTree)+1
    mainTree[[nid]] <- Node(FALSE, c(), c(), father, bigClas, nid)
    if(father != 0)
    {
      mainTree[[father]]$children <- c(mainTree[[father]]$children, nid)
      mainTree[[father]]$values <- c(mainTree[[father]]$values, give)
      mainTree[[nid]]$deep <- mainTree[[father]]$deep + 1
    }
    else
    {
      '<-'(
        levelsGen,
        function(c)
        {
          c <- levels(factor(c))
          return(c)
        }
      )
      As <- apply(Data, 2, levelsGen)
    }
    if(length(levels(clas))==1)
    {
      mainTree[[nid]]$isLeaf <- TRUE
      return(mainTree)
    }
    if(length(attris) == 0)
    {
      mainTree[[nid]]$isLeaf <- TRUE
      return(mainTree)
    }
    if(mainTree[[nid]]$deep %% 100 == 0)
    {
      print(mainTree[[nid]]$deep)
      print(Gain(names(Data[1,])[1], Data))
    }
    if((Gain(names(Data[1,])[1], Data) < 0.82 && mainTree[[nid]]$deep > 5) || mainTree[[nid]]$deep > 100)
    {
      mainTree[[nid]]$isLeaf <- TRUE
      return(mainTree)
    }
    gainArg <- lapply(attris, Gain, D=Data)
    aBest <- attris[which.max(gainArg)]
    mainTree[[nid]]$judge <- aBest
    abv <- levels(factor(Data[,aBest]))
    #abz <- unlist(As[aBest])
    abz <- c('有','无')
    '<-'(
      doEta,
      function(v)
      {
        Dv <- Data[which(Data[,aBest] %in% v),]
        if(is.null(dim(Dv)))
        {
          Dv <- t(as.matrix(Dv))
        }
        if(!(v %in% abv))
        {
          child <- Node(TRUE, c(), c(), nid, bigClas, length(mainTree)+1)
          mainTree[[child$id]] <<- child
          mainTree[[nid]]$children <<- c(mainTree[[nid]]$children, child$id)
          mainTree[[nid]]$values <<- c(mainTree[[nid]]$values, v)
          return()
        }
        else
        {
          Dv <- Dv[,-which(names(Dv[1,])==aBest)]
          if(is.null(dim(Dv)))
          {
            Dv <- t(as.matrix(Dv))
          }
          mainTree <<- TreeGenerate(Dv, mainTree, nid, v, As)
        }
      }
    )
    lapply(abz, doEta)
    return(mainTree)
  }
)

'<-'(
  TreeJudge,
  function(mainTree, datIn, now=1)
  {
    if(mainTree[[now]]$isLeaf)
    {
      return(mainTree[[now]]$label)
    }
    else
    {
      judge <- mainTree[[now]]$judge
      if(datIn[judge] %in% mainTree[[now]]$values)
      {
        nxt <- mainTree[[now]]$children[which(mainTree[[now]]$values==datIn[judge])]
        return(TreeJudge(mainTree, datIn, nxt))
      }
      else
      {
        return('')
      }
    }
  }
)

'<-'(
  TreeJudges,
  function(datIns, mainTree, by=1, test=TRUE)
  {
    if(test)
    {
      right <- as.vector(datIns[,1])
      datIns <- datIns[,-1]
    }
    temp <- unlist(apply(datIns, by, FUN = TreeJudge, mainTree=mainTree))
    if(test)
    {
      print(length(which(temp==right))/length(right))
      #print(length(which(temp==right))/length(temp[which(temp=='y')])) #召回率，对测试集无实际意义
    }
    return(temp)
  }
)

'<-'(
  delInputStr,
  function(s, li)
  {
    temp <- rep('无',length(li))
    names(temp) <- li
    temp[unlist(lapply(li, grepl, x = s))] <- '有'
    return(temp)
  }
)

'<-'(
  delInputStrs,
  function(ss)
  {
    tm <- matrix(unlist(lapply(ss, delInputStr, li=names(qq[1,-1]))), length(ss), length(names(qq[1,-1])), dimnames = list(1:length(ss), names(qq[1,-1])))
    return(tm)
  }
)