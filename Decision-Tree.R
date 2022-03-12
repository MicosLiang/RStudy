'<-'(
  Node,
  function(isLeaf, children, values, father, label, id)
  {
    temp <- list(isLeaf, children, father, label, id ,'', values)
    names(temp) <- c('isLeaf', 'children', 'father', 'label', 'id', 'judge', 'values')
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
  function(Data, mainTree, father, give)
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
    }
    if(length(levels(factor(Data[,1])))==1)
    {
      mainTree[[nid]]$isLeaf <- TRUE
      return(mainTree) 
    }
    if(length(attris) == 0)
    {
      mainTree[[nid]]$isLeaf <- TRUE
      return(mainTree)
    }
    gainArg <- lapply(attris, Gain, D=Data)
    aBest <- attris[which.max(gainArg)]
    mainTree[[nid]]$judge <- aBest
    abv <- levels(factor(Data[,aBest]))
    '<-'(
      doEta,
      function(v)
      {
        Dv <- Data[which(Data[,aBest] %in% v),]
        if(is.null(dim(Dv)))
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
            Dv <- as.matrix(Dv)
          }
          mainTree <<- TreeGenerate(Dv, mainTree, nid, v)
        }
      }
    )
    lapply(abv, doEta)
    return(mainTree)
  }
)