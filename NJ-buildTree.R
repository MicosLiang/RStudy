
source('./RStudy/Needleman-Wunsch.R')

'<-'(
  buildRanges,
  function(ss)
  {
    '<-'(
      JCK,
      function(seq1, seq2)
      {
        temp <- nw(seq1, seq2)
        '<-'(
          biJiao,
          function(x)
          {
            return(x[1]!=x[2])
          }
        )
        d <- sum(as.numeric(unlist(apply(temp, 2, biJiao)))) / dim(temp)[2]
        temp.d <- 4*d/3
        if(temp.d > 1)
        {
          temp.d <- 1
          warning('Liang@The relationship is too far!')
        }
        K <- -0.75 * log(1 - temp.d)
        return(K)
      }
    )
    '<-'(
      easyRange,
      function(seq1, seq2)
      {
        #temp <- nw(seq1, seq2)
        #temp <- liang.dsaSA(seq1, seq2)
        temp <- liang_NW(as.matrix(seq1), as.matrix(seq2), print)
        '<-'(
          biJiao,
          function(x)
          {
            return(x[1]!=x[2])
          }
        )
        d <- sum(as.numeric(unlist(apply(temp, 1, biJiao)))) / dim(temp)[2]
        return(d)
      }
    )
    rangeMatrix <- matrix(0, length(ss), length(ss))
    '<-'(
      calRanges,
      function(i)
      {
        '<-'(
          calRange,
          function(k)
          {
            if(i >= k)
            {
              return(0)
            }
            else
            {
              #return(JCK(ss[[i]], ss[[k]]))
              return(easyRange(ss[[i]], ss[[k]]))
            }
          }
        )
        return(unlist(lapply(1:length(ss),calRange)))
      }
    )
    rangeMatrix <- as.matrix(lapply(1:length(ss),calRanges))
    rangeMatrix <- apply(rangeMatrix, 1, unlist)
    rangeMatrix <- as.matrix(as.dist(rangeMatrix))
    rownames(rangeMatrix) <- names(ss)
    colnames(rangeMatrix) <- names(ss)
    return(rangeMatrix)
  }
)

'<-'(
  NJ,
  function(ss, rangeMatrix, times=50)
  {
    '<-'(
      meshRange,
      function(s1, s2)
      {
        tl <- list(s1, s2)
        tlb <- as.numeric(c(length(s1)>1, length(s2)>1))
        if(sum(tlb) > 0)
        {
          tlw <- which.max(tlb)
          dmi <- meshRange(tl[[tlw]][[1]], tl[[-tlw]])
          dni <- meshRange(tl[[tlw]][[2]], tl[[-tlw]])
          dmn <- meshRange(tl[[tlw]][[1]], tl[[tlw]][[2]])
          return((dmi + dni - dmn) / 2)
        }
        else
        {
          tr <- rangeMatrix[which(ssbak==unlist(s1[1])), which(ssbak==unlist(s2[1]))]
          return(tr)
        }
      }
    )
    '<-'(
      meshBranch,
      function(lst, mn)
      {
        m <- mn[1]
        n <- mn[2]
        nlst <- list(list(lst[[m]], lst[[n]]))
        '<-'(
          temp.op,
          function(i)
          {
            if(i != m && i != n)
            {
              nlst <<- c(nlst, list(lst[[i]]))
            }
          }
        )
        temp.no <- lapply(1:length(lst), temp.op)
        return(nlst)
      }
    )
    '<-'(
      getPairs,
      function(len, no=NULL)
      {
        lst <- list()
        mt <- matrix(0, len, len, dimnames = list(1:len, 1:len))
        if(!is.null(no))
        {
          mt <- mt[-no,-no]
          no = length(no)
        }
        else
        {
          no = 0
        }
        len = len - no
        if(len <= 1)
        {
          return(lst)
        }
        nums <- as.numeric(names(mt[1,]))
        '<-'(
          temp.do,
          function(i)
          {
            '<-'(
              temp.add,
              function(j)
              {
                lst <<- c(lst, list(c(nums[i], nums[j])))
              }
            )
            temp.no <- lapply((i+1):len, temp.add)
          }
        )
        temp.no <- lapply(1:(len-1), temp.do)
        return(lst)
      }
    )
    '<-'(
      sumBranch,
      function(mn, sqs)
      {
        m <- mn[1]
        m <- unlist(m)
        n <- mn[2]
        n <- unlist(n)
        sqslen <- length(sqs)
        if(sqslen <= 2)
        {
          return(Inf)
        }
        '<-'(
          cal1,
          function(i)
          {
            if(i==m || i==n)
            {
              return(0)
            }
            return(meshRange(sqs[[i]], sqs[[m]]) + meshRange(sqs[[i]], sqs[[n]]))
          }
        )
        p1 <- lapply(1:sqslen, cal1)
        p1 <- sum(unlist(p1)) / (2 * (sqslen - 2))
        p2 <- meshRange(sqs[[m]], sqs[[n]]) / 2
        nop <- getPairs(sqslen, c(m, n))
        '<-'(
          cal2,
          function(p)
          {
            return(meshRange(sqs[[p[1]]], sqs[[p[2]]]))
          }
        )
        p3 <- lapply(nop, cal2)
        p3 <- sum(unlist(p3)) / (sqslen - 2)
        ans <- p1 + p2 + p3
        return(ans)
      }
    )
    '<-'(
      nextTree,
      function(lst)
      {
        len <- length(lst)
        if(len==1)
        {
          return(lst)
        }
        allPairs <- getPairs(len)
        '<-'(
          cal,
          function(p)
          {
            return(meshRange(lst[[p[1]]], lst[[p[2]]]))
          }
        )
        s0 <- lapply(allPairs, cal)
        s0 <- sum(unlist(s0)) / (len - 1)
        
        sOthers <- unlist(lapply(allPairs, sumBranch, sqs=lst))
        sAll <- c(sOthers, s0)
        better <- which.min(sAll)
        if(better==length(sAll))
        {
          return(list(lst))
        }
        else
        {
          lock <<- TRUE
          betters <- which(sAll[-length(sAll)] == sAll[better])
          '<-'(
            temp.mesh,
            function(tb)
            {
              bp <- allPairs[[tb]]
              return(meshBranch(lst, bp))
            }
          )
          lsts <- lapply(betters, temp.mesh)
        }
        return(lsts)
      }
    )
    #seqs <- lapply(ss, delSeq)
    ssname <- names(ss)
    ssbak <- lapply(ss, paste, collapse='')
    
    seqs <- lapply(ssbak, list)
    alst <- list(seqs)
    '<-'(
      main.do,
      function()
      {
        talst <- list()
        '<-'(
          main.add,
          function(al)
          {
            talst <<- c(talst, nextTree(al))
          }
        )
        temp.no <- lapply(alst, main.add)
        return(talst)
      }
    )
    lock <- TRUE
    cnt <- 0
    while(lock && cnt < times)
    {
      lock <- FALSE
      cnt <- cnt + 1
      alst <- main.do()
    }
    '<-'(
      calS0,
      function(lst)
      {
        len <- length(lst)
        allPairs <- getPairs(len)
        '<-'(
          cal,
          function(p)
          {
            return(meshRange(lst[[p[1]]], lst[[p[2]]]))
          }
        )
        s0 <- lapply(allPairs, cal)
        s0 <- sum(unlist(s0)) / (len - 1)
        return(s0)
      }
    )
    s0s <- unlist(lapply(alst, calS0))
    mins <- which(s0s == s0s[which.min(s0s)])
    '<-'(
      findSmall,
      function(t)
      {
        return(length(alst[[t]]))
      }
    )
    lens <- unlist(lapply(mins, findSmall))
    ansTree <- alst[[mins[which.min(lens)]]]
    return(ansTree)
  }
)

'<-'(
  drawTree,
  function(tre, ss, dis, fname)
  {
    if(!(require(ggraph) && require(igraph)))
    {
      stop('Liang@You must install ggraph and their requirments!')
    }
    ss <- lapply(ss, paste, collapse='')
    #std <- list(ss[which.min(apply(dis, 1, sum))])
    stt <- tre
    '<-'(
      meshRange,
      function(s1, s2)
      {
        tl <- list(s1, s2)
        tlb <- as.numeric(c(length(s1)>1, length(s2)>1))
        if(sum(tlb) > 0)
        {
          tlw <- which.max(tlb)
          dmi <- meshRange(tl[[tlw]][[1]], tl[[-tlw]])
          dni <- meshRange(tl[[tlw]][[2]], tl[[-tlw]])
          dmn <- meshRange(tl[[tlw]][[1]], tl[[tlw]][[2]])
          return((dmi + dni - dmn) / 2)
        }
        else
        {
          tr <- dis[which(ss==unlist(s1[1])), which(ss==unlist(s2[1]))]
          return(tr)
        }
      }
    )
    nd <- data.frame(name=0, key='')
    '<-'(
      genTree,
      function(tre, ss, node=0, nl=0)
      {
        len <- length(tre)
        d <- data.frame()
        '<-'(
          temp.do,
          function(n)
          {
            d <<- rbind(d, data.frame(from=node, to=node+n+nl))
            if(length(tre[[n]])>1)
            {
              d <<- rbind(d, genTree(tre[[n]], ss, node+n+nl, len+n*10))
              #nd <<- rbind(nd, data.frame(name=node+n+nl, key=abs(meshRange(stt,tre[[n]]))))
              nd <<- rbind(nd, data.frame(name=node+n+nl, key=''))
            }
            else
            {
              #nd <<- rbind(nd, data.frame(name=node+n+nl, key=paste(names(ss)[which(ss == tre[[n]][[1]])], '(', abs(meshRange(stt,tre[[n]])), ')')))
              nd <<- rbind(nd, data.frame(name=node+n+nl, key=names(ss)[which(ss == tre[[n]][[1]])]))
            }
          }
        )
        temp.no <- lapply(1:len,temp.do)
        return(d)
      }
    )
    gt <- genTree(tre,ss)
    #print(nd)
    ggraph(graph_from_data_frame(gt, vertices = nd), layout = 'dendrogram', circular = FALSE) + 
      geom_edge_link() +
      geom_node_point() +
      geom_node_text(aes(label = key,), repel = TRUE) +
      labs(title = fname)+
      theme_void()
  }
)

#dat[which(grepl('com', names(dat)))]