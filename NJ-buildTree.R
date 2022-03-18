
#source('./RStudy/Needleman-Wunsch.R')

'<-'(
  NJ,
  function(ss, times=50)
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
          return(JCK(s1[[1]], s2[[1]]))
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
        sAll <- c(s0, sOthers)
        better <- which.min(sAll)
        if(better==1)
        {
          print('build compelete')
          print(sAll)
          return(lst)
        }
        else
        {
          bp <- allPairs[[(better-1)]]
          lst <- meshBranch(lst, bp)
        }
        return(lst)
      }
    )
    seqs <- lapply(ss, delSeq)
    seqs <- lapply(seqs, list)
    '<-'(
      main.do,
      function(u)
      {
        seqs <<- nextTree(seqs)
      }
    )
    temp.no <- lapply(1:times, main.do)
    return(seqs)
  }
)