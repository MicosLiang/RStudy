require(Rcpp)
sourceCpp('./RStudy/ML-buildTree.cpp')
source('./RStudy/sequence_alignment.R')
'<-'(
  MP,
  function(ss)
  {
    b <- round(dim(ss)[2]/2)
    '<-'(
      tmpfunc,
      function(x)
      {
        tmp <- ss[x,]=='-'
        ans <- length(which(tmp))<b
        return(ans)
      }
    )
    ss <- ss[which(sapply(1:dim(ss)[1],tmpfunc)),]
    return(ss)
    '<-'(
      judgeIs,
      function(h)
      {
        return(length(which(as.vector(table(h))>1)) > 1)
      }
    )
    ss <- ss[unlist(apply(ss,1,judgeIs)),]
    
    tree <- MP_build(ss, print)
    tree <- tree + 1
    return(tree)
  }
)

'<-'(
  MP2,
  function(ss,no_space=FALSE)
  {
    if(no_space){
      ss <- ss[which(sapply(1:dim(ss)[1],function(x)(return(all(ss[x,]!='-'))))),]
    }
    '<-'(
      judgeIs,
      function(h)
      {
        return(length(which(as.vector(table(h))>1)) > 1)
      }
    )
    ss <- ss[unlist(apply(ss,1,judgeIs)),]

    tree <- MP_build2(ss, print)
    tree <- tree + 1
    return(tree)
  }
)

'<-'(
  getLei,
  function(hou_zhui)
  {
    dezero <- max(qian_zhui_he(hou_zhui))
    if(dezero==0){
      return(list())
    }
    hou_zhui <- paste(hou_zhui, collapse = '_')
    tmp <- character(3)
    tmp[1] = '[_]*'
    tmp[2] = paste(rep(0, dezero),collapse = '_')
    tmp[3] = '[_]*'
    pat <- paste(tmp, collapse='')
    '<-'(
      getPart,
      function(part)
      {
        tmp <- as.numeric(strsplit(part,'_')[[1]])
        if(dezero>1){
          tmp1 <- getLei(tmp)
        }
        tmp <- sort(tmp[which(tmp>0)])
        if(dezero>1){
          tmp1[[length(tmp1)+1]] <- tmp
          tmp <- tmp1
        }
        return(tmp)
      }
    )
    leis <- strsplit(hou_zhui, pat)[[1]]
    leis <- lapply(leis, getPart)
    return(leis)
  }
)

'<-'(
  getLeis,
  function(tr)
  {
    leis <- getLei(tr)
    ans <- list()
    cnt <- 1
    for(each in leis){
      len <- length(each)
      ans[cnt:(cnt+len-1)] <- each
      cnt <- cnt + len
    }
    return(ans)
  }
)

'<-'(
  addList,
  function(li1, li2)
  {
    len1 <- length(li1)
    len2 <- length(li2)
    li1[len1:(len1+len2-1)] <- li2
    return(li1)
  }
)

'<-'(
  getLeis2,
  function(tr)
  {
    nnum <- max(tr)
    leis <- list()
    for(k in 1:nnum)
    {
      leis[k] <- k
    }
    cnt <- nnum + 1
    delTree <- paste(tr, collapse = '_')
    rg <- regexpr('[0-9]*_[0-9]*_0', delTree)
    fp <- rg[1]
    fl <- attr(rg,'match.length')
    while(fp>0){
      this_node <- as.numeric(strsplit(substring(delTree,fp,fp+fl),'_')[[1]])
      
      leis[[cnt]] <- sort(c(leis[[this_node[1]]],leis[[this_node[2]]]))
      
      if(fp == 1){
        delTree <- paste(c(cnt, substring(delTree,fp+fl)),collapse = '')
      } else {
        delTree <- paste(c(substring(delTree,1,fp-1), cnt, substring(delTree,fp+fl)),collapse = '')
      }
      cnt <-cnt + 1
      rg <- regexpr('[0-9]*_[0-9]*_0', delTree)
      fp <- rg[1]
      fl <- attr(rg,'match.length')
    }
    return(leis)
  }
)

'<-'(
  getTree,
  function(tr)
  {
    nnum <- max(tr)
    cnt <- nnum + 1
    edges <- NULL
    delTree <- paste(tr, collapse = '_')
    rg <- regexpr('[0-9]*_[0-9]*_0', delTree)
    fp <- rg[1]
    fl <- attr(rg,'match.length')
    while(fp>0){
      this_node <- as.numeric(strsplit(substring(delTree,fp,fp+fl),'_')[[1]])
      tmp <- matrix(0,2,2)
      tmp[1,2] <- this_node[1]
      tmp[2,2] <- this_node[2]
      tmp[,1] <- cnt
      edges <- rbind(tmp, edges)
      
      if(fp == 1){
        delTree <- paste(c(cnt, substring(delTree,fp+fl)),collapse = '')
      } else {
        delTree <- paste(c(substring(delTree,1,fp-1), cnt, substring(delTree,fp+fl)),collapse = '')
      }
      cnt <-cnt + 1
      rg <- regexpr('[0-9]*_[0-9]*_0', delTree)
      fp <- rg[1]
      fl <- attr(rg,'match.length')
    }
    return(edges) 
  }
)

'<-'(
  getEdgeColor,
  function(tr,real_class)
  {
    edges <- getTree(tr)
    bd <- max(edges)
    tc <- table(real_class)
    ntc <- names(tc)
    sumColor <- matrix(0,dim(edges)[1],length(tc))
    '<-'(
      find_up,
      function(lower, ad)
      {
        sumColor[which(edges[,2] == lower),ad] <<- sumColor[which(edges[,2] == lower),ad] + 1
        if(lower == bd){
          return()
        }
        sapply(edges[which(edges[,2] == lower),1], find_up, ad)
      }
    )
    for(i in 1:length(real_class)){
      find_up(i, which(ntc==real_class[i]))
    }
    #print(sumColor)
    ans <- character(dim(edges)[1])
    for(i in 1:dim(edges)[1]){
      ans[i] <- ntc[which.max(sumColor[i,])]
    }
    
    return(ans) 
  }
)

'<-'(
  drawTree,
  function(tree, bts, nms, real_class=NULL)
  {
    edges <- getTree(tree)
    edges <- data.frame(from=edges[,1],to=edges[,2])
    if(!is.null(real_class))
    {
      edges$类群 <- real_class
    } else {
      edges$类群 = '默认类群'
    }
    lb <- length(bts)
    bts <- as.character(bts)
    bts[which(bts=='-1')] <- ''
    bts[1:length(nms)] <- nms
    vertices <- data.frame(name =as.character(1:lb),support = bts)
    ggraph(graph_from_data_frame(edges, vertices = vertices), layout = 'dendrogram', circular = T) + 
      #geom_edge_link() +
      #geom_edge_elbow2(aes(colour = 类群)) +
      geom_edge_elbow(aes(colour = 类群)) +
      coord_fixed() +
      #geom_edge_bend2(aes(colour = 类群)) +
      #geom_edge_diagonal0(aes(colour = 类群)) +
      geom_node_point() +
      geom_node_text(aes(label = support,), repel = TRUE) +
      #ggtitle('最大简约法构建进化树 \n 代码实现：MicosLiang \n 数据来源：童宗中等（2002）') + 
      theme_void()
  }
)

'<-'(
  genSamples,
  function(ss, num=100, no_space=TRUE)
  {
    suc <- require(parallel)
    if(no_space){
      ss <- ss[which(sapply(1:dim(ss)[1],function(x)(return(all(ss[x,]!='-'))))),] 
    }
    '<-'(
      judgeIs,
      function(h)
      {
        return(length(which(as.vector(table(h))>1)) > 1)
      }
    )
    ss <- ss[unlist(apply(ss,1,judgeIs)),]
    di <- dim(ss)
    get <- round(di[1]/2)
    if(suc && num > 4)
    {
      '<-'(
        MP_test,
        function(s)
        {
          Rcpp::sourceCpp('./RStudy/ML-buildTree.cpp', rebuild = F)
          return(MP_build2(s, print))
        }
      )
      clist <- parallel::makeCluster(4)
      sss <- parLapply(clist, 1:num, function(x){return(ss[sample(1:di[1],get,replace = T),])})
      sapply(ls(".GlobalEnv"), (function(t){return(clusterExport(clist, t, envir = environment()))}))
      try(
        trees <- parLapply(clist, sss, MP_test)
      )
      trees <- parLapply(clist, trees, function(x){return(x+1)})
      try(
        anss <- parLapply(clist, trees, getLeis2)  
      )
      stopCluster(clist)
    }
    else
    {
      sss <- lapply(1:num, function(x){return(ss[sample(1:di[1],get,replace = T),])})
      trees <- lapply(sss, MP_build2, print)
      trees <- lapply(trees, function(x){return(x+1)})
      anss <- lapply(trees, getLeis2) 
    }
    return(anss)
  }
)

'<-'(
  compare_Lei,
  function(li1, li2)
  {
    len1 <- length(li1)
    len2 <- length(li2)
    if(len1 != len2)
    {
      return(F)  
    }
    return(all(li1==li2))
  }
)

'<-'(
  bootStrap,
  function(ss, tree, num=100, sams=NULL, no_space=TRUE)
  {
    #tree[which(tree==0)] = (1:length(which(tree==0))) + max(tree)
    delTree <- paste(tree, collapse = '_')
    mainLeis <- getLeis2(tree)
    nnum <- length(mainLeis)
    ans <- rep(0,nnum)
    if(is.null(sams))
    {
      sams <- genSamples(ss, num, no_space) 
    }
    for(k in 1:num)
    {
      jiao <- intersect(sams[[k]], mainLeis)
      len <- length(jiao)
      if(len){
        for(j in 1:len)
        {
          for(i in 1:nnum)
          {
            if(compare_Lei(jiao[[j]], mainLeis[[i]]))
            {
              ans[i] <- ans[i] + 1
            }
          }
        }
      }
    }
    return(ans)
  }
)

'<-'(
  MP_buildTrees,
  function(ss,num=100,seqs=NULL,use=0,no_space=TRUE,real_class=NULL,have_out=F)
  {
    if(!(require(ggraph) && require(igraph)))
    {
      stop('Liang@You must install ggraph and their requirments!')
    }

    if(!is.null(seqs))
    {
      dis <- getDistance(seqs,use=use)
      allDis <- rowSums(dis)
      new_order <- order(allDis,decreasing = TRUE)
      if(!is.null(real_class))
      {
        real_class <- real_class[new_order]
      }
      ss <- ss[,new_order]
    }
    nms <- names(ss[1,])
    tree <- MP2(ss,no_space)
    bts <- bootStrap(ss,tree,num,no_space=no_space)
    tr_ls <- getLeis2(tree)
    if(have_out)
    {
      for(k in 1:length(tr_ls)){
        if(any(tr_ls[[k]]==which(new_order==length(seqs))))
        {
          bts[[k]] <- -1
        }
      } 
    }
    if(!is.null(real_class))
    {
      real_class <- getEdgeColor(tree,real_class)
    }
    drawTree(tree, bts, nms,real_class)
  }
)

'<-'(
  del_alig,
  function(seqs, bs = 0){
    di <- dim(seqs)
    yuzhi <- round(di[2]*bs)
    
    '<-'(
      tmpfunc,
      function(x)
      {
        tmp <- seqs[x,]=='-'
        ans <- length(which(tmp))<=yuzhi
        return(ans)
      }
    )
    seqs <- seqs[which(sapply(1:dim(seqs)[1],tmpfunc)),]
    
    return(seqs)
  }
)

'<-'(
  getAligList,
  function(seqs)
  {
    seqs <- lapply(seqs, as.matrix)
    len <- length(seqs)
    ans <- list()
    for(i in 1:len){
      ans[[i]] <- list()
      for(k in i:len){
        if(i==k){
          next
        }
        ss <- liang_NW(seqs[[i]],seqs[[k]],print)
        ss <- tolower(ss)
        ans[[i]][[k-i]] <- ss
      }
    }
    return(ans)
  }
)

'<-'(
  getDistance,
  function(seqs,alig_list=NULL,use=0,otu=FALSE)
  {
    '<-'(
      JC,
      function(ss)
      {
        ss <- del_alig(ss,0)
        tmp <- 1 - length(which(ss[,1] != ss[,2])) / dim(ss)[1]
        if(tmp > 0){
          tmp <- -0.75 * log(tmp)
          tmp <- ifelse(tmp<0,-1/tmp,tmp)
          return(tmp)
        }
        warning('too far')
        return(Inf)
      }
    )
    '<-'(
      Kimura,
      function(ss)
      {
        ss <- del_alig(ss,0)
        d <- dim(ss)[1]
        #转换 a->g || c->t
        d1 <- length(which((ss[,1]=='a' & ss[,2]=='g') | (ss[,1]=='g' & ss[,2]=='a') | (ss[,1]=='c' & ss[,2]=='t') | (ss[,1]=='t' & ss[,2]=='c')))
        #颠换
        d2 <- length(which(ss[,1] != ss[,2])) - d1
        d1 <- d1 / d
        d2 <- d2 / d
        a <- 1 / (1- 2*d1 - d2)
        b <- 1 / (1 - 2*d2)
        if(a > 0 && b > 0)
        {
          tmp <- 0.5*log(a) + 0.25*log(b)
          tmp <- ifelse(tmp<0,-1/tmp,tmp)
          return(tmp)
        }
        warning('too far')
        return(Inf)
      }
    )
    len <- length(seqs)
    dis <- matrix(0,len,len,dimnames = list(names(seqs),names(seqs)))
    if(is.null(alig_list))
    {
      alig_list <- getAligList(seqs) 
    }
    for(i in 1:len){
      for(k in i:len){
        if(i==k){
          next
        }
        ss <- alig_list[[i]][[k-i]]
        if(use==0){
          dis[k,i] <- JC(ss)
        }
        else {
          dis[k,i] <- Kimura(ss)
        }
      }
    }
    return(as.matrix(as.dist(dis)))
  }
)


'<-'(
  NJ,
  function(distance)
  {
    tree <- NJ_build(distance, print) + 1
    return(tree)
  }
)

'<-'(
  NJ_bootstrap,
  function(tree, alig_list,seqs,num=500)
  {
    delTree <- paste(tree, collapse = '_')
    mainLeis <- getLeis2(tree)
    nnum <- length(mainLeis)
    ans <- rep(0,nnum)
    
    len <- length(alig_list)
    '<-'(
      getSmps,
      function(x)
      {
        tans <- list()
        for(i in 1:len){
          tans[[i]] <- list()
          for(k in i:len){
            if(i==k){
              next
            }
            ss <- alig_list[[i]][[k-i]]
            tlen <- dim(ss)[1]
            ss <- ss[sample(1:tlen,round(tlen/2),replace = TRUE),]
            tans[[i]][[k-i]] <- ss
          }
        }
        return(tans)
      }
    )
    smps <- lapply(1:num,getSmps)
    for(i in 1:num)
    {
      tr <- NJ(getDistance(seqs,smps[[i]]))
      jiao <- intersect(getLeis2(tr), mainLeis)
      tlen <- length(jiao)
      if(tlen){
        for(j in 1:tlen)
        {
          for(k in 1:nnum)
          {
            if(compare_Lei(jiao[[j]], mainLeis[[k]]))
            {
              ans[k] <- ans[k] + 1
            }
          }
        }
      }
    }
    return(ans)
  }
)

'<-'(
  NJ_buildTrees,
  function(seqs,otu=F,use=0,num=100,real_class=NULL)
  {
    if(!(require(ggraph) && require(igraph)))
    {
      stop('Liang@You must install ggraph and their requirments!')
    }
    nms <- names(seqs)
    alig_list <- getAligList(seqs)
    distance <- getDistance(seqs,alig_list,use,otu)
    tree <- NJ(distance)
    bts <- NJ_bootstrap(tree, alig_list, seqs,num)
    tr_ls <- getLeis2(tree)
    if(otu)
    {
      for(k in 1:length(tr_ls)){
        if(any(tr_ls[[k]]==length(seqs)))
        {
          bts[[k]] <- -1
        }
      } 
    }
    if(!is.null(real_class))
    {
      real_class <- getEdgeColor(tree,real_class)
    }
    drawTree(tree, bts, nms, real_class)
  }
)
