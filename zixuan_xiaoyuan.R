jianqiuluo <- as.matrix(read.csv('./RStudy/data/school_plants/jianqiuluo.csv'))
w1 <- my_lm2(jianqiuluo[,2:4],jianqiuluo[,1])

lm_R(jianqiuluo[,2:4],jianqiuluo[,1],w1)
tmp_t(jianqiuluo[,2:4],jianqiuluo[,1],w1)

conghua <- as.matrix(read.csv('./RStudy/data/school_plants/songhua.csv'))
w2 <- my_lm2(conghua[,2:4],conghua[,1])
lm_R(conghua[,2:4],conghua[,1],w2)
tmp_t(conghua[,2:4],conghua[,1],w2)

require(gt)
my_iris <- as.matrix(read.csv('./RStudy/data/school_plants/iris.csv'))
'<-'(
  group_difference,
  function(g1, g2)
  {
    g1_mean <- mean(g1)
    g2_maen <- mean(g2)
    ga  <- c(g1,g2)
    all_mean <- mean(ga)
    msa <- ((g1_mean-all_mean)^2 * length(g1) + (g2_maen-all_mean)^2 * length(g2)) / (length(g1)+length(g2))
    mse <- sum(((ga - all_mean)^2)) / length(ga)
    f <- msa / mse
    return(df(f, length(g1), length(g2)))
  }
)
fenxi <- function()
{
  p1 <- matrix(NA,6,4,dimnames = list(NULL,c('row','num','mean','var')))
  for(i in 1:6)
  {
    p1[i,] <- c(names(my_iris[1,])[ifelse(i%%3==0,1,i%%3)],length(my_iris[,i]),mean(my_iris[,i]),var(my_iris[,i]))
  }
  p1 <- as.data.frame(p1)
  p1$group <- c(rep('leaf_length',3),rep('leaf_wide',3))
  gt(p1,rowname_col = 'row',groupname_col = 'group') %>%
  tab_header(title = 'Three type iris difference summary')
  p2 <- matrix('-',6,3,dimnames = list(NULL,c('purple','red','yellow')))
  for(i in 1:3)
  {
    for(k in 1:3)
    {
      if(i >= k)
      {
        next
      }
      p2[i,k] <- group_difference(my_iris[,i],my_iris[,k])
      p2[i,k] <- ifelse(as.numeric(p2[i,k])<0.05,paste(p2[i,k],'*'),p2[i,k])
    }
  }
  for(i in 4:6)
  {
    for(k in 4:6)
    {
      if(i >= k)
      {
        next
      }
      p2[i,k-3] <- group_difference(my_iris[,i],my_iris[,k])
      p2[i,k-3] <- ifelse(p2[i,k-3]<0.05,paste(p2[i,k-3],'*'),p2[i,k-3])
    }
  }
  p2 <- as.data.frame(p2)
  p2$group <- c(rep('leaf_length',3),rep('leaf_wide',3))
  p2$row <- names(my_iris[1,])[1:3]
  gt(p2,rowname_col = 'row',groupname_col = 'group') %>%
  tab_header(title = 'Three type iris difference p-value')
  
}
fenxi()
