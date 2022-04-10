#Author : liang
#info   : partices

'<-'(
  matrixAddup,
  function(m)
  {
    temp.c <- c()
    apply(m,2,function(mr) temp.c <<- c(temp.c, mr))
    return(temp.c)
  }
)

main.tx1 <- getFitData(iris3[,,1])
main.tx2 <- getFitData(iris3[,,2])
main.tx3 <- getFitData(iris3[,,3])

main.tt1 <- iris3[,,1][main.tx1,]
main.tt2 <- iris3[,,2][main.tx2,]
main.tt3 <- iris3[,,3][main.tx3,]

main.td1 <- iris3[,,1][main.tx1*-1,]
main.td2 <- iris3[,,2][main.tx2*-1,]
main.td3 <- iris3[,,3][main.tx3*-1,]

main.testdata <- array(c(matrixAddup(main.tt1),matrixAddup(main.tt2),matrixAddup(main.tt3)), dim=c(10,4,3))
main.traindata <- array(c(matrixAddup(main.td1),matrixAddup(main.td2),matrixAddup(main.td3)), dim=c(40,4,3))

main.wpn <- rep(10,3)

'<-'(main.data, main.traindata)
'<-'(main.margin, dim(main.data))
'<-'(main.weights, t(matrix(rnorm(1 + main.margin[2], 0, 0.01), main.margin[3], 1 + main.margin[2], byrow = TRUE)))
'<-'(main.t_weights, matrix(0, main.margin[3], 1 + main.margin[1]))
'<-'(main.type, matrix(c(1, rep(-1, main.margin[3])), main.margin[3], main.margin[3]))
'<-'(main.eta, 0.001)
'<-'(main.n_iter, 2000)
'<-'(main.wrong, rep(0, main.margin[3]))
'<-'(main.record, matrix(0, main.n_iter, main.margin[3]))

'<-'(getFitData,
    function(population, number = 0.2, byrow = TRUE)
    {
      if(byrow)
      {
        temp.dim = dim(population)[1]
        '<-'(temp.c,
             sample
             (
                1:temp.dim,
                floor
                (
                  '*'(temp.dim, number)
                ),
                replace = FALSE
             )
        )
        #return(list(population[temp.c,], population[temp.c*-1,]))
        return(temp.c)
      }
    }
)

'<-'(
  thre,
  function(z)
  {
    return((z<=0)*-1 + (z>0)*1)
  }
)

'<-'(
  getDw,
  function(k, xis,way)
  {
    #print(main.t_weights[k,] + main.eta * out * xis)
    #print(main.eta * out * xis)
    if(way[k]!=0)
    {
      main.wrong[k] <<- main.wrong[k] + 1
    }
    "<<-"(main.t_weights[k,], main.t_weights[k,] + main.eta * way[k] * xis)
    #print(main.t_weights[k,])
  }
)

'<-'(
  learn,
  function(i, way, data)
  {
    lapply(1:length(way[i,]),getDw,xis = data[i,],way=way[i,])
  }
)

'<-'(
  getZ,
  function(x_)
  {
    #print(main.t_weights)
    '<-'(temp.data, cbind(1, t(x_)))
    temp.z <- temp.data %*% main.weights
    temp.thr <- thre(temp.z)
    temp.way <- main.type - temp.thr
    #print(temp.way)
    lapply(1:dim(temp.way)[1],learn,way=temp.way,data=temp.data)
  }
)

'<-'(
  fit,
  function(v){
    '<<-'(main.t_weights, matrix(0, main.margin[3], 1 + main.margin[2]))
    apply(main.data,1,getZ)
    main.weights <<- t(main.t_weights) + main.weights
    '<<-'(main.record[v,], main.wrong)
    main.wrong <<- rep(0, main.margin[3])
  }
)

'<-'(
  guess,
  function(t,j)
  {
    temp.cnt <- 0
    temp.func <- function(dat)
    {
      dat <- c(1,dat)
      if((dat %*% main.weights[,j])[1,1] > 0)
      {
        temp.cnt <<- temp.cnt + 1
      }
    }
    apply(main.testdata[,,t],1,temp.func)
    print(temp.cnt)
  }
)

'<-'(
  cacul,
  function(t1,i)
  {
    t1 <- c(1,t1)
    temp.ans <- t1 %*% main.weights
    if(temp.ans[1] > 0){
      temp.type <- 1
    }
    else if(temp.ans[3] > 0)
    {
      temp.type <- 3
    }
    else
    {
      temp.type <- 2
    }
    if(temp.type != i)
    {
      main.wpn[i] <<- main.wpn[i] - 1
    }
  }
)

'<-'(
  test,
  function()
  {
    lapply(1:main.margin[3],function(i) apply(main.testdata[,,i],1,cacul,i=i))
  }
)

'<-'(
  drawLine,
  function(dat){
    plot(x=1:length(dat),y=dat,'o')
  }
)
  
'<-'(
  draw,
  function(dat)
  {
    apply(dat, 2, drawLine)
  }
)

#train
lapply(1:main.n_iter, fit)

#
dev.new()
plot(x=1:length(main.record[,1]),y=main.record[,1],'l',col='DarkTurquoise',xlab="n of iter",ylab='update times',main="parameters's n-upt")
points(x=1:length(main.record[,2]),y=main.record[,2],'l',col='RosyBrown')
points(x=1:length(main.record[,3]),y=main.record[,3],'l',col='DeepPink')
lines(x=1:length(main.record[,2]),y=main.record[,2],col='RosyBrown')
lines(x=1:length(main.record[,3]),y=main.record[,3],col='DeepPink')
legend(1240,80,c("Setosa","Versicolor","Virginica"),col=c('DarkTurquoise','RosyBrown','DeepPink'),text.col=c('DarkTurquoise','RosyBrown','DeepPink'))

#
dev.new()
test()
drawIris<-function(line, colr)
{
  points(x=c(1,2,3,4),y=line,col=colr)
  lines(x=1:4,y=line,col=colr)
}
plot(x=rep(1,4),y=rep(-10,4),xlim=c(1,4),ylim=c(0,10),xlab="features",ylab="value",main="iris3's feature space")
apply(iris3[,,1],1,drawIris,colr='DarkTurquoise')
apply(iris3[,,2],1,drawIris,colr='RosyBrown')
apply(iris3[,,3],1,drawIris,colr='DeepPink')
legend(3.2,9.6,c("Setosa","Versicolor","Virginica"),col=c('DarkTurquoise','RosyBrown','DeepPink'),text.col=c('DarkTurquoise','RosyBrown','DeepPink'))

#
dev.new()
barplot(main.wpn,names.arg = c("Setosa","Versicolor","Virginica"),col=c('DarkTurquoise','RosyBrown','DeepPink'),main="Correct classification rate")