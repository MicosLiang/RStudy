'<-'(
  sigmoid,
  function(x)
  {
    return(1/(1+exp(-1*x)))
  }
)

'<-'(
  sigmoid_loss,
  function(x)
  {
    tmp <- sigmoid(x)
    return(tmp*(1-tmp))
  }
)

'<-'(
  turn_back,
  function(x)
  {
    return(x)
  }
)

'<-'(
  turn_back_loss,
  function(x)
  {
    return(1)
  }
)

'<-'(
  new.layer,
  function(info, input_num, neuron_num,active_function,loss_function,eta=0.01)
  {
    layer <- info
    attr(layer, 'connection') <- matrix(rnorm(input_num*neuron_num),input_num,neuron_num)
    attr(layer, 'threshold') <- rnorm(neuron_num)
    attr(layer, 'active') <- active_function
    attr(layer, 'loss') <- loss_function
    attr(layer, 'eta') <- eta
    '<-'(
      attr(layer, 'predict'),
      function(layer, x)
      {
        return(attr(layer, 'active')((x %*% attr(layer, 'connection')) - attr(layer, 'threshold')))
      }
    )
    '<-'(
      attr(layer, 'fix'),
      function(layer, input, output, error)
      {
        dt <- attr(layer, 'eta')*attr(layer, 'loss')(output)*error
        attr(layer, 'threshold') <- attr(layer, 'threshold') + dt
        dw <- input %*% t(dt)
        attr(layer, 'connection') <- attr(layer, 'connection') + dw
        return(layer)
      }
    )
    return(layer)
  }
)

'<-'(
  new.network,
  function(info, input_dim, output_dim, layers_dim, is_div = T)
  {
    network <- info
    attr(network, 'deep') <- length(layers_dim)
    #attr(network, 'errors') <- 'error'
    #attr(attr(network, 'errors'), 'out') <- rep(0, output_dim)
    in_num <- input_dim
    cnt <- 1
    for(i in layers_dim)
    {
      attr(network, as.character(cnt)) <- new.layer(cnt, in_num, i, sigmoid, sigmoid_loss, 0.01)
      #attr(attr(network, 'errors'), as.character(cnt)) <- rep(0, i)
      cnt <- cnt + 1
      in_num <- i
    }
    if(is_div)
    {
      attr(network, 'out')  <- new.layer('out', in_num, output_dim, sigmoid, sigmoid_loss, 0.01)
    }
    else
    {
      attr(network, 'out')  <- new.layer('out', in_num, output_dim, turn_back, turn_back_loss, 0.01) 
    }
    '<-'(
      attr(network, 'predict'),
      function(network, x)
      {
        ps <- 'errors'
        attr(ps, '0') <- x
        for(i in as.character(1:attr(network, 'deep')))
        {
          layer <- attr(network, i)
          x <- attr(layer, 'active')((x %*% attr(layer, 'connection')) - attr(layer, 'threshold'))
          attr(ps, i) <- x
        }
        layer <- attr(network, 'out')
        tmp <- ps
        ps <- attr(layer, 'active')((x %*% attr(layer, 'connection')) - attr(layer, 'threshold'))
        for(i in as.character(0:attr(network, 'deep')))
        {
          attr(ps, i) <- attr(tmp, i)
        }
        return(ps)
      }
    )
    '<-'(
      attr(network, 'fix'),
      function(network, xs, ys, iter=100)
      {
        #pre-training
        x_num <- dim(xs)[1]
        deep <- attr(network, 'deep')
        x_dim <- dim(xs)[2]
        for(each in as.character(1:deep))
        {
          for(it in 1:iter)
          {
            for(k in 1:x_num)
            {
              x <- xs[k,]
              y <- ys[k,]
              layer <- attr(network, each)
              vna_layer <- new.layer('vna',dim(attr(layer,'connection'))[2],x_dim,turn_back,turn_back_loss)
              layer <- attr(network, each)
              ps <- attr(network, 'predict')(network, x)
              tx <- attr(ps, each)
              error <- x - attr(vna_layer, 'predict')(vna_layer, tx)
              if(each != '1')
              {
                ox <- attr(ps, as.character(as.numeric(each)-1))
              }
              else
              {
                ox <- x
              }
              attr(network, each) <- attr(layer, 'fix')(layer, ox, tx, error)
              x_dim <- length(tx) 
            }
          }
        }
        attr(network, 'out') <- vna_layer
        return(network)
        
        hsy <- numeric(iter)
        #fine-tuning
        for(it in 1:iter)
        {
          for(k in 1:x_num)
          {
            x <- xs[k,]
            y <- ys[k,]
            yp <- attr(network, 'predict')(network, x)
            error <- y - yp
            hsy[it] <- sum(error^2)
            layer <- attr(network, 'out')
            connect_weights <- attr(layer, 'connection')
            attr(network, 'out') <- attr(layer, 'fix')(layer, attr(yp, as.character(deep)), yp, error)
            error <- attr(layer, 'loss')(yp) * error
            for(i in as.character(deep:1))
            {
              layer <- attr(network, i)
              in_x <- attr(yp, as.character(as.numeric(i)-1))
              error <- as.vector(connect_weights %*% error)
              connect_weights <- attr(layer, 'connection')
              attr(network, i) <- attr(layer, 'fix')(layer, in_x, attr(yp, i), error)
              error <- attr(layer, 'loss')(yp) * error
            }
          }
        }
        plot(1:iter, hsy,'l')
        return(network)
      }
    )
    return(network)
  }
)

dat1 <- c(10.9,14.1,12.6,13.2,11.6)
dat2 <- c(6.1,4.5,5.5,6.0,4.7)
plot(dat1,rep(0,5),col=2,cex=3,pch=16,xlim = c(0,15),ylim=c(-5,5),xlab='',ylab='length')
lines(dat2, rep(0,5),col=3,cex=3,pch=16,type='p')
legend(legend = c('月季','玫瑰'),'topleft',pch=c(16,16),col=c(2,3))