'<-'(
  chi_square,
  function(xv, xe, typ=0)
  {
    if(typ==1)
    {
      xv <- xv / sum(xv)
    }
    else if(typ==2)
    {
      xe <- xe * sum(xv)
    }
    return(sum((xv-xe)^2/xe))
  }
)

'<-'(
  x2_live,
  function(x2,df)
  {
    k <- df/2
    return(0.5^k/my_gamma(k) * x2 ^ (k-1) * exp(1)^(-x2/2))
  }
)

'<-'(
  ingamma,
  function(z,s,isdown=F)
  {
    b <- 1/s*z^s*exp(1)^-z
    fs <- 1
    M <- 1
    for(i in 1:100)
    {
      fs <- fs*(s+i)
      M <- M + z^i/fs
    }
    if(isdown)
    {
      return(my_gamma(s)-b*M)
    }
    return(b*M)
  }
)

'<-'(
  my_gamma,
  function(z)
  {
    #it seems well when dof<100
    return(sqrt(2*pi/z)*((1/exp(1))*(z+1/((12*z)-(1/10*z))))^z)
  }
)

'<-'(
  chi_p,
  function(x,k)
  {
    return(ingamma(x/2,k/2,T)/my_gamma(x/2))
  }
)