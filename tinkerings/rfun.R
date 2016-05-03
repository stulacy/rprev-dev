.First=function(){
   require("ctest", quietly=TRUE)
   require(survival)
cat("\n")
}  

"weisim"  =   function(n=100,p=1,b=1,shp=1,rancens=T,rate=1,trunc=F,tau=2)
{ 
  x=matrix(0,n,p)
  for(it in 1:p){
    x[,it]=rnorm(n)
  }
  bx=x%*%as.vector(b)
  t=(-log(runif(n)))^(1/shp)/exp(bx)
  ic=rep(1,n)
  if(rancens==T){
    c1=(-log(runif(n))/rate)
    i=c1<t
    t[i]=c1[i]
    ic[i]=0
  }
  if(trunc==T){
   i=t>tau
    t[i]=tau
    ic[i]=0
 }
  i=order(t)
  time=t[i]
  cens=ic[i]
  x=x[i,]
  cat(100*sum(1-cens)/n, "percent censored\n")
data.frame(time,cens,x)    
}

weisim2 <- function(n=100, beta=NULL, mu=NULL, sigma=NULL, pbinomial=NULL, scale=1, shape=1, rancens=T, rate=1, trunc=F, tau=2) { 
    
    if (! length(mu) == length(sigma))
        stop("Error: mu and sigma must have the same length.")
    
    if (! length(beta) == length(mu) + length(pbinomial))
        stop("Error: Number of betas must be equal to length of mu + pbinomial")
    
    # Always have intercept, which is just scale if no other covariates are supplied
    x = cbind(rep(1, n))
    # Add continuous variables
    if (!is.null(mu) & !is.null(sigma)) 
        x = cbind(x, mapply(function(m, s) rnorm(n, m, s), mu, sigma))
    # Add binary variables
    if (!is.null(pbinomial)) 
        x = cbind(x, sapply(pbinomial, function(p) rbinom(n, 1, p)))
    
    bx = x %*% c(scale, beta)
    
    t = (-log(runif(n))) ^ (1/shape) * exp(bx)  # This is rweibull
    ic = rep(1,n)
    if (rancens == T) {
        c1 = (-log(runif(n))/rate)  # This is rexp
        i = c1 < t
        t[i] = c1[i]
        ic[i] = 0
    }
    if (trunc == T) {
        i = t > tau
        t[i] = tau
        ic[i] = 0
    }
    i = order(t)
    time = round(t[i])
    cens = ic[i]
    x = x[i,-1]
    cat(100 * sum(1-cens) / n, "percent censored\n")
    data.frame(time=time, status=cens, x)    
}

"simexpmart"  =function(n = 100, lambda = 1, tau = 2, rancens = F,
                         gamma=1,ngrid = 100)
{
	time =  - log(runif(n))/lambda
	cens = time * 0 + 1
	i = time > tau
	time[i] = tau
	cens[i] = 0
        if(rancens==T){
          c=-log(runif(n))/gamma
         i = time > c
	time[i] = c[i]
	cens[i] = 0
        }
	tgrid = seq(0, tau, length = ngrid)
	base = tgrid * lambda
	mart = matrix(0, n, ngrid)
      
		for(it in 1:n) {
		mart[it,  ] =  - base
		j = sum(tgrid < time[it])
		if(j < ngrid) {
			mart[it, (j + 1):ngrid] = mart[it, j]
			if(cens[it] == 1) {
				mart[it, (j + 1):ngrid] = mart[it, (j + 1):
					ngrid] + 1
				}
		}
		
	}
    
list(mart=mart,tgrid=tgrid) 
} 
 "simexpmart.plot" =  function(mart)
{ tau=max(mart$tgrid)
 
  
  plot(c(0, 0), c(0, 0), pch = " ", xlim = c(0, tau),
       ylim=1.1*range(as.vector(mart$mart)), xlab = "Time", ylab = "M(t)")
  for(it in 1:length(mart$mart[,1])) {
        lines(mart$tgrid, mart$mart[it,  ])
      }
  cat("\n")
}
"simexpmart.cs"  =function(mart,cstime)
{ j=sum(mart$tgrid<=cstime)
  m=mart$mart[,j]
  hist(m)
  cat("At cross-section time ",cstime,"\n")
  cat("  Mean M(t)",round(mean(m),3),"  (se=",round((var(m)/length(m))^.5,3),")\n")
  cat("  Variance M(t)",round(var(m),3),"\n")
  
 cat("\n") 
}

"expmartvar"  =   function(cstime,lambda = 1,rancens = F,gamma=1)
{if(rancens==F){ev=1-exp(-lambda*cstime)}else {
  ev=lambda*(1-exp(-(lambda+gamma)*cstime))/(lambda+gamma)}
print(round(ev,3))
 cat("\n")
}

"stratplot" =
  function(datafrm,cols,strtcol){
    usedat<<-datafrm
  xmat<<-as.matrix(usedat[,cols])
    srt<<-usedat[,strtcol]
    ss = survfit(coxph(Surv(usedat$time, usedat$cens) ~ xmat+ strata(srt)))
    ns = length(unique(ss$strata))
    n0 = 0
    a =  log(- log(ss$surv))
    plot(c(0, 0), c(0, 0), pch = " ", xlab = "Time", ylab = "log A(t)", xlim = 
         range(ss$time), ylim = range(a))
    for(it in 1:ns) {
      nu = ss$strata[it]
      j = n0 + (1:nu)
      lines(ss$time[j], a[j])
      n0 = n0 + nu
    }
    cat("\n")
  }


"piplot"  =   function(datafrm,cols)
{ usedat<<-datafrm
  xmat<<-as.matrix(usedat[,cols])
    if(length(cols)==1) { xmat=matrix(xmat,ncol=length(cols))}
  cc=coxph(Surv(usedat$time,usedat$cens)~xmat)
  ss=survfit(cc)
    cat("Solid line is fitted, dotted observed\n")
   if(length(cols)>1){
  xbar=apply(xmat,2,"mean")
}else { xbar=mean(as.vector(xmat))}
  b=cc$coefficients
  bb=b%*%xbar
   bx=(xmat)%*%b
  j=order(bx)  
  n=length(bx)
  bx1=bx[j][round(n/3)]
   bx2=bx[j][round(2*n/3)]
 
    i1=bx<=bx1
    i2=(bx>bx1)&(bx<=bx2)
    i3=(bx>bx2)
    b1=mean(bx[i1])
    b2=mean(bx[i2])
    b3=mean(bx[i3])
    s=ss$surv
    s1=s^{exp(b1-bb)}
    s2=s^{exp(b2-bb)}
       s3=s^{exp(b3-bb)}
  plot(ss$time,s1,type="s",ylim=c(0,1),xlim=c(0,max(usedat$time)),xlab="Time",ylab="Survival")
    lines(ss$time,s2,type="s")
    lines(ss$time,s3,type="s")
    ss1=survfit(Surv(usedat$time[i1],usedat$cens[i1]))
    lines(ss1$time,ss1$surv,lty=2,type="s")
      ss1=survfit(Surv(usedat$time[i2],usedat$cens[i2]))
    lines(ss1$time,ss1$surv,lty=2,type="s")
        ss1=survfit(Surv(usedat$time[i3],usedat$cens[i3]))
    lines(ss1$time,ss1$surv,lty=2,type="s")
cat("\n")
}
"plaal" =   function(datafrm, cols, maxtime, mpl = c(2, 2))
{
  xmat = as.matrix(datafrm[, cols])
  tt = datafrm$time
  ic = datafrm$cens
  i=order(tt)
  tt=tt[i]
  ic=ic[i]
  xmat=xmat[i,]
  n = length(tt)
  x = cbind(rep(1, length = n), xmat)
  p = dim(x)[2]
  aa = matrix(0, n, p + 1)
  ta = rep(0, n)
  nc = 0
  nuse = sum(tt < maxtime)
  for(it in 1:nuse) {
    cat(it, " out of ", nuse, "\n")
    if(ic[it] == 1) {
      nc = nc + 1
      i = tt >= tt[it]
      xx = x[i,  ]
      dn = rep(0, length = sum(i))
      dn[1] = 1
      aa[nc, 2:(p + 1)] = solve(t(xx) %*% xx) %*% t(xx) %*%
        dn
      aa[nc, 1] = tt[it]
    }
  }
  aalhaz = aa[1:nc,  ]
  par(mfrow = mpl)
  for(it in 1:(p - 1)) {
    plot(aalhaz[, 1], cumsum(aalhaz[, it + 2]), type = "s", xlab = 
         "Time", ylab = "A(t)", cex = 0.59999999999999998)
  }
  par(mfrow=c(1,1))
  cat("\n")
}
"weisimcp" =function(n = 100, p = 1, b1 = 1,b2=1,cptau=1, shp = 1, rancens = T, rate = 1, trunc = F,
           maxtau = 2)
{
  x = matrix(0, n, p)
  for(it in 1:p) {
    x[, it] = rnorm(n)
  }
  bx1 = x %*% as.vector(b1)
  bx2 = x %*% as.vector(b2)
  t = ( - log(runif(n))/exp(bx1))^(1/shp)
  i=t>cptau
  t[i]=(cptau^shp - log(runif(sum(i)))/exp(bx2[i]))^(1/shp)

  
  ic = rep(1, n)
  if(rancens == T) {
    c1 = ( - log(runif(n))/rate)
    i = c1 < t
    t[i] = c1[i]
    ic[i] = 0
  }
  if(trunc == T) {
    i = t > maxtau
    t[i] = maxtau
    ic[i] = 0
  }
  i = order(t)
  time = t[i]
  cens = ic[i]
  x = x[i,  ]
  cat((100 * sum(1 - cens))/n, "percent censored\n")
  data.frame(time, cens, x)
}
"infplot"  =  function(datafrm,cols)
{   usedat <<- datafrm
    xmat <<- as.matrix(usedat[, cols])
  tt=usedat$time
  ic=usedat$cens
  x=xmat
  cc=coxph(Surv(tt,ic)~x)
  b=cc$coefficient
  v=solve(cc$var)
  n=length(tt)
  d=rep(0,length=length(tt))
  for(it in 1:n) {
  cat(it," out of ", n,"\n")  
  tt1=tt[-it]
  ic1=ic[-it]
  x1=x[-it,]
  b1=coxph(Surv(tt1,ic1)~x1)$coefficient
  d[it]=t(b1-b)%*%v%*%(b1-b)
}
  dinf<<-d
  plot(1:length(d),d,xlab="Subject",ylab="D",cex=0.7)
  cbind(1:n,d)
}

"infplot1" <-
  function(datafrm,cols)
{   usedat <<- datafrm
    xmat <<- as.matrix(usedat[, cols])
    tt<-usedat$time
    ic<-usedat$cens
    x<-xmat
    cc<-coxph(Surv(tt,ic)~x)
    dfbeta<-residuals(cc,type="dfbeta")
    if(length(cols)==1){dfbeta=matrix(dfbeta,ncol=1)}
    b<-cc$coefficient
    
    v<-solve(cc$var)
     n<-length(tt)
    d1<-rep(0,length=length(tt))
    for(it in 1:n) {
      cat(it," out of ", n,"\n")  
                  d1[it]<-t(dfbeta[it,])%*%v%*%dfbeta[it,]
    }
   
    plot(1:length(d1),d1,xlab="Subject",ylab="D",cex=0.7)
    cbind(1:n,d1)
  }



"coxphcp"  =   function(datafrm,cols,cptau)
{ usedat = datafrm
    x = as.matrix(usedat[, cols])
  tt=usedat$time
  ic=usedat$cens
  tt1=tt
  ic1=ic
  i=tt1>cptau
  tt1[i]=cptau
  ic1[i]=0
  tt2=tt
  ic2=ic
  i=tt2<=cptau
  ic2[i]=0
  cc=coxph(Surv(tt,ic)~x)
  l0=cc$loglik[2]
  cat("No change\n")
  print(cc)
  cat("\n-----------------------------------------------------\n")
   cc1=coxph(Surv(tt1,ic1)~x)
  l1=cc1$loglik[2]
  cat("Before change\n")
  print(cc1)
  cat("\n-----------------------------------------------------\n")
   cc2=coxph(Surv(tt2,ic2)~x)
  l2=cc2$loglik[2]
  cat("After change\n")
  print(cc2)
  cat("\n-----------------------------------------------------\n")
  cat("Log-likelihoods\n\nNo change\n")
  print(l0)
  cat("Before change\n")
  print(l1)
  cat("After change\n")
  print(l2)
  cat("Total with change\n")
  print(l1+l2)
  
list(all=cc,before=cc1,after=cc2)  
}

"piplotgamfrail" =  function(datafrm,cols)
{ usedat<<-datafrm
  xmat<<-as.matrix(usedat[,cols])
  n=length(usedat$time)
   if(length(cols)==1) { xmat=matrix(xmat,ncol=length(cols))}
  cc=coxph(Surv(usedat$time,usedat$cens)~xmat+frailty(1:length(usedat$time)))
  var=cc$history$"frailty(1:length(usedat$time)"$theta
  cat("Frailty variance\n")
  print(var)
    cat("Solid line is fitted, dotted observed\n")
  ss=survfit(cc)
  cumhaz=-log(ss$surv)
   if(length(cols)>1){
  xbar=apply(xmat,2,"mean")
}else { xbar=mean(as.vector(xmat))}
  b=cc$coefficients
  bb=b%*%xbar
  bx=(xmat)%*%b
  j=order(bx)  
  n=length(bx)
  bx1=bx[j][round(n/3)]
  bx2=bx[j][round(2*n/3)]
  
  i1=bx<=bx1
  i2=(bx>bx1)&(bx<=bx2)
  i3=(bx>bx2)
  b1=mean(bx[i1])
  b2=mean(bx[i2])
  b3=mean(bx[i3])
  s=ss$surv
  s1=(1/(1+var*exp(b1-bb)*cumhaz))^(1/var)
  s2=(1/(1+var*exp(b2-bb)*cumhaz))^(1/var)      
  s3=(1/(1+var*exp(b3-bb)*cumhaz))^(1/var)
  plot(ss$time,s1,type="s",ylim=c(0,1),xlim=c(0,max(usedat$time)),xlab="Time",ylab="Survival")
  lines(ss$time,s2,type="s")
  lines(ss$time,s3,type="s")
  ss1=survfit(Surv(usedat$time[i1],usedat$cens[i1]))
  lines(ss1$time,ss1$surv,lty=2,type="s")
  ss1=survfit(Surv(usedat$time[i2],usedat$cens[i2]))
  lines(ss1$time,ss1$surv,lty=2,type="s")
  ss1=survfit(Surv(usedat$time[i3],usedat$cens[i3]))
  lines(ss1$time,ss1$surv,lty=2,type="s")
  cat("\n")
}
"weisimfrail" =   function(n = 100, frvar=1,p = 1, b = 1, shp = 1, rancens = T, rate = 1, trunc = F,
           tau = 2)
{
  x = matrix(0, n, p)
  for(it in 1:p) {
    x[, it] = rnorm(n)
  }
  bx = x %*% as.vector(b)
   if(frvar>0){
  z=rgamma(n,1/frvar)*frvar} else {
      z=rep(1,length=n)}
  t = ( - log(runif(n))/(z*exp(bx)))^(1/shp)
  ic = rep(1, n)
  if(rancens == T) {
    c1 = ( - log(runif(n))/rate)
    i = c1 < t
    t[i] = c1[i]
    ic[i] = 0
  }
  if(trunc == T) {
    i = t > tau
    t[i] = tau
    ic[i] = 0
  }
  i = order(t)
  time = t[i]
  cens = ic[i]
  x = x[i,  ]
  cat((100 * sum(1 - cens))/n, "percent censored\n")
  data.frame(time, cens, x)
}

"weisimfrail2" =
  function(n = 100, frvar = 1, p = 1, b = 1, shp = 1, rancens = F, rate = 1, trunc = F, tau = 2,plt=T)
{
  x = matrix(0, n, p)
  for(it in 1:p) {
    x[, it] = rnorm(n)
  }
  bx = x %*% as.vector(b)
  if(frvar>0){
    z = rgamma(n, 1/frvar) * frvar}else{
      z=rep(1,length=n)}
  t1 = ( - log(runif(n))/(z * exp(bx)))^(1/shp)
  t2=( - log(runif(n))/(z * exp(bx)))^(1/shp)
  
  ic1 = rep(1, n)
  ic2=rep(1,n)
  if(rancens == T) {
    c1 = ( - log(runif(n))/rate)
    i = c1 < t1
    t1[i] = c1[i]
    ic1[i] = 0
    i=c1<t2
    t2[i] = c1[i]
    ic2[i] = 0
  }
  if(trunc == T) {
    i = t1 > tau
    t1[i] = tau
    ic1[i] = 0
    i = t2 > tau
    t2[i] = tau
    ic2[i] = 0
  }
  t=c(t1,t2)
  ic=c(ic1,ic2)
  x=rbind(x,x)
  sub=c(1:n,1:n)
  time=t
  cens = ic
  subject=sub
  cat((100 * sum(1 - cens))/(2*n), "percent censored\n")
  if((plt==T)&(rancens==F)&(trunc==F)){
    par(mfrow=c(1,2)) 
    plot(t1,t2,xlim=c(0,max(t1)),ylim=c(0,max(t2)),xlab="T1",ylab="T2")
    text(0.1*max(t1),0.9*max(t2),"Original")
    if(frvar>0){
      u1=1-(1+frvar*exp(bx)*t1^shp)^(-1/frvar)
      u2=1-(1+frvar*exp(bx)*t2^shp)^(-1/frvar)}
    else{ u1=1-exp(-exp(bx)*t1^shp)
          u2=1-exp(-exp(bx)*t2^shp)}
    plot(u1,u2,xlim=c(0,1),ylim=c(0,1),xlab="U1",ylab="U2")
    text(0.22,0.9,"Uniformly transformed")
    u=c(u1,u2)
    cat("Correlations",round(cor(t1,t2),2), " for T and ",round(cor(u1,u2),2)," for U\n")      } 
  data.frame(time, cens, x,subject)
}


"compden"  =   function(v1,v2)
{ par(mfrow=c(1,2))
  xm=qgamma(0.995,1/v1)*v1
  
  x=seq(0.001,xm,length=100)
  a=1/v1
  d1=a^a*x^(a-1)*exp(-a*x)/gamma(a)
  d2=exp(-log(x)^2/(2*v2))/(x*sqrt(2*pi*v2))
plot(x,d1,type="l",xlab="Original scale",ylab="Pdf",ylim=c(0,1.2*max(c(d1,d2))))
  
  lines(x,d2,lty=2)
  legend(0.5*xm,max(c(d1,d2)),c("Gamma","Lognormal"),lty=1:2,cex=0.7)
  xm=2*log(xm)
   x=seq(-xm,xm,length=100)
   d1=a^a*exp(x)^a*exp(-a*exp(x))/gamma(a)
   plot(x,d1,type="l",xlab="Log scale",ylab="Pdf")
  d2=exp(-x^2/(2*v2))/sqrt(2*pi*v2)
  lines(x,d2,lty=2)  
  cat("\n")
}







