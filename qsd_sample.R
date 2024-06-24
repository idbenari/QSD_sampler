# Simulation of QSD for a process 
# Offspring distribution parameters
library(compiler)
library(parallel)

# Geometric model
p_0 = 0.6 # Probability of dying
p =  0.3  # Geometric 
## Real is Geom p/p_0
q=1-p/p_0
m = (1-p_0)/(1-p) 
mi=1/m
## sampling from offspring, starting with population of size x
TV_calc = function(mu)
{
  mu[,2]=mu[,2]/sum(mu[,2])
  diff=0.5*sum(abs(dgeom(mu[,1]-1,q)-mu[,2]))+0.5*(1-sum(dgeom(mu[,1]-1,q)))
  return(diff)
}

# Polynomial model
## offspring distribution
OffSpring = function (n=20)
  # Polynomial in paper, r = 1.052907
#  return(c(0.838,0.008,0.031,0.011,0.021,0.029,0.019,0.014,0.029))
 return (c(p_0,(1-p_0)*(1-p)*p^(0:(n-2))))
## Set our distribution
mu=OffSpring()
## Degree of polynomial
L=length(mu)-1
## Expectation
m=sum(mu*(0:L))
mi=1/m

## Heavy tail 
m=(1-p_0)*pi^2/6
mi=1/m


# General sampler command
## cur = current state
## x = where to return if zero
psample = function (cur,x=1)
{
  if (cur==0) return(x)
## Geometric
  return(sum(rbinom(cur,1,1-p_0)*(rgeom(cur,1-p)+1)))
## Heavy tail
#  return(sum(rbinom(cur,1,1-p_0)*(floor(runif(cur)^(-0.5)))))## Use a fixed distribution mu 
## Fixed distribution 
#   return(sum(sample(0:L,cur,prob=mu,replace=TRUE)))
}
psample_cmp=cmpfun(psample)


# Our data collection and processing 
# data frame containing data collected

# Cycle path sampler one step.
CycSample = function(x=1,T=10^6) # x serves as an anchor point, T is number of cycles from x to 0 or x
{
path=c(x)
t=1
powers=c(0)
cycles=0
while(cycles<T)
{ 
# Check if cycle is complete
  if (path[t] %in% c(0,x) & powers[t]>0 ) 
  {
   # mark completion 
   cycles=cycles+1
   # record the time
#   powers[t]=-powers[t]
   # start new
   path[t+1]=x
   powers[t+1]=0
  }
  else 
  {  nxt=psample_cmp(path[t],x)
     path[t+1]=nxt
     powers[t+1]=powers[t]+1
  }
t=t+1
}
return(data.frame(state=path[1:(t-1)],weight=powers[1:(t-1)]))
##summary=aggregate(by=list("state"=path),x=list("weight"=powers),FUN=c)
}
CycSample_cmp = cmpfun(CycSample)

# Return map approach (when absorbed jump back according to empirical distribution)
RetSample = function(x=3,T=10^6)
{
  path=c()
  path[1]=x
  t=1
  while (t<=T)
  {
    cur=path[t]
    if (cur==0) 
    {
      nozero=which(path!=0)
      path[t+1]=sample(path[nozero],1)
    }
    else
    {
    path[t+1]=psample_cmp(cur)
    }
    t=t+1
  }
  nozero=which(path!=0)
  summary = aggregate(by = list("state"=path[nozero]), x = list("weight"=rep(1,length(nozero))), FUN = sum)
  return(summary)
}
RetSample_cmp=cmpfun(RetSample)

# Running the simulations
# and processing the data

seed=10

## Multicore 
set.seed(seed)
T=10000
cores=12
x=1
ctime= proc.time()
# Getting our samples
cQSDtemp=mclapply(1:cores,function (t) CycSample_cmp(x,T))
# Collecting all states visited
state=unlist(sapply(1:cores, function(n) cQSDtemp[[n]]$state))
# Collecting all weights
weight=unlist(sapply(1:cores, function(n) cQSDtemp[[n]]$weight))
cQSD=aggregate(x=list("weight"=weight),by=list("state"=state),FUN=table)
rm(cQSDtemp)

nst=length(cQSD[,1])
# cQSD=cQSD[2:nst,]
nst=nst-1
for (t in 1:nst) 
     cQSD[t,3]=sum(mi^as.integer(names(cQSD[t,2][[1]]))*cQSD[t,2][[1]])

#cQSD[,3]=cQSD[,3]/sum(cQSD[,3])
#cQSD=cQSD[,c(1,3)]
ctime=proc.time()-ctime
ctime


set.seed(seed)
rtime= proc.time()
  rQSD=RetSample_cmp(T=5*10^5)
rtime=proc.time()-rtime
rtime

# Numerical method
P = function (z)
{ 
  # Geom with parameter 1-p
  #p_0+(1-p_0)*(1-p)*z /(1-z*p)
  return(sum(mu * z^(0:L)))
}
P_cmp=cmpfun(P)


# Returns vector with numerical approximation of probabilities
# for 1,...,n 
# Must have:
## P and m computed beforehand
##  r maximizing absolute value of r-P(r) 
## (denominator in A not too small. Otherwise: completely off)

nQSD=function(n=200,r)
{
  un=sapply(0:(n-1),function (k) exp (2 *pi *1i * k/n ))
  A=matrix(nrow=n,ncol=n)
  for (i in 1:n) 
    A[i,]= sapply(1:n, function (j) r*un[j] / (r*un[j]-P_cmp(r*un[i])))
  A=A-n*diag(m,n)
  ev=eigen(A)
  l=length(ev$values)
  # You want this number REALLY small (e.g. 10^-10)
  print(ev$values[l])
  vec=fft(ev$vectors[,l])
  vec=vec*r^(-c(0:(n-1)))
  vec=-vec/vec[1]
  vec=vec[2:n]
  return(vec)
}
nQSD_cmp=cmpfun(nQSD)

# List Cyc Samples
LCycSample = function(x=1,N=100) # x serves as an anchor point, N is number of cycles from x to 0 or x
{
  path=list() # list of all paths of the process
  cycle=1 # number of current cycle
  while(cycle<=N) # N cycles 
  { # Initiate a cycle
    cur=x # current state 
    path[[cycle]]=c(cur) 
    t1=1 # time (+1) 
    while(t1==1 | !(cur %in% c(0,x)))
    {
    cur=psample_cmp(cur,x)
    t1=t1+1
    path[[cycle]][t1]=cur
    }
    cycle=cycle+1
  }
  return(path)
}
LCycSample_cmp = cmpfun(LCycSample)

# Get estimate from sampling
x=3
N=10^6
path = LCycSample_cmp(x,N) # Get our sample
cyclength=sapply(path,length) # Get the lengths of each cycle
laststate=sapply(1:N,function(n) path[[n]][cyclength[n]]) # record last state of each cycle
cyclength=cyclength-1 #Change to "real" time
path=sapply(1:N,function(n) path[[n]][1:cyclength[n]]) # remove last state in each path (last is recorded anyway)
xending = which(laststate==x) # all those ending at x
m_est = function (eps=0.001,c=0.2)
{ 
 vals = cyclength[xending]
 wrkl=0
 onest=0
 while(abs(onest-1)>eps)
  {
    oneder = sum(sapply(vals,function(n) n*exp(wrkl*n)))/N
    wrkl = wrkl + min((1-onest)/oneder,c)
    onest = sum(sapply(vals,function(n) exp(wrkl*n)))/N
    if (onest>1) print("exceeding!")
 }
 print(onest)
 return(exp(-wrkl))
}


close2one=sum(mi^cyclength[xending])/N
close2one # theoretically, this is 1
mux_est = mean(mi^cyclength-1)
mux_est = (mi-1)/mux_est   # Using formula for estimating mu (x)
mux_est # Our estimate for mu(x)



# Estimating all others 

powers = list("powers"=unlist(sapply(1:N, function (n) 0:(cyclength[n]-1))))
state_powers = aggregate (x = powers, by =list("state"=unlist(path)),FUN=c)

states=state_powers[[1]]
mu_est = sapply(1:length(state_powers[[2]]), function(n) sum(mi^state_powers[[2]][[n]]))
mu_est = mu_est / sum(mu_est)

plot(states,log(mu_est))
points(states,log((1-q)*q^(states-1)),col='red',pch="*")






