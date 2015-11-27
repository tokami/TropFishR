#Curves to be added include:
#tt.richards, for richards fit to trouser trawl data
#gamma, for net selectivity.

selncurves=function(rtype) {
  switch(rtype,
   "norm.loc"={
     r=function(lens,Meshsize,th) {
       relsize=Meshsize/Meshsize[1]
       seln=exp(-(lens-th[1]*relsize)^2/(2*th[2]^2)) 
       return(seln) } },
   "norm.sca"={
     r=function(lens,Meshsize,th) {
       relsize=Meshsize/Meshsize[1]
       seln=exp(-(lens-th[1]*relsize)^2/(2*th[2]^2*relsize^2)) 
       return(seln) } },
   "lognorm"={
     r=function(lens,Meshsize,th) {
       relsize=Meshsize/Meshsize[1]
       seln=(relsize/lens)*exp(th[1]-th[2]^2/2)
       seln=seln*exp( -(log(lens)-th[1]-log(relsize))^2/(2*th[2]^2) )
       return(seln) } },
   "binorm.sca"={
     r=function(lens,Meshsize,th) {
       relsize=Meshsize/Meshsize[1]
       seln1=exp(-(lens-th[1]*relsize)^2/(2*th[2]^2*relsize^2))
       seln2=exp(-(lens-th[3]*relsize)^2/(2*th[4]^2*relsize^2))
       p=exp(th[5])/(1+exp(th[5])) #i.e., th[5]=logit(p)
       seln=p*seln1+(1-p)*seln2
       return(seln) } }, 
   "bilognorm"={
     r=function(lens,Meshsize,th) {
       relsize=Meshsize/Meshsize[1]
       seln1=(relsize/lens)*exp(th[1]-th[2]^2/2)
       seln1=seln1*exp( -(log(lens)-th[1]-log(relsize))^2/(2*th[2]^2) )
       seln2=(relsize/lens)*exp(th[3]-th[4]^2/2)
       seln2=seln2*exp( -(log(lens)-th[3]-log(relsize))^2/(2*th[4]^2) )
       p=exp(th[5])/(1+exp(th[5])) #i.e., th[5]=logit(p)
       seln=p*seln1+(1-p)*seln2
       return(seln) } },
   "tt.logistic"={
     r=function(lens,Meshsize,th) {
       control=(Meshsize==Meshsize[1])
       p=exp(th[3])/(1+exp(th[3])) #i.e., th[3]=logit(p)
       wk=exp(th[1]+th[2]*lens)
       lselect=wk/(1+wk)
       seln=(1-p)*control+p*lselect*(1-control)
       return(seln) } },    
   stop(paste("\n",rtype, "not recognised, possible curve types are \n", 
        "\"norm.loc\", \"norm.sca\", \"lognorm\" \n", 
        "\"binorm.sca\", \"bilognorm\", and \"tt.logistic\"")) 
  )#End of switch
  return(r) }

