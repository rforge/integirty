Pi <-
function(th,it,D=1){
a<-it[,1]
b<-it[,2]
c<-it[,3]
d<-it[,4]
e<-exp(D*a*(th-b))
Pi<-c+(d-c)*e/(1+e)
Pi[Pi==0]<-1e-10
Pi[Pi==1]<-1-1e-10
dPi<-D*a*e*(d-c)/(1+e)^2
d2Pi<-D^2*a^2*e*(1-e)*(d-c)/(1+e)^3
d3Pi<-D^3*a^3*e*(d-c)*(e^2-4*e+1)/(1+e)^4
res<-list(Pi=Pi,dPi=dPi,d2Pi=d2Pi,d3Pi=d3Pi)
return(res)}

