data(DengueSimR02)

r.max<-seq(20,1000,20)
r.min<-seq(0,980,20)

#Lets see if there's a difference in spatial dependence by time case occurs
type<-2-(DengueSimR02[,"time"]<75)
tmp<-cbind(DengueSimR02,type=type)

typed.pi<-get.pi.typed(tmp,typeA=1,typeB=2,r=r.max,r.low=r.min)
typed.pi.type.null<-get.pi.typed.permute(tmp,typeA=1,typeB=2,r=r.max,r.low=r.min,permutations=100)
