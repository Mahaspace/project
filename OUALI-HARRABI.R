setwd("C:/Users/maha/Documents/ensta/cours/stat/project/project")
rm(list=objects())
graphics.off()
##*******************************************************************
##**************************partie 2**********************************
#question 1
set.seed(1)
# simulation d'une loi de pareto (tau)
#génération d'un échantillon de taille 50 d'une loi uniforme
x=runif(50,0,1)
tau0=3
y=x^(-1/tau0)
hist(y, proba=TRUE)
curve(tau0/x^(tau0+1),add=TRUE)
#question2
B=100
n=50
alpha=0.05
Tchap=replicate(B,1/mean(log(1/runif(n)^(1/tau0))))

hist(Tchap,proba=TRUE)
##*********************************************************************
##estiamtion  du niveau du test 
##***************pour le test gamma
test=n*tau0/Tchap
est=mean(test > qgamma(0.95,n,1))
est
# pour le test normal 
test2=sqrt(n)*(Tchap/tau0-1)
est2=mean(test2 > qnorm(0.95,0,1))
est2
##***************************************************************
##Courbe de puissance 

curve(pgamma((x*qgamma(0.05,n,1))/tau0,n,1), xlim=c(2,5))
abline(h=0.05)
# estimateur biaisé
##******************************************************************

#estimation de la puissance 
tau=seq(2,5,0.01)
# estimation de la puissance
tau=seq(2,5,0.01)
y1=seq(2,5,0.01)
y2=seq(2,5,0.01)
for (i in 1:length(tau))
{
  y1[i]=mean((replicate(B,1/mean(log(1/runif(n)^(1/tau[i]))))>n*tau0/qgamma(0.05,n,1)));
  y2[i]=mean((replicate(B,1/mean(log(1/runif(n)^(1/tau[i]))))>tau0*qnorm(0.95,0,1)/sqrt(n)+tau0))
}
points(tau,y1,col=blues9)
points(tau,y2,col="red")

# pour n=50  le nuage de points de l'estimation de la puissance entoure la courbe théorique , les estimations et la valeur théorique sont très proche.
# lorsque n diminue n=5, le nuage de point s'éligne de la courbe théorique
# le premier estimateur qui suit la loi gamma reste plus  precis 
# c'est exactement la loi des grand nombre , il faut augmenter l'echantillon pour s'approcher de la valeur réelle 
##**************************************************
#pour n=5
# courbe de puissance 

# curve(pgamma((x*qgamma(0.05,n,1))/tau0,n,1), xlim=c(2,5))
# abline(h=0.05)
# 
# #estimation de la puissance 
# tau=seq(2,5,0.01)
# # estimation de la puissance
# tau=seq(2,5,0.01)
# y1=seq(2,5,0.01)
# y2=seq(2,5,0.01)
# for (i in 1:length(tau))
# {
#   y1[i]=mean((replicate(B,1/mean(log(1/runif(n)^(1/tau[i]))))>n*tau0/qgamma(0.05,n,1)));
#   y2[i]=mean((replicate(B,1/mean(log(1/runif(n)^(1/tau[i]))))>tau0*qnorm(0.95,0,1)/sqrt(n)+tau0))
# }
# points(tau,y1,col=blues9)
# points(tau,y2,col="red")
##******************************************
#question 4
  f=function(q,n,alpha=0.05)
  {
    g=q^n*exp(-q)-(qgamma(pgamma(q,n,1)+1-alpha,n,1)^n*exp(-qgamma(pgamma(q,n,1)+1-alpha,n,1)))
    return(g)
  }
 
curve(f(x,n,alpha=0.05),from=qgamma(0.01,n,1), to=qgamma(0.05,n,1))
abline(0,0)
sol=uniroot(f,c(37,38),n) # d'après la courbe, la solution se situe entre 37 et 38
q1=sol$root
q1

q2=qgamma(pgamma(q1,n,1)+1-alpha,n,1)
q2
##************************************
# pour n=5
# 
# f=function(q,n,alpha=0.05)
# {
#   g=q^n*exp(-q)-(qgamma(pgamma(q,n,1)+1-alpha,n,1)^n*exp(-qgamma(pgamma(q,n,1)+1-alpha,n,1)))
#   return(g)
# }
# 
# curve(f(x,n,alpha=0.05),from=qgamma(0.01,n,1), to=qgamma(0.05,n,1))
# abline(0,0)
# sol=uniroot(f,c(qgamma(0.01,n,1),qgamma(0.049,n,1)),n)
# q1=sol$root
# q1
 # n=5

#test bilatéral

curve(1-pgamma(qgamma(1-0.05/2,n,1)*x/tau0,n,1)+pgamma(qgamma(0.05/2,n,1)*x/tau0,n,1),from=2, to=5)
# test non  biasé
curve(1-pgamma(q2*x/tau0,n,1)+pgamma(q1*x/tau0,1),from=2, to=5)


#superposition des deux courbes 


#####
###partie 3


##***************************************************
# les plus grandes villes françaises
data=read.delim2("VillesFrançaises.txt",comment.char = "#",fileEncoding = "UTF8")
V=data$p2018
sort(V,decreasing = TRUE)

# les observations renormalisées 

VF=V[1:39]/V[40]
nVF=length(VF)
hist(VF,proba=TRUE)

# estimationde tau
tau_F=1/mean(log(VF))
curve(tau_F/x^(tau_F+1), add=TRUE)
q_emp <- quantile(VF, seq(1,nVF,1)/(nVF+1), names=FALSE)
qpareto<-function(x, tau_F){1/(1-x)^(1/tau_F)*(x>=0)}
quantiles_pareto <- qpareto(seq(1,nVF,1)/(nVF+1), tau_F)
plot(quantiles_pareto, q_emp, main="Courbe quantiles empiriques = f(quantiles th?oriques de la loi de Pareto) \n pour les villes en France", type='l')
abline(0, 1)
#appliquer le test de Kolmogorov-Smirnov ? l'?chantillon
ppareto<-function(x, tau){1-1/x^tau*(x>=1)}
ks.test(VF,ppareto,tau=tau_F)
##**************************************************
# les plus grandes villes allemandes 
data2=read.delim2("VillesAllemandes.txt",comment.char = "#",fileEncoding = "UTF8")
V2=data2$p2020
sort(V2,decreasing = TRUE)

  VA=V2[1:39]/V2[40]
nVA=length(VA)
hist(VA,proba=TRUE)


# estimationde tau
tau_A=1/mean(log(VA))
curve(tau_A/x^(tau_A+1), add=TRUE)
q_emp <- quantile(VA, seq(1,nVA,1)/(nVA+1), names=FALSE)
qpareto<-function(x, tau){1/(1-x)^(1/tau)*(x>=0)}
quantiles_pareto <- qpareto(seq(1,nVA,1)/(nVA+1), tau_A)
plot(quantiles_pareto, q_emp, main="Courbe quantiles empiriques en fonction des quantiles theoriques de la loi de Pareto) \n pour les villes en Allemagne", type='l')
abline(0, 1)
#appliquer le test de Kolmogorov-Smirnov ? l'?chantillon
ppareto<-function(x, tau){1-1/x^tau*(x>=1)}
ks.test(VA,ppareto,tau=tau_A)

##******************************
# pour les villes du Royaume Uni 
  data3=read.delim2("VillesUK.txt",comment.char = "#",fileEncoding = "UTF8")
V3=data3$p2019
sort(V3,decreasing = TRUE)

VUK=V3[1:39]/V3[40]
nVUK=length(VUK)
hist(VUK,proba=TRUE)


# estimationde tau
tau_UK=1/mean(log(VUK))
curve(tau_UK/x^(tau_UK+1), add=TRUE)
q_emp <- quantile(VUK, seq(1,nVUK,1)/(nVUK+1), names=FALSE)
qpareto<-function(x, tau){1/(1-x)^(1/tau)*(x>=0)}
quantiles_pareto <- qpareto(seq(1,nVUK,1)/(nVUK+1), tau_UK)
plot(quantiles_pareto, q_emp, main="Courbe quantiles empiriques en fonction des quantiles theoriques de la loi de Pareto) \n pour les villes de l'UK", type='l')
abline(0, 1)
#appliquer le test de Kolmogorov-Smirnov ? l'?chantillon
ppareto<-function(x, tau){1-1/x^tau*(x>=1)}
ks.test(VA,ppareto,tau=tau_UK)
