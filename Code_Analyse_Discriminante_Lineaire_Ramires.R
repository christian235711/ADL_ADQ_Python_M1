
################################ ADL et ADQ ###################################

library(matlib)
library(ade4)
library(adegraphics)
library(rgl) 
library(tidyverse)
library(cowplot)
library(MASS)
library(klaR)


## La base de données est Iris
data_iris<-iris

## Pour chaque classe de la variable "Species", on attribue un numéro afin de faciliter le programmme

# 1:setosa ; 2:versicolor; 3:virginica

val_123<-data.frame(as.numeric(data_iris$Species))

##  Maintenant, on modifie "Species" dans la base de données
data_iris<-data.frame(iris[,-5] , val_123)

## On prend 110 individus parmi 150 de manière aléatoire
echa <- sample(1:150,110)

## Les variables
val_predicteur<-data_iris[echa,] # 110 individus

val_non_predicteur <-data_iris[-echa,] # 40 individus

## Représentation graphique des 110 individus

b1<-ggplot(iris[echa,], aes(x= Sepal.Width, y=Sepal.Length, col= Species) ) + geom_point(aes(shape=Species))
b2<-ggplot(iris[echa,], aes(x= Petal.Width, y=Petal.Length, col= Species) ) + geom_point(aes(shape=Species))
b3<-ggplot(iris[echa,], aes(x= Petal.Width, y=Sepal.Length, col= Species) ) + geom_point(aes(shape=Species))
b4<-ggplot(iris[echa,], aes(x= Sepal.Width, y=Petal.Width, col= Species) ) + geom_point(aes(shape=Species))

plot_grid(b1, b2, b3,b4, ncol =2, nrow = 2)+ggtitle("Représentation des fleurs par rapport à ses variables explicatives \n")+
  theme(plot.title = element_text(hjust = 0.5))


a1<-ggplot(iris[echa,])+geom_density(aes( x=Petal.Length, fill=Species ),alpha=0.25) 
a3<-ggplot(iris[echa,])+geom_density(aes( x=Petal.Width, fill=Species ),alpha=0.25)
a2<-ggplot(iris[echa,])+geom_density(aes( x=Sepal.Length, fill=Species ),alpha=0.25)
a4<-ggplot(iris[echa,])+geom_density(aes( x=Sepal.Width, fill=Species ),alpha=0.25)

plot_grid(a1, a2, a3,a4, ncol =2, nrow = 2)+ggtitle("Distribution des fleurs par rapport à ses variables explicatives \n")+
  theme(plot.title = element_text(hjust = 0.5))


partimat(Species~., data=iris[echa,], method='lda',main="Frontières de décision pour l'ADL" )
partimat(Species~., data=iris[echa,], method='qda',main="Frontières de décision pour l'ADQ" )


################## Il faut estimer pi^, u^,sigma^ 

n= length(echa) # nombre d'individus pris pour faire la prédiction 
K=3 # nombre de classes

####### n_k
n_k<-rep(0,K) # vecteur vide
for(j in 1:length(n_k) ){ #  1:3
  
  for( k in 1:length(echa)) { # parmi les 110 valeurs 
    if (val_predicteur[k,5] ==j)
        n_k[j]=n_k[j]+1
  }
}
n_k
######## pi_k
pi_k<-rep(0,K) # vecteur vide
for(j in 1:length(n_k) ){
    pi_k[j] = n_k[j]/n
}
pi_k

####### u_k
u_k<-matrix(data=0,nrow=K,ncol=4) # matrice vide
S_i <- matrix(data = 0, nrow = K, ncol = 4) # matrice vide
X<- matrix(data = 0, nrow = 1, ncol = 4)  # matrice vide

for(j in 1:length(n_k) ){     
    for( k in 1:length(echa)) { # parmi les 110 valeurs 
      if (val_predicteur[k,5] ==j)
        X<-matrix(as.numeric(unlist(val_predicteur[k,1:4])), nrow=nrow(val_predicteur[k,1:4]))
        S_i[j,]=S_i[j,]+X
    }
  u_k[j,]<-S_i[j,]/n_k[j]
}
S_i 
u_k

##############  ANALYSE DISCRIMINANTE LINEAIRE ####################

########### sigma

M<-matrix(data = 0, nrow = 4, ncol = 4) # matrice vide
sigma<-matrix(data = 0, nrow = 4, ncol = 4) # matrice vide

for(j in 1:length(n_k)){    
   for(k in 1:length(echa) ){
      if(val_predicteur[k,5] == j) {    
        X<-matrix(as.numeric(unlist(val_predicteur[k,1:4])), nrow=nrow(val_predicteur[k,1:4]))
        M <- M +  t(X - u_k[j,])%*%(X - u_k[j,]) 
      }
   }
}

sigma <- M/(n-K)

###### La fonction linéaire discriminante (Delta)

NR <- matrix(data=0,ncol=4, nrow = 150- length(echa) ) # petit x
delta_k <- matrix(data=0,ncol=3, nrow = 150- length(echa) )
Delta <- matrix(data=0,ncol=1, nrow = 150- length(echa) )


for(i in 1:(150- length(echa)) ) {  # 40
    
  for(j in 1:length(n_k) ) {   # 3
      NR[i,]<-matrix(as.numeric(unlist(val_non_predicteur[i,1:4])), nrow=nrow(val_non_predicteur[i,1:4]))

      delta_k[i,j]<-log(pi_k[j])+ NR[i,]%*%inv(sigma)%*%( as.matrix(u_k[j,]) ) - 0.5*t( as.matrix( u_k[j,] ) )%*%inv(sigma)%*%( as.matrix( u_k[j,] ) )
  }
  
  # On choisit le maximum de vraisemblance
  for(j in 1:length(n_k) ) {  # 3
    if ( max(delta_k[i,]) == delta_k[i,j] )
      Delta[i] <-j
  }    
}

Delta

## Nombre de prévisions correctes par classe

vec_bien_pred <- rep(0, length(n_k))
for( k in 1: length(n_k)){ # 3
  for (i in 1:(150-length(echa)) ){ # 40
  
    if( Delta[i] == val_non_predicteur[i,5] && val_non_predicteur[i,5]==k )
      vec_bien_pred[k]<- vec_bien_pred[k] + 1 
  }
}

vec_bien_pred # nombre d'individus par classe qu'on trouve à l'aide de l'échantillon d'apprentissage
c(sum(val_non_predicteur[,5]==1),sum(val_non_predicteur[,5]==2),sum(val_non_predicteur[,5]==3)) # nombre d'individus par classe qu'on devrait trouver


## Erreur globale

rr=1- sum (vec_bien_pred)/sum(c(sum(val_non_predicteur[,5]==1),sum(val_non_predicteur[,5]==2),sum(val_non_predicteur[,5]==3)))
rr


######################  ANALYSE DISCRIMINANTE QUADRATIQUE ############

###### sigma_k

M<-matrix(data = 0, nrow = 4, ncol = 4) # matrice vide

sigma_1<-matrix(data = 0, nrow = 4, ncol = 4) # matrice vide
sigma_2<-matrix(data = 0, nrow = 4, ncol = 4)
sigma_3<-matrix(data = 0, nrow = 4, ncol = 4)

for(j in 1:length(n_k)){    
  for(k in 1:length(echa) ){
    if(val_predicteur[k,5] == j) {    
      X<-matrix(as.numeric(unlist(val_predicteur[k,1:4])), nrow=nrow(val_predicteur[k,1:4]))
      M <-M+ t(X - u_k[j,])%*%(X - u_k[j,]) 
    }
  }
  
  if(j==1){ # classe 1
     sigma_1 <- M/(n_k[j]-1)
     M<-0
  }
  if(j==2){ # classe 2
    sigma_2 <- M/(n_k[j]-1)
    M<-0 
  }
  if(j==3){ # classe 3
    sigma_3 <- M/(n_k[j]-1)
    M<-0 
  }
}

###### La fonction quadratique discriminante (Delta)

NR <- matrix(data=0,ncol=4, nrow = 150- length(echa) ) # vecteur vide
delta_k <- matrix(data=0,ncol=3, nrow = 150- length(echa) )
Delta <- matrix(data=0,ncol=1, nrow = 150- length(echa) )


for(i in 1:(150- length(echa)) ) {  # 40
  
  for(j in 1:length(n_k) ) {   # 3
    NR[i,]<-matrix(as.numeric(unlist(val_non_predicteur[i,1:4])), nrow=nrow(val_non_predicteur[i,1:4]))
    
    if (j==1){ # classe 1
    delta_k[i,j]<-log(pi_k[j]) -0.5*(log(det(sigma_1)))  -0.5*( NR[i,]-u_k[j,] )%*%inv(sigma_1)%*%( as.matrix( NR[i,]-u_k[j,]) )
    }
    if (j==2){ # classe 2
      delta_k[i,j]<-log(pi_k[j]) -0.5*(log(det(sigma_2))) - 0.5*( NR[i,]-u_k[j,] )%*%inv(sigma_2)%*%( as.matrix( NR[i,]-u_k[j,]) )
    }
    if (j==3){ # classe 3
      delta_k[i,j]<-log(pi_k[j]) -0.5*(log(det(sigma_3))) -0.5*( NR[i,]-u_k[j,] )%*%inv(sigma_3)%*%( as.matrix( NR[i,]-u_k[j,]) )
    }
  }
  
  # On choisit le maximum de vraisemblance
  for(j in 1:length(n_k) ) {  # 3
    if ( max(delta_k[i,]) == delta_k[i,j] )
      Delta[i] <-j
  }    
}
Delta

## Nombre de prévisions correctes par classe

vec_bien_pred <- rep(0, length(n_k))
for( k in 1: length(n_k)){ # 3
  for (i in 1:(150-length(echa)) ){ # 40
    
    if( Delta[i] == val_non_predicteur[i,5] && val_non_predicteur[i,5]==k )
      vec_bien_pred[k]<- vec_bien_pred[k] + 1 
  }
}

vec_bien_pred # nombre d'individus par classe qu'on trouve à l'aide de l'échantillon d'apprentissage
c(sum(val_non_predicteur[,5]==1),sum(val_non_predicteur[,5]==2),sum(val_non_predicteur[,5]==3)) # nombre d'individus par classe qu'on devrait trouver

### Erreur globale

rr=1- sum (vec_bien_pred)/sum(c(sum(val_non_predicteur[,5]==1),sum(val_non_predicteur[,5]==2),sum(val_non_predicteur[,5]==3)))
rr

