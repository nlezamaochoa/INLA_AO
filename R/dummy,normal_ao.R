a
#################INLA TO PREDICT THE PRESENCE OF MOBULA MOBULAR IN THE EASTERN PACIFIC OCEAN######
####PRESENCES AND 10000 RANDOM ABSENCES FROM THE REST OF THE BYCATCH DATABASE#########
#COMPARE RESULTS WITH GAMS####
graphics.off()
par("mar")
par(mar=c(1,1,1,1))
#CARGAR LIBRERIAS
rm(list=ls())
data(wrld_simpl)
library(dismo)
library(maptools) 
library(ncdf4)
library(raster)
library(rasterVis)
library(RColorBrewer)
library(zoo)
library(INLA)
require(vegan)
library(raster)
library(lattice)
library(maptools)
require(gridExtra)
require(rworldmap)
require(rworldxtra)
library(maps)
library(mapdata)
library(mapproj)
library (rgrs)
library(rgeos)
library(sp)
library(rgdal)
library(gdata)
library(RColorBrewer)
library(rgr)
library(fields)
library(geoR)
library(gstat)
library(splancs)
library(RandomFields)
library(raster)
library(gstat)
library(dummies)

Year
month
type
depth
chl
nit
sal
ssh

###############Read data###########################
dados <- read.csv("C://Users//Nerea//Desktop//NEREA.csv")
names(dados)
str(dados)

dados$school_type
levels(dados$school_type)
dados$school_type3=dummy(dados$school_type)
#dados$Type <- as.numeric(dados$Type)
dados$school_type3


#With this function you split the dataset in tow
#80% =training to run the model again
#20% = test to validate the model
#do this k times
group <- kfold (dados,5)

##presences (occ_P)
training<-dados [group !=1,]
test<- dados [group ==1,]
dim(training)
dim(test)

######Build a coords vector##########
names(dados)
coords <- cbind(training[,8:9])
colnames(coords) <- c("longitude","latitude")
head(coords)
str(coords)
plot(coords)

##############extension#######
ext <- extent(-50, 15,-25, 25)
par(mfrow=c(1,1))
plot(ext)
############Environmental data###################

#A partir de tus datos ambientales observados en cada localidad vamos a
#hacer un ?nico predictor con la media. Seria eso:
#creamos los raster para hacer luego la prediction
r<-raster(xmn=-50,xmx=15,ymn=-25,ymx=25,nrows=61,ncols=64)#extent de tu area
xy=cbind(training$lon,training$lat)#la y lon de todos tus datos
xy
xy<-as.data.frame(xy)
head(xy)


#Chl
chl=rasterize(xy, r, dados$chl, fun=mean)#la variable que quieremos
chl

O2=rasterize(xy, r, dados$O2, fun=mean)#la variable que quieremos
O2


#Ni
Ni=rasterize(xy, r, dados$Ni,fun=mean)#la variable que quieremos
Ni
#SSH
SSH=rasterize(xy, r, dados$SSH,fun=mean)#la variable que quieremos
SSH

plot(SSH, col=tim.colors(100)[1:100])
#month
Month=rasterize(xy, r, dados$Month,fun=mean)#la variable que quieremos
Month
#type
school_type=rasterize(xy, r, dados$school_type3, fun=mean)#la variable que quieremos
school_type
plot(school_type, col=tim.colors(100)[1:100])

#sal
Sal=rasterize(xy, r, dados$Sal, fun=mean)#la variable que quieremos
Sal

#depth
depth=rasterize(xy, r, dados$depth, fun=mean)#la variable que quieremos
depth

#year
Year=rasterize(xy, r, dados$Year, fun=mean)#la variable que quieremos
Year

#Crop
chl<-crop(chl, ext)
O2<-crop(O2, ext)
Ni<-crop(Ni,ext)
SSH<-crop(SSH,ext)
Month<-crop(Month, ext)
Sal<-crop(Sal, ext)
depth<-crop(depth, ext)
Year<-crop(Year, ext)
school_type<-crop(school_type, ext)


###########renombrar variables#####################
chl_val_val <- training$chl
chl_val_val

O2_val_val <- training$O2
O2_val_val

Ni_val_val <- training$Ni
Ni_val_val

SSH_val_val <- training$SSH
SSH_val_val

Month_val_val <- training$Month
Month_val_val

school_type_val_val <- training$school_type
school_type_val_val

Year_val_val <- training$Year
Year_val_val

Sal_val_val <- training$Sal
Sal_val_val

depth_val_val <- training$depth
depth_val_val


############################################
##A PARTIR DE AQU? CAMBIAR POR EL SCRIPT DE IOSU
##########Build mesh#######################
######################################
#NUESTRO MESH
#####################################
br.sp <- getMap(resolution = "high") 
#a partir de nuestras posiciones, hacer el mesh
coords <- as.matrix(training[,c("longitude","latitude")])
boundary=inla.nonconvex.hull(coords,convex=-0.01)#circulo azul
boundary
mesh<-inla.mesh.2d(boundary=boundary,max.edge=c(2,2),offset=c(2,2))
plot(mesh)
#max.edge: el tama?o de los traiangulso fuera de lo azul y dentro
#offset: el espacio dentro del circulo azul (primer valor, ) y denentre el circulo negro y el circulo azul. Dejar mucho espacio o no

##########Plot mesh##############
#x11()
##tiff("C:/Users/Nerea/Desktop/mesh.tif", res = 600, units = "cm", width = 20, height = 20)



##dev.off()

windows(20,20)
plot(mesh)
plot(br.sp,add=T,xlim=c(-45, 15),ylim=c(-25, 25),col="gray75")
points(coords,pch=20,col="red")
box()


###########Create a SPDE model object for a Matern-like spatial covariance function#####################
#Build the SPDE model on the mesh. SPDE define the model on the mesh nodes. To model the spatial autocorrelation
spde <- inla.spde2.matern(mesh, alpha=1) #smoothness parameter of the matern correlation matrix, alpha=1 as default
spde
#Define the nxm projector matrix to project the process at the mesh nodes to locations
A.est<- inla.spde.make.A(mesh=mesh, loc=coords, ncol=2)
A.est

###################Estimation matriz  ############
names(dados)
presence<-training$PA
presence
plot(density(presence))

#stack the model to be able to have all the covariates in the same format, ya que n? de nodos es diferente a n? de locations
est<-inla.stack(data=list(y=presence),
                A=list(A.est, 1),
                effects=list(spat=1:spde$n.spde,
                             list(b0=1,O2=O2_val_val, school_type=school_type_val_val,
                                  chl=chl_val_val, Ni=Ni_val_val, 
                                  SSH=SSH_val_val, 
                                  Month=Month_val_val)))
est

#######################Prediction matrix#############

A.pred = inla.spde.make.A(mesh)
A.pred
str(mesh)
meshcoords <- cbind(mesh$loc[,1], mesh$loc[,2])
meshcoords


chl_m<- extract(chl, meshcoords)
summary(chl_m)
#chl_m <- decostand(chl_m, "standardize", na.rm=T) #padronização

O2_m<- extract(O2, meshcoords)
summary(O2_m)
#O2_m <- decostand(O2_m, "standardize", na.rm=T) #padronização

Ni_m<- extract(Ni, meshcoords)
summary(Ni_m)
#Ni_m <- decostand(Ni_m, "standardize", na.rm=T) #padronização

SSH_m<- extract(SSH, meshcoords)
summary(SSH_m)
#SSH_m <- decostand(SSH_m, "standardize", na.rm=T) #padronização

Month_m<- extract(Month, meshcoords)
summary(Month_m)
#Month_m <- decostand(Month_m, "standardize", na.rm=T) #padronização

Type_m<- extract(Type, meshcoords)
summary(Type_m)
#Type_m <- decostand(Type_m, "standardize", na.rm=T) #padronização

#importante! no incluimos el type (porque es dolphin, school and fad)
pred<-inla.stack(data=list(y=rep(NA,spde$n.spde)),
                 A=list(A.pred, 1),
                 effects=list(spat=1:spde$n.spde,
                              list(O2=O2_m,
                                    chl=chl_m,  Ni=Ni_m,
                                    SSH=SSH_m)), tag='pred')
pred
stk=inla.stack(est, pred)#probar sin pred si tarda mucho, pero el pred va dentro (est, pred)
stk

#############Run the model######################
#fit the posterior marginal distributions for all model parameters
formula4<-y ~ 0 + b0+ f(inla.group(chl), model="rw2")                        
formula6<-y ~ 0 + b0+f(inla.group(Ni), model="rw2") 
formula7<-y ~ 0 + b0+f(inla.group(O2), model="rw2")  
formula9<-y ~ 0 + b0+  f(inla.group(SSH), model="rw2")
formula14<-y ~ 0 + b0+Type
formula15<-y ~ 0 + b0+f(inla.group(Month), model="rw2")


res4<- inla(formula4,
            family="binomial", 
            data=inla.stack.data(stk), 
            keep=FALSE,
            #control.family=list(hyper=list(prec=list(param=c(1,0.5)))),
            control.predictor=list(A=inla.stack.A(stk), compute=TRUE),
            control.inla=list(tolerance=1e-5,numint.maxfeval= 10e6),
            #control.fixed = list(expand.factor.strategy='inla'),
            control.compute = list(return.marginals=TRUE, dic=TRUE, cpo=TRUE)
)
date()
summary(res4)


res6<- inla(formula6,
            family="binomial", 
            data=inla.stack.data(stk), 
            keep=FALSE,
            #control.family=list(hyper=list(prec=list(param=c(1,0.5)))),
            control.predictor=list(A=inla.stack.A(stk), compute=TRUE),
            control.inla=list(tolerance=1e-5,numint.maxfeval= 10e6),
            #control.fixed = list(expand.factor.strategy='inla'),
            control.compute = list(return.marginals=TRUE, dic=TRUE, cpo=TRUE)
)
date()
summary(res6)


res7<- inla(formula7,
            family="binomial", 
            data=inla.stack.data(stk), 
            keep=FALSE,
            #control.family=list(hyper=list(prec=list(param=c(1,0.5)))),
            control.predictor=list(A=inla.stack.A(stk), compute=TRUE),
            control.inla=list(tolerance=1e-5,numint.maxfeval= 10e6),
            #control.fixed = list(expand.factor.strategy='inla'),
            control.compute = list(return.marginals=TRUE, dic=TRUE, cpo=TRUE)
)
date()
summary(res7)#no significant


res9<- inla(formula9,
            family="binomial", 
            data=inla.stack.data(stk), 
            keep=FALSE,
            #control.family=list(hyper=list(prec=list(param=c(1,0.5)))),
            control.predictor=list(A=inla.stack.A(stk), compute=TRUE),
            control.inla=list(tolerance=1e-5,numint.maxfeval= 10e6),
            #control.fixed = list(expand.factor.strategy='inla'),
            control.compute = list(return.marginals=TRUE, dic=TRUE, cpo=TRUE)
)
date()
summary(res9)

res14<- inla(formula14,
             family="binomial", 
             data=inla.stack.data(stk), 
             keep=FALSE,
             #control.family=list(hyper=list(prec=list(param=c(1,0.5)))),
             control.predictor=list(A=inla.stack.A(stk), compute=TRUE),
             control.inla=list(tolerance=1e-5,numint.maxfeval= 10e6),
             #control.fixed = list(expand.factor.strategy='inla'),
             control.compute = list(return.marginals=TRUE, dic=TRUE, cpo=TRUE)
)
date()
summary(res14)

res15<- inla(formula15,
             family="binomial", 
             data=inla.stack.data(stk), 
             keep=FALSE,
             #control.family=list(hyper=list(prec=list(param=c(1,0.5)))),
             control.predictor=list(A=inla.stack.A(stk), compute=TRUE),
             control.inla=list(tolerance=1e-5,numint.maxfeval= 10e6),
             #control.fixed = list(expand.factor.strategy='inla'),
             control.compute = list(return.marginals=TRUE, dic=TRUE, cpo=TRUE)
)
date()
summary(res15)


#all significant


#with all like in gams
formula18b<-y ~  0 + b0  + f(inla.group(O2), model="rw2")  +f(inla.group(chl), model="rw2") +
  f(inla.group(Ni), model="rw2")   +
  f(inla.group(SSH), model="rw2")+   f(inla.group(Month), model="rw2")+f(spat,model=spde)

res18b<- inla(formula18b,
              family="binomial", 
              data=inla.stack.data(stk), 
              keep=FALSE,
              #control.family=list(hyper=list(prec=list(param=c(1,0.5)))),
              control.predictor=list(A=inla.stack.A(stk), compute=TRUE, link=1),#link=1
              control.inla=list(tolerance=1e-5,numint.maxfeval= 10e6),
              #control.fixed = list(expand.factor.strategy='inla'),
              control.compute = list(return.marginals=TRUE,dic=TRUE, cpo=TRUE)
)
date()
summary(res18)
par(mar=c(1,1,1,1))
plot(res18)
table(dados$Type)
#DIC: 3616.90
#CPO: -2138.12

#with all like in gams
formula18c<-y ~  0 + b0  +O2  +chl +school_type + Month + Ni+f(spat,model=spde)

res18c<- inla(formula18c,
              family="binomial", 
              data=inla.stack.data(stk), 
              keep=FALSE,
              #control.family=list(hyper=list(prec=list(param=c(1,0.5)))),
              control.predictor=list(A=inla.stack.A(stk), compute=TRUE, link=1),#link=1
              control.inla=list(tolerance=1e-5,numint.maxfeval= 10e6),
              #control.fixed = list(expand.factor.strategy='inla'),
              control.compute = list(return.marginals=TRUE,dic=TRUE, cpo=TRUE)
)
date()
summary(res18c)


formula20<-y ~  0 + b0  +   school_type +   f(inla.group(Month), model="rw2")
res20<- inla(formula20,
              family="binomial", 
              data=inla.stack.data(stk), 
              keep=FALSE,
              #control.family=list(hyper=list(prec=list(param=c(1,0.5)))),
              control.predictor=list(A=inla.stack.A(stk), compute=TRUE, link=1),#link=1
              control.inla=list(tolerance=1e-5,numint.maxfeval= 10e6),
              #control.fixed = list(expand.factor.strategy='inla'),
              control.compute = list(return.marginals=TRUE,dic=TRUE, cpo=TRUE)
)
date()
summary(res20)


formula21<-y ~  0 + b0  +   school_type +   f(inla.group(Month), model="rw2")+ 
  f(inla.group(chl), model="rw2")
res21<- inla(formula21,
             family="binomial", 
             data=inla.stack.data(stk), 
             keep=FALSE,
             #control.family=list(hyper=list(prec=list(param=c(1,0.5)))),
             control.predictor=list(A=inla.stack.A(stk), compute=TRUE, link=1),#link=1
             control.inla=list(tolerance=1e-5,numint.maxfeval= 10e6),
             #control.fixed = list(expand.factor.strategy='inla'),
             control.compute = list(return.marginals=TRUE,dic=TRUE, cpo=TRUE)
)
date()
summary(res21)


formula22<-y ~  0 + b0  +   school_type +   f(inla.group(Month), model="rw2")+ 
  f(inla.group(chl), model="rw2") + f(inla.group(Ni), model="rw2")
res22<- inla(formula22,
             family="binomial", 
             data=inla.stack.data(stk), 
             keep=FALSE,
             #control.family=list(hyper=list(prec=list(param=c(1,0.5)))),
             control.predictor=list(A=inla.stack.A(stk), compute=TRUE, link=1),#link=1
             control.inla=list(tolerance=1e-5,numint.maxfeval= 10e6),
             #control.fixed = list(expand.factor.strategy='inla'),
             control.compute = list(return.marginals=TRUE,dic=TRUE, cpo=TRUE)
)
date()
summary(res22)


formula23<-y ~  0 + b0  +   school_type +   f(inla.group(Month), model="rw2")+ 
  f(inla.group(chl), model="rw2") + f(inla.group(Ni), model="rw2")+f(inla.group(O2), model="rw2")
res23<- inla(formula23,
             family="binomial", 
             data=inla.stack.data(stk), 
             keep=FALSE,
             #control.family=list(hyper=list(prec=list(param=c(1,0.5)))),
             control.predictor=list(A=inla.stack.A(stk), compute=TRUE, link=1),#link=1
             control.inla=list(tolerance=1e-5,numint.maxfeval= 10e6),
             #control.fixed = list(expand.factor.strategy='inla'),
             control.compute = list(return.marginals=TRUE,dic=TRUE, cpo=TRUE)
)
date()
summary(res23)

formula24<-y ~  0 + b0  +   school_type +   f(inla.group(Month), model="rw2")+ 
  f(inla.group(chl), model="rw2") + f(inla.group(Ni), model="rw2")+f(inla.group(O2), model="rw2")+
  f(inla.group(SSH), model="rw2")
res24<- inla(formula24,
             family="binomial", 
             data=inla.stack.data(stk), 
             keep=FALSE,
             #control.family=list(hyper=list(prec=list(param=c(1,0.5)))),
             control.predictor=list(A=inla.stack.A(stk), compute=TRUE, link=1),#link=1
             control.inla=list(tolerance=1e-5,numint.maxfeval= 10e6),
             #control.fixed = list(expand.factor.strategy='inla'),
             control.compute = list(return.marginals=TRUE,dic=TRUE, cpo=TRUE)
)
date()
summary(res24)


##SACAMOS LA PREDICCION de type of set NO ME SALE!
igr <- inla.stack.index(stk, "pred")$data
ee=exp(res18c$summary.fitted.values[igr,1])
ee

levels (school_type)

p=ee/(1+ee)
p

Type_m=as.factor(school_type_m)
df <- data.frame(school_type_m, p)
df <- na.omit(df) 
library(ggplot2)
windows(20,20)
ggplot(df)+geom_boxplot(aes(school_type_m,p))


#prediction resto de variables
results<-as.data.frame(meshcoords)
results$p<-p
dados$SSH
##SACAMOS LA VARIBALE: chl, ssh, Ni, Type, Month, O2
cc<-as.data.frame(meshcoords)
head(cc)
CHL_v<-extract(SSH,cc)#CAMBIAR POR CADA VARIABLE
head(CHL_v)
CHL_v<-extract(chl,cc)#CAMBIAR POR CADA VARIABLE
head(CHL_v)
CHL_v<-extract(Ni,cc)#CAMBIAR POR CADA VARIABLE
head(CHL_v)
CHL_v<-extract(O2,cc)#CAMBIAR POR CADA VARIABLE
head(CHL_v)
#chl, type, month, O2, Ni, ssh, 

##GRAFICA
library(ggplot2)
#jpeg("CHL_tuna.jpeg",  width = 2200, height = 1600, res = 300)
windows(20,20)
qplot(CHL_v, results$p, geom='smooth', span =0.4)+ 
  theme(axis.text.x = element_text(size=14),axis.title=element_text(size=14),
        axis.text.y = element_text(size=14))+ xlab("SSH")+ylab("Probability of presence")

windows(20,20)
qplot(CHL_v, results$p, geom='smooth', span =0.4)+ 
  theme(axis.text.x = element_text(size=14),axis.title=element_text(size=14),
        axis.text.y = element_text(size=14))+ xlab("Chl")+ylab("Probability of presence")

windows(20,20)
qplot(CHL_v, results$p, geom='smooth', span =0.4)+ 
  theme(axis.text.x = element_text(size=14),axis.title=element_text(size=14),
        axis.text.y = element_text(size=14))+ xlab("Ni")+ylab("Probability of presence")

windows(20,20)
qplot(CHL_v, results$p, geom='smooth', span =0.4)+ 
  theme(axis.text.x = element_text(size=14),axis.title=element_text(size=14),
        axis.text.y = element_text(size=14))+ xlab("O2")+ylab("Probability of presence")





tiff("C:/Users/Nerea/Desktop/chl.tif", res = 600, units = "cm", width = 20, height = 20)
tiff("C:/Users/Nerea/Desktop/Ni.tif", res = 600, units = "cm", width = 20, height = 20)
tiff("C:/Users/Nerea/Desktop/O2.tif", res = 600, units = "cm", width = 20, height = 20)
tiff("C:/Users/Nerea/Desktop/SSH.tif", res = 600, units = "cm", width = 20, height = 20)
tiff("C:/Users/Nerea/Desktop/Month.tif", res = 600, units = "cm", width = 20, height = 20)



dev.off()

#spatial effect
pgrid<-inla.mesh.projector(mesh,dims = c(600,600))
prd.m<-inla.mesh.project(pgrid,res18c$summary.ran$spat$mean)#mean 
prd.m
prd.s<-inla.mesh.project(pgrid,res18b$summary.ran$spat$sd)#standard deviation

matrix<-cbind(res18$summary.ran$spat$mean,meshcoords)#ir cambiando con el resto de valores
head(matrix)
matrix<-as.data.frame(matrix)
colnames(matrix)<-c("spatial","Lon","Lat")
head(matrix)

#save
write.csv(matrix,"C:/Users//modelpcuser2//Desktop/spat_nerea.csv")

table(xy.in<-inout(pgrid$lattice$loc,coords))
prd.m[!xy.in]<-prd.s[!xy.in]<-NA

#spatial effect mean
rast0p <- raster(prd.m)
rotate <- function(rast, angle=90, resolution=res(rast)) {
  y <- rast; crs(y) <- "+proj=aeqd +ellps=sphere +lat_0=90 +lon_0=0"
  projectRaster(y, res=resolution, 
                crs=paste0("+proj=aeqd +ellps=sphere +lat_0=90 +lon_0=", -angle))
}


rast0p<-rotate(rast0p)

crs(rast0p)<- CRS("+proj=longlat +datum=WGS84")
ext<-extent(-50, 15,-20, 32)#change extent
rast0p<-setExtent(rast0p, ext)#now you have a raster


#MAP
windows(20,20)
plot(rast0p, col=tim.colors(100)[1:100], main="M. mobular - Mean Spatial effect", axes=T)#you will need library(rworldxtra)
coast<- getMap(resolution = "high") #you will need library(rworldmap)
plot(coast, col='dark grey',add=T)
box()


#response value
pgrid<-inla.mesh.projector(mesh,dims = c(600,600))
igr <- inla.stack.index(stk, "pred")$data
head(igr)
ms.m2<-inla.mesh.project(pgrid,(res18c$summary.fitted.values[igr,1]))#mean of random field
ms.m3<-inla.mesh.project(pgrid,(res18b$summary.fitted.values[igr,2]))#sd
ms.m4<-inla.mesh.project(pgrid,(res18b$summary.fitted.values[igr,3]))#25%
ms.m5<-inla.mesh.project(pgrid,(res18b$summary.fitted.values[igr,4]))#97%

matrix<-cbind(res18$summary.fitted.values[igr,1],meshcoords)#ir cambiando con el resto de valores
head(matrix)
matrix<-as.data.frame(matrix)
colnames(matrix)<-c("spatial","Lon","Lat")
head(matrix)
write.csv(matrix,"C:/Users//modelpcuser2//Desktop/lineal.csv")
#SCRIPT MARIA PARA PLOTEAR LA PREDICTION CON MAP
rast02 <- raster(ms.m2)
rast03 <- raster(ms.m3)
rast04 <- raster(ms.m4)
rast05 <- raster(ms.m5)


rotate <- function(rast, angle=90, resolution=res(rast)) {
  y <- rast; crs(y) <- "+proj=aeqd +ellps=sphere +lat_0=90 +lon_0=0"
  projectRaster(y, res=resolution, 
                crs=paste0("+proj=aeqd +ellps=sphere +lat_0=90 +lon_0=", -angle))
}


rast0p<-rotate(rast02)#rotar 4 veces

crs(rast0p)<- CRS("+proj=longlat +datum=WGS84")
ext<-extent(-50, 15,-20, 32)#change extent
rast0<-setExtent(rast0p, ext)#now you have a raster


#MAP
windows(20,20)
plot(rast0, col=tim.colors(100)[1:100] ,main="Mean prediction M. mobular", axes=T)#you will need library(rworldxtra)
coast<- getMap(resolution = "high") #you will need library(rworldmap)
plot(coast, col='dark grey',add=T)
box()





rotate <- function(rast, angle=90, resolution=res(rast)) {
  y <- rast; crs(y) <- "+proj=aeqd +ellps=sphere +lat_0=90 +lon_0=0"
  projectRaster(y, res=resolution, 
                crs=paste0("+proj=aeqd +ellps=sphere +lat_0=90 +lon_0=", -angle))
}

rast02<-rotate(rast02)
rast03<-rotate(rast03)
rast04<-rotate(rast04)
rast05<-rotate(rast05)



rast0p<-rotate(rast02)#

crs(rast02)<- CRS("+proj=longlat +datum=WGS84")
ext<-extent(-45,30,-25,25)#change extent
rast02<-setExtent(rast02, ext)#now you have a raster

crs(rast03)<- CRS("+proj=longlat +datum=WGS84")
ext<-extent(-150, -70,-20, 32)#change extent
rast03<-setExtent(rast03, ext)#now you have a raster

crs(rast04)<- CRS("+proj=longlat +datum=WGS84")
ext<-extent(-150, -70,-20, 32)#change extent
rast04<-setExtent(rast04, ext)#now you have a raster

crs(rast05)<- CRS("+proj=longlat +datum=WGS84")
ext<-extent(-150, -70,-20, 32)#change extent
rast05<-setExtent(rast05, ext)#now you have a raster

#MAP
windows(20,20)
plot(predict_present, col=tim.colors(100)[1:100] ,main="Mean prediction M. mobular", axes=T)#you will need library(rworldxtra)
coast<- getMap(resolution = "high") #you will need library(rworldmap)
plot(coast, col='dark grey',add=T)
box()

windows(20,20)
plot(rast03, col=tim.colors(100)[1:100], main="SD M. mobular", axes=T)#you will need library(rworldxtra)
coast<- getMap(resolution = "high") #you will need library(rworldmap)
plot(coast, col='dark grey',add=T)
box()

windows(20,20)
plot(rast04, col=tim.colors(100)[1:100], main="25% prediction M. mobular", axes=T)#you will need library(rworldxtra)
coast<- getMap(resolution = "high") #you will need library(rworldmap)
plot(coast, col='dark grey',add=T)
box()

windows(20,20)
plot(rast05, col=tim.colors(100)[1:100], zlim=c(0,1),  main="97% prediction M. mobular", axes=T)#you will need library(rworldxtra)
coast<- getMap(resolution = "high") #you will need library(rworldmap)
plot(coast, col='dark grey',add=T)
box()

############# EXTRAEMOS PREDICCION ##############
igr <- inla.stack.index(stk, "pred")$data

ee=exp(res18c$summary.fitted.values[igr,1])
ee

p=ee/(1+ee)
p

meshcoords2<-as.data.frame(meshcoords)

meshcoords2$p<-p

pred<-meshcoords2
pred

head(pred)
head(test)
head(training)
############ MERGE ENTRE PUNTOS (vecinos) ###############

library(RANN)

closest <-nn2(pred[, 1:2],test[, 8:9], 1)
closest
closest<-as.data.frame(closest)
head(closest)



closest$id_test<-rep(1:23480)

pred$nn.idx<-rep(1:1421)

test$id_test<-rep(1:23480)



union<-merge(test, closest, by= "id_test")
union
union2<-merge(union, pred, by= "nn.idx")
union2
union3<-aggregate(union2$nn.dists, by=list(union2$nn.idx), FUN=min, na.rm=TRUE)
union3

names (union3) = c("nn.idx", "nn.dists")

uniones<-merge(union3, union2, by=c("nn.idx","nn.dists"), all.x=T)
uniones


m=uniones
head(m)


##After run the model and have the predicton

#you validate it with the test dataset that you have do before

cor(m$p,m$PA)

#windows(20,20)
plot(m$p,m$PA)


ID <-  as.numeric(dimnames(m)[[1]])  # Id

verif<-as.data.frame(cbind(ID,m$PA,m$p))  #los uno en una tabla en forma: ID, Observados, predichos


library(PresenceAbsence)

tiff("C:/Users/Nerea/Desktop/auc.tif", res = 600, units = "cm", width = 20, height = 20)
tiff("C:/Users/Nerea/Desktop/plots.tif", res = 600, units = "cm", width = 20, height = 20)


dev.off()

windows(20,20)
auc.roc.plot(verif, main="AUC plot")
windows(20,20)
presence.absence.summary(verif,which.model=1)

auc(verif)

threshold<-optimal.thresholds(verif,threshold=101,which.model=1,opt.methods=4);threshold

CMX<-cmx(verif,threshold=threshold[,2]);CMX

Kappa(CMX)

sensitivity(CMX)

specificity(CMX)

threshold


prob.mean.raster<-raster(list(x = pgrid$x,
                              y = pgrid$y,
                              z = ms.m2))



pred= as.data.frame(cbind(coordinates(prob.mean.raster),getValues(prob.mean.raster)))
head(pred)
pred.mean.daily <- aggregate(V3 ~  x + y, pred, function(x) c(Mean = mean(x), CV = cv(x)))
head(pred.mean.daily)


#write.csv(pred,"C:/Users/Nerea/Desktop/p1.csv")

#require(plyr)
# make your data frame
#I<-pred.mean.daily

# make an adjustment grid
#k<-expand.grid(c(0.5,0.5),c(0.5,0.5),0)

# use plyr:ddply() to expand out each entry into the correponding 4 entries
#new_I<-ddply(I,.(y,x),function(x)as.list(x)+k)
#colnames(new_I)<-c("lon","lat","newlat","newlon","tmp")

#head(new_I)
#dim(new_I)
#write.csv(new_I,"C:/Users/Nerea/Desktop/p2.csv")


#rast <- raster(xmn= -150, ymn= -70, xmx = -25, ymx = 35, resolution = 0.5,
#crs = '+proj=utm +zone=33 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ')

#rast
#head(pred.mean.daily)
coordinates(pred.mean.daily) = ~x + y #set spatial coordinates to create a Spatial object

plot(pred.mean.daily)



x.range <- as.numeric(c(-50.25, 15.25)) # map extent
y.range <- as.numeric(c(-25.25, 25.25))
grd <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = 0.5), 
                   y = seq(from = y.range[1], to = y.range[2], by = 0.5)) # expand points to grid
coordinates(grd) <- ~x + y 
gridded(grd) <- TRUE
plot(grd)
head(grd)
idw <- idw(V3.Mean ~ 1, locations = pred.mean.daily, newdata = grd) # apply idw model for the data
idw.output = as.data.frame(idw)
head(idw.output)
write.csv(idw.output,"C:/Users//Nerea//Desktop/pff.csv")


