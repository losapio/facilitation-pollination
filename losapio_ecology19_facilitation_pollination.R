### code by Gianalberto Losapio (losapiog@ethz.ch gianalbertolosapio@gmail.com)
### Citation: Plant interactions shape pollination networks via nonadditive effects
### Losapio et al. 2019. Ecology
library(vegan)
library(effects)
library(bipartite)
library(igraph)
library(nlme)
library(betalink)

load("losapio_ecology19_facilitation_pollination.RData")

# This database includes the following items

insect.record	# data of plant-insect interactions
additive		# summary data
				# tr = treatment: e = additive; o = control
				# fs = foundation species: h = Are; s = Hor
				# microh = microhabitat: a = alone; o = open; w = together
spint			# summary data of species interacitons in shared networks

snet			# networks of each treatments combination
bip.are.obs		# network with Arenaria/control
bip.are.exp		# network with Arenaria/additive
bip.hor.obs		# network with Hormathophylla/control
bip.hor.exp		# network with Hormathophylla/additive

null.are.exp	# null models of Arenaria/control
null.are.obs	# null models of Arenaria/additive
null.hor.obs	# null models of Hormathophylla/control
null.hor.exp	# null models of Hormathophylla/additive

nst.are.obs.tot # nestedness of Arenaria/control
nst.are.exp.tot # nestedness of Arenaria/additive
nst.hor.obs.tot # nestedness of Hormathophylla/control
nst.hor.exp.tot # nestedness of Hormathophylla/additive

nulls			# nestedness of null models

pval.l.are.obs	# p-value nestedness of Arenaria/control
pval.l.are.exp	# p-value nestedness of Arenaria/additive
pval.l.hor.obs	# p-value nestedness of Hormathophylla/control
pval.h.hor.exp	# p-value nestedness of Hormathophylla/additive

z.are.obs		# z-score nestedness of Arenaria/control
z.are.exp		# z-score nestedness of Arenaria/additive
z.hor.obs		# z-score nestedness of Hormathophylla/control
z.hor.exp		# z-score nestedness of Hormathophylla/additive

mod.addi.wdiv	# mixed-effects model of weighted pollinator diversity 
mod.visits		# mixed-effects model of weighted pollinator visitation rate
mod.intd		# mixed-effects model of shared-species interaction diversity 

### workflow
# nestedness & network str
# example with are.obs

web=ifelse(empty(bip.are.obs) >0,1,0)

## plants
nst = NA
z=0

for(i in 1:(nrow(web)-1)){for(j in (i+1):nrow(web)){
		print(c(i,j))
		z=z+1
		num=0
		for(k in 1:ncol(web)){
			if(web[i,k]==web[j,k]&web[i,k]==1) num=num+1
			}
		nst[z]= num/min(sum(web[i,]),sum(web[j,]))
}}

nst.are.obs.pl = mean(nst)

## poll
nst = NA
z=0

for(i in 1:(ncol(web)-1)){for(j in (i+1):ncol(web)){
		print(c(i,j))
		z=z+1
		num=0
		for(k in 1:nrow(web)){
			if(web[k,i]==web[k,j]&web[k,i]==1) num=num+1
			}
		nst[z]= num/min(sum(web[,i]),sum(web[,j]))
}}

nst.are.obs.in = mean(nst)

## tot
nst.are.obs.tot = (nst.are.obs.pl+nst.are.obs.in)/2

### observed P value
ip <- rep(0,nm)
for(i in 1:nm){
	if(nst.are.obs.tot > nulls$nst.are.obs.tot1[i]) ip[i] <- 1
}

pval.h.are.obs <- 1-(sum(ip)/(nm+1))

ip <- rep(0,nm)
for(i in 1:nm){
	if(nst.are.obs.tot < nulls$nst.are.obs.tot1[i]) ip[i] <- 1
}

pval.l.are.obs <- 1-(sum(ip)/(nm+1))

#z-score
z.are.obs <- (nst.are.obs.tot-mean(nulls$nst.are.obs.tot1))/sd(nulls$nst.are.obs.tot1)
z.are.exp <- (nst.are.exp.tot-mean(nulls$nst.are.exp.tot1))/sd(nulls$nst.are.exp.tot1)
z.hor.obs <- (nst.hor.obs.tot-mean(nulls$nst.hor.obs.tot1))/sd(nulls$nst.hor.obs.tot1)
z.hor.exp <- (nst.hor.exp.tot-mean(nulls$nst.hor.exp.tot1))/sd(nulls$nst.hor.exp.tot1)

### stats

# diversity of insects: $wdiv
mod.addi.wdiv = lme(fixed=div ~ tr*fs,
	random= ~ 1| plot,
	method="REML",
	data=additive)
anova(mod.addi.wdiv, test="Chisq")
summary(mod.addi.wdiv)

lsmeans(mod.addi.wdiv, pairwise ~tr*fs)

# visitation rate: $infl
mod.visits = lme(fixed=infl ~ tr * fs,
	random= ~ 1| sp/plot,
	method="REML",
	data=visits)
anova(mod.visits, test="Chisq")
summary(mod.visits)

lsmeans(mod.visits, pairwise ~tr*fs)

# interaction diversity of species for shared networks: $intd

mod.intd = lme(fixed=intd ~ gui + tf * fs,
	random= ~ 1| sp,
	method="REML",
	data=spint)
anova(mod.intd)
summary(mod.intd)

lsmeans(mod.intd, pairwise ~gui + fs*tf)

####### plots
## diversity of insects
ci1 <- summary(lsmeans(mod.addi.wdiv, ~tr*fs))

plot(0,0,type="n",xlim=c(1,4), ylim=c(min(ci1$lsmean - ci1$SE),max(ci1$lsmean + ci1$SE)), xaxt="n",xlab="",ylab="",yaxt="n")
axis(2,at=c(.8,1,1.2,1.4,1.6), las=1)

segments(1, ci1$lsmean[1]-ci1$SE[1], 1, ci1$lsmean[1]+ci1$SE[1])
segments(2, ci1$lsmean[2]-ci1$SE[2], 2, ci1$lsmean[2]+ci1$SE[2])
segments(3, ci1$lsmean[3]-ci1$SE[3], 3, ci1$lsmean[3]+ci1$SE[3])
segments(4, ci1$lsmean[4]-ci1$SE[4], 4, ci1$lsmean[4]+ci1$SE[4])

segments(1, ci1$lsmean[1], 2, ci1$lsmean[2], lty=2)
segments(3, ci1$lsmean[3], 4, ci1$lsmean[4], lty=2)

points(1, ci1$lsmean[1], cex=1.5, pch=19, col="black", bg="white")
points(2, ci1$lsmean[2], cex=1.5, pch=17, col="black", bg="white")
points(3, ci1$lsmean[3], cex=1.5, pch=19, col="black", bg="white")
points(4, ci1$lsmean[4], cex=1.5, pch=17, col="black", bg="white")

# interaction diversity of species

ci1 <- summary(lsmeans(mod.intd, ~tr*fs))

plot(0,0,type="n",ylim=c(1,4), xlim=c(min(ci1[,1]),max(ci1[,3])), xaxt="n",xlab="",ylab="",yaxt="n")
axis(1,at=c(-1,-0.5,0,0.5), las=1)

segments(ci1[1,1], 4, ci1[1,3], 4)
segments(ci1[2,1], 3, ci1[2,3], 3)
segments(ci1[3,1], 2, ci1[3,3], 2)
segments(ci1[4,1], 1, ci1[4,3], 1)

points(ci1[1,2], 4)
points(ci1[2,2], 3)
points(ci1[3,2], 2)
points(ci1[4,2], 1)

# visitation rate

ci1 <- summary(lsmeans(mod.visits, ~tr*fs))

plot(0,0,type="n",xlim=c(1,4), ylim=c(min(ci1$lsmean - ci1$SE),max(ci1$lsmean + ci1$SE)), xaxt="n",xlab="",ylab="",yaxt="n")
axis(2,at=c(log(0.001), log(0.002), log(0.005), log(0.01)), las=1)

segments(1, ci1$lsmean[1]-ci1$SE[1], 1, ci1$lsmean[1]+ci1$SE[1])
segments(2, ci1$lsmean[2]-ci1$SE[2], 2, ci1$lsmean[2]+ci1$SE[2])
segments(3, ci1$lsmean[3]-ci1$SE[3], 3, ci1$lsmean[3]+ci1$SE[3])
segments(4, ci1$lsmean[4]-ci1$SE[4], 4, ci1$lsmean[4]+ci1$SE[4])

segments(1, ci1$lsmean[1], 2, ci1$lsmean[2], lty=2)
segments(3, ci1$lsmean[3], 4, ci1$lsmean[4], lty=2)

points(1, ci1$lsmean[1], cex=1.5, pch=19, col="black", bg="white")
points(2, ci1$lsmean[2], cex=1.5, pch=17, col="black", bg="white")
points(3, ci1$lsmean[3], cex=1.5, pch=19, col="black", bg="white")
points(4, ci1$lsmean[4], cex=1.5, pch=17, col="black", bg="white")

dev.off()

# nestedness and density lines

plot(0,0,type="n",xlim=c(0.19,0.38), ylim=c(0,380), xaxt="n",xlab="",ylab="",yaxt="n")
axis(1,at=c(0.2,0.25,0.3,0.35,0.4), las=1)

#
lines(density(nulls$nst.are.exp.tot1, adjust=2), col="red", lwd=3)
abline(v=nst.are.exp.tot, col="red", lwd=3)

#
lines(density(nulls$nst.are.obs.tot1, adjust=2), col="green", lwd=3)
abline(v=nst.are.obs.tot, col="green", lwd=3)

#
lines(density(nulls$nst.hor.exp.tot1, adjust=2), col="blue", lwd=3)
abline(v=nst.hor.exp.tot, col="blue", lwd=3)

#
lines(density(nulls$nst.hor.obs.tot1, adjust=2), col="black", lwd=3)
abline(v=nst.hor.obs.tot, col="black", lwd=3)

### SI
# plant species similarity
plant.spmat

plantmat <- as.matrix(plant.spmat[,c(5,9,10,17,18,19,20,23)])
plantmatcov<-plant.spmat[,1:2]
plantmatcov$treatm<-paste(plantmatcov$foundation, plantmatcov$microh, sep="")

plantmatcov$treatm[plantmatcov$treatm=="ho"]<-"At1"
plantmatcov$treatm[plantmatcov$treatm=="hw"]<-"At2"
plantmatcov$treatm[plantmatcov$treatm=="so"]<-"Hs1"
plantmatcov$treatm[plantmatcov$treatm=="sw"]<-"Hs2"

ord<-metaMDS(plantmat)

plot(ord, disp="sites",type="n")
ordiellipse(ord, treatm, col=2, kind = "ehull", lwd=3) 
points(ord, disp="sites", pch=21, col="black", bg="green", cex=1.3)
ordiellipse(ord, treatm, col=1:4, draw="polygon") 
ordispider(ord, treatm, col=1:4, label = T) 

#