#### R SCRIPT FOR SIMULATING SUMMARY STATISTICS UNDER ABC MODEL 2 - DISPERSAL KERNEL BASED ON DISTANCE
set.seed(100)

library(fields)
library(igraph)
library(qgraph)
library(msm)
library(gtools)
library(abc)
  
##### POPULATION TRENDS AND RECRUIRMENT
pop.trends<-read.csv("d:/CONTAIN/Mink Border Example/Emulator/Data Exporting/Sim0-PopTrends.csv", header=TRUE, sep=",")

pop.trends$Year<-(pop.trends$Year-min(pop.trends$Year))+1       ### For whenever the first year in the simulations is not = 1

summary(pop.trends)

### Sort the data by cell id and year
pop.trends<-pop.trends[order(pop.trends$id.cell, pop.trends$Year),] 

#### TRANSFORM INTO MATRIX
juv<-matrix(pop.trends$NJuvs, ncol=max(pop.trends$Year), nrow=max(pop.trends$id.cell), byrow=TRUE)      #### Number of juveniles per year (columns) and cell (rows)

ad<-matrix(pop.trends$NInd_stage1, ncol=max(pop.trends$Year), nrow=max(pop.trends$id.cell), byrow=TRUE)        #### Number of adults per year (columns) and cell (rows)

##### Occasions and others
n.occ<-max(pop.trends$Year)         ### Number of years

n.cell<-max(pop.trends$id.cell)     ### Number of cells

cell.id<-pop.trends$id.cell         ### Id of each cell

### Coordinates of each cell
cell.coord<-unique(data.frame(x=pop.trends$x, y=pop.trends$y, cell.id=pop.trends$id.cell))

#### Put together the patch centroids
coord.patches<-matrix(c(cell.coord$x, cell.coord$y), nrow=n.cell, ncol=2)

#### Pairwise euclidean distance matrix calculated using function rdist from package fields (times the original cell resolution in metres; 150 metres)
D<-rdist(coord.patches)*150					

#### DATA GENERATING FUNCTION

mink.m2<-function(x){

        #### Survival
	    phi.ad<-x[1]

	    phi.juv<-x[2]

        #### Recruitment
        #### Litter size per female
        f<-x[3]

        ########## Probabilities of immigration
        pd.ad<-x[4]      ##### Adult probabiliy of dispersal

        pd.juv<-x[5]     ### Juvenile probability of dispersal

        #### Initial population size
        juv.init<-x[6]

        ad.init<-x[7] 

        ##### Dispersal kernel
        c<-1/(x[8])

        dk<-exp(-c*D)

  	    diag(dk)<-0

        ### Normalisation of the kernel by row
        norm.par<-rowSums(dk)

        norm.k<-apply(dk, 1, function(x) x/norm.par)

        ##### CREATE EMPTY MATRICES
        juv.est<-ad.est<-matrix(0, ncol=ncol(juv), nrow=nrow(juv))

        juv.est[ , 1]<-rpois(n.cell, juv.init)

        ad.est[, 1]<-rpois(n.cell, ad.init)

    #### Loop to generate adult and juvenile population sizes per cell and year

    for (i in 1:n.cell){

             for (j in 2:n.occ){

			    #### Immigration & litter size

                immi.ad<-rbinom(1, size=sum(ad.est[-i, j-1]), prob=pd.ad*mean(norm.k[i, -i]))     ##### Mean total adults dispersing in the landscape that will arrive in cell i

                immi.juv<-rbinom(1, size=sum(juv.est[-i, j-1]), prob=pd.juv*mean(norm.k[i, -i]))       #### Mean total juveniles dispersing in the landscape that will arrive in cell i

                tot.immi<-immi.ad+immi.juv                                                ##### Total immigration into cell i

                ##### Population sizes
                #### Apparent survival adults
                surv.mink<-rbinom(1, size=ad.est[i, j-1], prob=phi.ad*(1-pd.ad))
                             
                #### Recruitment juveniles
                rec.mink<-rbinom(1, size=juv.est[i, j-1], prob=phi.juv*(1-pd.juv))
    
                ### Adult population
                ad.est[i, j]<-rpois(1, surv.mink+rec.mink+tot.immi)

                #### Number born
                juv.est[i, j]<-rpois(1, ad.est[i, j]*f)
                
                     }
        }

    #### Output - summary statistic per year
    out.sim<-c(colMeans(ad.est), apply(ad.est, 2, sd), apply(ad.est, 2, function(x) quantile(x, 0.0275)),
            apply(ad.est, 2, function(x) quantile(x, 0.975)), colMeans(juv.est), apply(juv.est, 2, sd),
            apply(juv.est, 2, function(x) quantile(x, 0.0275)), apply(juv.est, 2, function(x) quantile(x, 0.975)))

    return(out.sim)

	}

### Define the number of simulations to run
n.its<-20000

### First set a progress bar
pb<-txtProgressBar(min=0, max=n.its, style = 3)

### Empty matrices to store the results of the simulations
pred.vector<-matrix(0, nrow=n.occ*4*2, ncol=n.its)

prior.vector<-matrix(0, nrow=8, ncol=n.its)

#### Timing the function
start.time<-Sys.time()

### Run the simulations
								      
for (i in 1:n.its){
     
    prior.vector[, i]<-c(runif(1, 0.5, 1), runif(1, 0.5, 1),
                    runif(1, 1, 3), runif(1, 0, 0.5), runif(1, 0.2, 0.8), runif(1, 1, 5), runif(1, 1, 5), runif(1, 10000, 50000))  #### Priors

    pred.vector[ , i]<-tryCatch(mink.m2(prior.vector[, i]), error=function(e) NA)                ### Storing the summary statistics

    setTxtProgressBar(pb, i)

}

close(pb)


end.time<-Sys.time()

time.taken<-end.time-start.time

time.taken

#### Identifying combinations of prior values that produced NAs in the summary statistics

sum(is.na(colMeans(pred.vector)))/n.its                     ### Proportion not plausible

id.na<-which(is.na(colMeans(pred.vector)))

pred2<-pred.vector[ , -id.na]                               ### Exclude simulations that produce NAs

ncol(pred2)

prior2<-prior.vector[ , -id.na]

##### Exporting

write.table(pred.vector, "d:/CONTAIN/Mink Border Example/ABC Full Modelling/Outputs/M2-SumStats1.csv", sep=",")

write.table(prior.vector, "d:/CONTAIN/Mink Border Example/ABC Full Modelling/Outputs/M2-PriorValues.csv", sep=",")



