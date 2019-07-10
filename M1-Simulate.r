#### R SCRIPT FOR SIMULATING SUMMARY STATISTICS UNDER ABC MODEL 1 - Constant dispersal across lanscape
set.seed(100)

library(fields)
library(igraph)
library(qgraph)
library(msm)
library(gtools)
library(abc)

##### POPULATION TRENDS AND RECRUITMENT
pop.trends<-read.csv("d:/CONTAIN/Mink Border Example/Emulator/Data Exporting/Sim0-PopTrends.csv", header=TRUE, sep=",")

summary(pop.trends)

pop.trends$Year<-(pop.trends$Year-min(pop.trends$Year))+1 ### For whenever the first year in the simulations is not = 1

### Sort the data by cell id and year
pop.trends<-pop.trends[order(pop.trends$id.cell, pop.trends$Year),] 

#### TRANSFORM INTO MATRIX
juv<-matrix(pop.trends$NJuvs, ncol=max(pop.trends$Year), nrow=max(pop.trends$id.cell), byrow=TRUE)      #### Number of juveniles per year (columns) and cell (rows)

ad<-matrix(pop.trends$NInd_stage1, ncol=max(pop.trends$Year), nrow=max(pop.trends$id.cell), byrow=TRUE)        #### Number of adults per year (columns) and cell (rows)

##### Occasions and others
n.occ<-max(pop.trends$Year)         ### Number of years

n.cell<-max(pop.trends$id.cell)     ### Number of cells

cell.id<-pop.trends$id.cell         ### Id of each cell

#### DATA GENERATING FUNCTION

mink.m1<-function(x){

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

        ##### CREATE EMPTY MATRICES
        juv.est<-ad.est<-matrix(0, ncol=ncol(juv), nrow=nrow(juv))

        ### Initial abundances of adults and juveniles

        juv.est[ , 1]<-rpois(n.cell, juv.init)

        ad.est[, 1]<-rpois(n.cell, ad.init)


    #### Loop to generate adult and juvenile population sizes per cell and year

    for (i in 1:n.cell){

        for (j in 2:n.occ){

			    #### Immigration & litter size

                immi.ad<-rbinom(1, size=sum(ad.est[-i, j-1]), prob=(1-pd.ad)^(n.cell-2))               ##### Total adults dispersing in the landscape & arriving in cell i

                immi.juv<-rbinom(1, size=sum(juv.est[-i, j-1]), prob=(1-pd.juv)^(n.cell-2))            #### Total juveniles dispersing in the landscape & arriving in cell i

                tot.immi<-immi.ad+immi.juv                                                ##### Total immigration into cell i

                ### Mean litter size
                #f<-rpois(1, mean.litter)

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
pb<-txtProgressBar(min = 0, max = n.its, style = 3)

### Empty matrices to store the results of the simulations
pred.vector<-matrix(0, nrow= n.occ*4*2, ncol=n.its)

prior.vector<-matrix(0, nrow=7, ncol=n.its)

								      
#### Timing the function
start.time<-Sys.time()

### Run the simulations
								      
for (i in 1:n.its){
     
    prior.vector[, i]<-c(runif(1, 0.5, 1), runif(1, 0.5, 1),
                    runif(1, 1, 3), runif(1, 0, 0.5), runif(1, 0.2, 0.8), runif(1, 1, 5), runif(1, 1, 5))       #### Priors

    pred.vector[ , i]<-tryCatch(mink.m1(prior.vector[, i]), error=function(e) NA)                            ### Storing the summary statistics

    setTxtProgressBar(pb, i)

}

close(pb)

end.time<-Sys.time()

time.taken<-end.time-start.time

time.taken

#### Identifying combinations of prior values that produced NAs in the summary statistics

sum(is.na(colMeans(pred.vector)))/n.its                     ### Proportion not plausible

id.na<-which(is.na(colMeans(pred.vector)))

pred2<-pred.vector[ , -id.na]                                  ### Exclude simulations that produce NAs

ncol(pred2)

prior2<-prior.vector[ , -id.na]

##### Exporting

write.table(pred2, "d:/CONTAIN/Mink Border Example/ABC Full Modelling/Outputs/M1-SumStats1.csv", sep=",")

write.table(prior2, "d:/CONTAIN/Mink Border Example/ABC Full Modelling/Outputs/M1-PriorValues.csv", sep=",")

