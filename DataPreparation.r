#### R script for converting and exporting the Range Shifter output into a comma-separated-value file to use in the parameter simulation functions

set.seed(2012)

sim0<-read.delim("d:/CONTAIN/Mink Border Example/Range Shifter Models/Outputs/Sim0_Pop.txt", header=TRUE)   #### Read the populaiton file

summary(sim0)

###### Extract one of the replicated runs of the RS model

rep1<-sim0[sim0$Rep==0, ]

#### Reorder the data, modify years, extract coordinates, and assign individual cell values (from the raster and derived from the x and y coordinates)
rep1$Year<-rep1$Year+1

coord<-table(data.frame(x=rep1$x, y=rep1$y))

patch.id<-data.frame(id.cell=c(1:(ncol(coord)*nrow(coord))), y=rep(as.numeric(colnames(coord)), nrow(coord)), x=rep(as.numeric(rownames(coord)), each=ncol(coord)))


### Create new dataset
tot.rep1<-merge(rep1, patch.id, by=c("x", "y"))

id.cell<-as.factor(tot.rep1$id.cell)            ### Raster cell ID

tot.rep1$id.cell<-as.numeric(id.cell)

n.cell<-max(tot.rep1$id.cell)                     ### Number of cells inhabited by the species during the RS simulatoins

n.cell

###### change variable NInd to total population size after breeding to avoid confusion

tot.rep1$NInd<-tot.rep1$NInd_stage1+tot.rep1$NJuvs

summary(tot.rep1)

#### Export
write.table(tot.rep1, "d:/CONTAIN/Mink Border Example/ABC Full Modelling/Input data/Sim0-PopTrends.csv", sep=",")

