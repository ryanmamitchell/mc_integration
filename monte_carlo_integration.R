library(randtoolbox)
library(gtools)
library(DiceDesign)
library(plotly)
library(cubature)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)



N_p=100000 #Number of points to generate in each pseudo-random sequence


### Generating pseudo-random sequences ###
sobol=sobol(n=N_p, dim=10)
halton=halton(n=N_p, di=10)
faure=runif.faure(n=N_p, dim=10)$design


### Function applied to each sequence ###
basic_function<-function(seq) {abs(3*seq-1)} 


### Interactive plot of the function in 3D domain ###
points=cbind(seq(0, 1, length=1000), seq(0, 1, length=1000))
points=cbind(seq(-10, 10, length=1000), seq(-10, 10, length=1000))
z=outer(basic_function(points[,1]), basic_function(points[,2]))
t<-list(family="Palatino", size=13, color="black")
Fig<-plot_ly(x=points[,1], y=points[,2], z=z, showscale = FALSE, colors=colorRampPalette(brewer.pal(11,"Spectral"))(41))%>%add_surface
Fig<-Fig %>% layout(scene=list( 
  xaxis=list(title="X", family="Palatino", Size=15, color="black"),
  yaxis=list(title="Y", family="Palatino", Size=15, color="black"),
  zaxis=list(title="Z", family="Palatino", Size=15, color="black")), font=t)
Fig



### Quasi Monte Carlo integration (for pseudo-random sequences) ###
integral_QMC = function(sequence, k) {
  space <- apply(sequence, 2, basic_function)
  dep_var=cbind(space[,1], matrix(nrow=nrow(space), ncol=k-1))
  for (i in 1:(k-1)) {            
    dep_var[,i+1] = dep_var[,i] * space[, i+1]                  #Calculating the resulting sequence (dependent variable)
  }                                             
  span = function (seq) {                                       #Finding the span of each sequence
    max(seq) - min(seq)
  }                                              
  seq_span = apply(sequence, 2, span)
  hypervolume=prod(seq_span)                                    #Calculating the volume of a hypercube
  integration_result_QMC = mean(dep_var[,k]) * hypervolume      #Integration: the mean of our resulting sequence multiplied by the volume of a hypercube 
}



#### Standard Monte Carlo integration (for random sequence) ###
rep=1000 #Number of repetitions


iteration_MC=vector(length=length(rep))


#Monte Carlo integral takes number of dimensions, number of generated random points
#and number of repetitions as arguments

integral_MC = function(k, N, reps) {
  for (r in 1:reps){
    unif=matrix(nrow=N, ncol=k)                                    #Generating new uniformly distributed sequence in each run
    for (i in 1:k){
      unif[,i]=runif(n=N)
      cbind(unif[,i-1], unif[,i])
    }
    space <- apply(unif, 2, basic_function)
    dep_var=cbind(space[,1], matrix(nrow=nrow(space), ncol=k-1))
    for (i in 1:(k-1)) {            
      dep_var[,i+1] = dep_var[,i] * space[, i+1]                 
    }                                             
    span = function (seq) {                                       
      max(seq) - min(seq)
    }   
    seq_span = apply(unif, 2, span)
    hypervolume=prod(seq_span)                                    
    iteration_MC[r] = mean(dep_var[,k]) * hypervolume
  }
  integration_result_MC=mean(iteration_MC)
}



### Integration results for different sequences and dimensions ###
dimensions=c(5, 7, 10)

sequences_QMC=c("sobol", "halton", "faure")
sequence_MC=c("unif")


integration_results_QMC=data.frame(matrix(nrow=length(dimensions), ncol=length(sequences_QMC)))
colnames(integration_results_QMC)=sequences_QMC
rownames(integration_results_QMC)=dimensions


for (s in sequences_QMC) {
  for (row in 1:length(dimensions)){
    integration_results_QMC[row, s] = integral_QMC(get(s), dimensions[row])}
}

integration_results_QMC



### Estimation results for different number of generated points ###
N_range=c(100, 1000, 5000 ,10000, seq(50000, 500000, by=50000))


## Pseudo-random ##
integration_N_5D_QMC=data.frame(matrix(ncol=length(sequences_QMC)))
integration_N_7D_QMC=data.frame(matrix(ncol=length(sequences_QMC)))
integration_N_10D_QMC=data.frame(matrix(ncol=length(sequences_QMC)))

colnames(integration_N_5D_QMC)=sequences_QMC
colnames(integration_N_7D_QMC)=sequences_QMC
colnames(integration_N_10D_QMC)=sequences_QMC


#5D
for (s in sequences_QMC) {
  for (n in N_range){
    sobol=sobol(n=n, dim=5)
    halton=halton(n=n, dim=5)
    faure=runif.faure(n=n, dim=5)$design
    integration_N_5D_QMC[n, s]=integral_QMC(get(s)[1:n,], 5)
  }
}

integration_N_5D_QMC=na.omit(integration_N_5D_QMC)


#7D
for (s in sequences_QMC) {
  for (n in N_range){
    sobol=sobol(n=n, dim=7)
    halton=halton(n=n, dim=7)
    faure=runif.faure(n=n, dim=7)$design
    integration_N_7D_QMC[n, s]=integral_QMC(get(s)[1:n,], 7)
  }
}

integration_N_7D_QMC=na.omit(integration_N_7D_QMC)


#10D
for (s in sequences_QMC) {
  for (n in N_range){
    sobol=sobol(n=n, dim=10)
    halton=halton(n=n, dim=10)
    faure=runif.faure(n=n, dim=10)$design
    integration_N_10D_QMC[n, s]=integral_QMC(get(s)[1:n,], 10)
  }
}

integration_N_10D_QMC=na.omit(integration_N_10D_QMC)



## Random ##
integration_N_5D_MC=data.frame(matrix(ncol=length(sequence_MC)))
integration_N_7D_MC=data.frame(matrix(ncol=length(sequence_MC)))
integration_N_10D_MC=data.frame(matrix(ncol=length(sequence_MC)))

colnames(integration_N_5D_MC)=sequence_MC
colnames(integration_N_7D_MC)=sequence_MC
colnames(integration_N_10D_MC)=sequence_MC


#5D
for (n in N_range) {
  integration_N_5D_MC[n,]=integral_MC(5, n, rep)
}

integration_N_5D_MC=na.omit(integration_N_5D_MC)


#7D
for (n in N_range) {
  integration_N_7D_MC[n,]=integral_MC(7, n, rep)
}

integration_N_7D_MC=na.omit(integration_N_7D_MC)

#10D
for (n in N_range) {
  integration_N_10D_MC[n,]=integral_MC(10, n, rep)
}

integration_N_10D_MC=na.omit(integration_N_10D_MC)




### True value of the integral ###
integral_num_5D <- function(x) {abs(3*x[1]-1) * abs(3*x[2]-1)*abs(3*x[3]-1)*abs(3*x[4]-1)*abs(3*x[5]-1) }
true_value_5D=adaptIntegrate(integral_num_5D, lowerLimit = rep(0, 5), upperLimit = rep(1, 5))$integral


integral_num_7D <- function(x) {abs(3*x[1]-1) * abs(3*x[2]-1)*abs(3*x[3]-1) * abs(3*x[4]-1) * abs(3*x[5]-1)*abs(3*x[6]-1) * abs(3*x[7]-1)}
true_value_7D=adaptIntegrate(integral_num_7D, lowerLimit = rep(0, 7), upperLimit = rep(1, 7))$integral


integral_num_10D <- function(x) {abs(3*x[1]-1) * abs(3*x[2]-1)*abs(3*x[3]-1) * abs(3*x[4]-1) * abs(3*x[5]-1)*abs(3*x[6]-1) * abs(3*x[7]-1)*abs(3*x[8]-1) * abs(3*x[9]-1) * abs(3*x[10]-1)}
#true_value10D=adaptIntegrate(integral_num_10D, lowerLimit = rep(0, 10), upperLimit = rep(1, 10))$integral
true_value_10D=0.16150567521


true_values=c(true_value_5D, true_value_7D, true_value_10D)



#### Combining results from QMC and MC integration into a single dataframe ###
integration_N_5D=cbind(integration_N_5D_QMC, integration_N_5D_MC)
integration_N_7D=cbind(integration_N_7D_QMC, integration_N_7D_MC)
integration_N_10D=cbind(integration_N_10D_QMC, integration_N_10D_MC)



#### Integration error ###

#5D
integration_error5D=data.frame(matrix(nrow=nrow(integration_N_5D), ncol=ncol(integration_N_5D))) 
colnames(integration_error5D)=colnames(integration_N_5D)

for (row in 1:nrow(integration_error5D)) {
  integration_error5D[row,]=abs(integration_N_5D[row,]-true_values[1])/true_values[1]
  integration_error5D
}


#7D
integration_error7D=data.frame(matrix(nrow=nrow(integration_N_7D), ncol=ncol(integration_N_7D))) 
colnames(integration_error7D)=colnames(integration_N_7D)                                                      

for (row in 1:nrow(integration_error7D)) {
  integration_error7D[row,]=abs(integration_N_7D[row,]-true_values[2])/true_values[2]          
}                                                                                                                                    


#10D
integration_error10D=data.frame(matrix(nrow=nrow(integration_N_10D), ncol=ncol(integration_N_10D)))
colnames(integration_error10D)=colnames(integration_N_10D)

for (row in 1:nrow(integration_error10D)) {
  integration_error10D[row,]=abs(integration_N_10D[row,]-true_values[3])/true_values[3]           
}                                                                                               



### Absolute error  ### 
SE5D=abs(matrix(rep(true_value_5D, 44), nrow=length(N_range), ncol=ncol(integration_N_5D))
      -(integration_N_5D))



SE7D=abs(matrix(rep(true_value_7D, 44), nrow=length(N_range), ncol=ncol(integration_N_7D))
      -(integration_N_7D))



SE10D=abs(matrix(rep(true_value_10D, 44), nrow=length(N_range), ncol=ncol(integration_N_10D))
       -(integration_N_10D))



##### GRAPHS #####


#Generic QMC and Unif 
library(cowplot)

m <- 1000
prime <- c(2,3)

#executing sequences
halton_gen <- halton(m,2)
halton_gen <- as.data.frame(halton_gen)
colnames(halton_gen) <- c("X1","X2")


faure_gen <- (runif.faure(m,2))
faure_gen <- faure_gen[1]
faure_gen <- as.data.frame(faure_gen)
colnames(faure_gen) <- c("X1","X2")

sobol_gen <- sobol(m,2,1,prime)
sobol_gen <- as.data.frame(sobol_gen)
colnames(sobol_gen) <- c("X1","X2")


r1 <- runif(m,0,1)
r2 <- runif(m,0,1)
random_gen <- data.frame(r1,r2)
random_gen
colnames(random_gen) <- c("X1","X2")


#plotting graphs

random_plot_gen <-  ggplot(data = random_gen) + geom_point(mapping = aes(x = X1, y = X2), color='steelblue',size=.5) + ggtitle("Uniform") +
  theme(plot.title = element_text(family="Palatino",face="bold",color='black',hjust=0.5)) +
  theme(axis.title.x = element_text(family="Palatino",color='black')) +
  theme(axis.title.y = element_text(family = "Palatino",color='black')) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

halton_plot_gen <- ggplot(data = halton_gen) + geom_point(mapping = aes(x = X1, y = X2,),color='steelblue',size=.5) + ggtitle("Halton") +
  theme(plot.title = element_text(family="Palatino",face="bold",color='black',hjust=0.5)) +
  theme(axis.title.x = element_text(family="Palatino",color='black')) +
  theme(axis.title.y = element_text(family = "Palatino",color='black')) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

faure_plot_gen <- ggplot(data = faure_gen) + geom_point(mapping = aes(x = X1, y = X2),color='steelblue',size=.5) + ggtitle("Faure") +
  theme(plot.title = element_text(family="Palatino",face="bold",color='black',hjust=0.5)) +
  theme(axis.title.x = element_text(family="Palatino",color='black')) +
  theme(axis.title.y = element_text(family = "Palatino",color='black')) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

sobol_plot_gen <- ggplot(data = sobol_gen) + geom_point(mapping = aes(x = X1, y = X2),color='steelblue',size=.5) + ggtitle("Sobol") +
  theme(plot.title = element_text(family="Palatino",face="bold",color='black',hjust=0.5)) +
  theme(axis.title.x = element_text(family="Palatino",color='black')) +
  theme(axis.title.y = element_text(family = "Palatino",color='black')) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
plot_grid(random_plot_gen,halton_plot_gen,faure_plot_gen,sobol_plot_gen, labels =NULL)



#Final Result Graphs


# adding dimension number and consolidating SE data frames

SE5D$dim <- c(rep("5",14))
SE7D$dim <- c(rep("7",14))
SE10D$dim <- c(rep("10",14))

g1 <- rbind(SE5D,SE7D,SE10D)

g1$N <- c(rep(N_range,3))

g5 <- g1[4:14,]
g7 <- g1[18:28,]
g10 <- g1[32:42,]
g2 <- rbind(g5,g7,g10)



#adding log2 scale
g1$logN <- log2(g1$N)


#plots by dimension
# s = 5
colors <- c("Halton" = "steelblue","Faure" = "darkred", "Sobol" = "darkgreen", "Unif" = "black")
plot5 <- ggplot(data = g2, aes(x = logN)) + 
  
  geom_line(data = subset(g2,dim==5),mapping=aes(y = halton, color = "Halton")) +
  geom_line(data = subset(g2,dim==5),mapping=aes(y = faure, color = "Faure")) +
  geom_line(data = subset(g2,dim==5),mapping=aes(y = sobol, color = "Sobol")) +
  geom_line(data = subset(g2,dim==5),mapping=aes(y = unif, color = "Unif")) +
  
  
  #Title + theme
  ggtitle("Sequences Absolute Error", subtitle = "S = 5") +
  labs(x = " Log2 of N", y= "Absolute Error", color = "Sequences") + 
  scale_color_manual(values = colors, breaks = c("Halton","Faure","Sobol", "Unif")) +
  
  theme(plot.title = element_text(family="Palatino",face="bold",color='black',hjust=0.5)) +
  theme(plot.subtitle = element_text( family="Palatino",hjust = 0.5)) +
  theme(axis.title.x = element_text(family="Palatino",color='black')) +
  theme(axis.title.y = element_text(family = "Palatino",color='black')) +
  theme(legend.title = element_text(family = "Palatino")) +
  theme(legend.position = "bottom")

#S7
plot7 <- ggplot(data = g2, aes(x = logN)) + 
  
  geom_line(data = subset(g2,dim==7),mapping=aes(y = halton, color = "Halton")) +
  geom_line(data = subset(g2,dim==7),mapping=aes(y = faure, color = "Faure")) +
  geom_line(data = subset(g2,dim==7),mapping=aes(y = sobol, color = "Sobol")) +
  geom_line(data = subset(g2,dim==7),mapping=aes(y = unif, color = "Unif")) +
  
  
  #Titles and theme
  ggtitle("Sequences Absolute Error", subtitle = "S = 7") +
  labs(x = " Log2 of N", y= "Absolute Error", color = "Sequences") + 
  scale_color_manual(values = colors, breaks = c("Halton","Faure","Sobol", "Unif")) +
  
  theme(plot.title = element_text(family="Palatino",face="bold",color='black',hjust=0.5)) +
  theme(plot.subtitle = element_text( family="Palatino",hjust = 0.5)) +
  theme(axis.title.x = element_text(family="Palatino",color='black')) +
  theme(axis.title.y = element_text(family = "Palatino",color='black')) +
  theme(legend.title = element_text(family = "Palatino")) +
  theme(legend.position = "bottom")


#S10
plot10 <- ggplot(data = g2, aes(x = logN)) + 
  
  geom_line(data = subset(g2,dim==10),mapping=aes(y = halton, color = "Halton")) +
  geom_line(data = subset(g2,dim==10),mapping=aes(y = faure, color = "Faure")) +
  geom_line(data = subset(g2,dim==10),mapping=aes(y = sobol, color = "Sobol")) +
  geom_line(data = subset(g2,dim==10),mapping=aes(y = unif, color = "Unif")) +
  
  
  #Title and theme
  ggtitle("Sequences Absolute Error", subtitle = "S = 10") +
  labs(x = " Log2 of N", y= "Absolute Error", color = "Sequences") + 
  scale_color_manual(values = colors, breaks = c("Halton","Faure","Sobol", "Unif")) +
  
  theme(plot.title = element_text(family="Palatino",face="bold",color='black',hjust=0.5)) +
  theme(plot.subtitle = element_text( family="Palatino",hjust = 0.5)) +
  theme(axis.title.x = element_text(family="Palatino",color='black')) +
  theme(axis.title.y = element_text(family = "Palatino",color='black')) +
  theme(legend.title = element_text(family = "Palatino")) +
  theme(legend.position = "bottom")


ggarrange(plot5,plot7,plot10, labels = NULL, ncol = 3, nrow = 1)
