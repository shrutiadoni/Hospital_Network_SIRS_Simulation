library(dplyr)

## loading data
setwd("C:/Users/shrut/OneDrive/Desktop/Shruti/UIC/IDS 564/Project/Epidemics") 

#setwd("D:/UIC/Spring 2020/IDS 564Social Med and Net Analysis/Labs/project/Final File/Hospital") 

node_data <- read.csv("nodes.csv", head=T, as.is=T) 
edge_data <- read.csv("edges.csv", head=T, as.is=T) 

head(node_data)
head(edge_data)

#Day 1 Data

edge_data$time_start <- edge_data$time_start / 20 - 5   # 120 seconds(1:02pm, Mon) -> 1, unit : 20 seconds
edge_data$time_end <- edge_data$time_end / 20 - 5   # 120 seconds -> 1


# Node color based on the node type
node_data$color <- "darkgray"      # PAT
node_data$color[node_data$type == "NUR"] <- "blue"
node_data$color[node_data$type == "MED"] <- "cyan"
node_data$color[node_data$type == "ADM"] <- "green"

# create spell columns for nodes 0 - Inf
# Onset and terminus decide for how long the node will pe present in the network. We have fixed the onset and terminus making all nodes active throughout the time.
node_data$onset <- 0
node_data$terminus <- Inf

# create edge attributes

edge_data$attrs <- NA


for(i in 1 :  nrow(edge_data)){
  n1 <- node_data$type[node_data$Id == edge_data$Source[i]]
  n2 <- node_data$type[node_data$Id == edge_data$Target[i]]
  n3 <- sort(c(n1,n2))
  edge_data$attrs[i] <- paste(n3[1],"-",n3[2], sep = "")
}

head(edge_data)

###combining dynamic edges

# a contact is divided by multiple time units = increase in no. of edges
edge_data[edge_data$attrs == "PAT-PAT",][1:10,] 


test<-edge_data %>% group_by(Source,Target) %>%mutate(temp_id = row_number())      
#row_number() : numbering by same group

test
test<-as.data.frame(test)

test$temp_ing <- 1

for (i in 1:nrow(test)){          
  if ( i < nrow(test) & ((test$temp_id[i] - test$temp_id[i + 1]) == -1 & test$Source[i] == test$Source[i+1] ) ){ 
    if (test$temp_ing[i] != 0 ){start <- i}    
    test$time_end[i] <- NA   
    test$temp_ing[i + 1] <- 0                
  }
  else{ if(test$temp_ing[i] == 0){test$time_end[start] <- test$time_end[i]   
  test$time_end[i] <- NA}
  }
  
}
test

# remove rows with NA in time_end
edge_data <- test[complete.cases(test$time_end), ]
head(edge_data)

node_data$newid <- 1:75 
edge_data$from <- match(edge_data$Source,node_data$Id)
edge_data$to <- match(edge_data$Target,node_data$Id)


# SIRS function
library(igraph)

simulate_sir <- function(network.i, simlength=15, p.t=0.2, display_net=TRUE, removeafter=2, susceptibleafter=3) {
  
  all_time_infected_nodes = list()
  all_time_susceptible_again_nodes = list()
  
  links <-get.edgelist(network.i)
  N<- vcount(network.i)
  time_stats<-list()
  
  # Initialize time stats. 
  # Number of nodes in S, I, or R status in each round of time
  time_stats$infected_t<-rep(1,simlength)
  time_stats$removed_t<-rep(0,simlength)
  #susceptible is total that are not removed or infected
  time_stats$susceptible_t<-rep(N-1,simlength)
  
  infected <- logical(N) # initialize infection status
  susceptible <- rep(TRUE, N) # initialize susceptible status
  removed<-logical(N)
  patientzero <- sample(N,1) # select 'patient zero'
  
  # Initialize a vector that keeps track of the time of infection, susceptible and removed for each node. 
  infected_time<-rep(0, N)
  removed_time<-rep(0, N)
  susceptible_time<-rep(0, N)
  
  #patient zero  
  infected[patientzero] <- TRUE
  susceptible[patientzero] <-FALSE
  
  # Used to count towards a removal; after a certain number of periods, the node will be immune (i.e. removed)
  infected_time[patientzero] <- 1
  
  if (N > 50) {
    V(network.i)$size <- 7
    V(network.i)$label <- ""
  }
  if (display_net) {
    
    fixlayout <- layout.kamada.kawai(network.i)  # store a fixed layout for the graph
    node.colour <- rep("SkyBlue2",N) # initialize node colours (SkyBlue2 is also the default node colour in igraph)
    node.colour[patientzero] <- "red" # infected nodes will be coloured red
    plot(network.i,layout=fixlayout, main="Time = 0", vertex.color=node.colour)
  }
  for (i in 1:simlength) {
    
    # Original spreading mechanism, that did not account for removed nodes
   discordant.links <- which(xor(infected[as.integer(links[,1])],infected[as.integer(links[,2])])) # find the indices of links that connect an infected individual to an uninfected
    
    transmit <- rbinom(length(discordant.links),1,p.t) # determine randomly which of the discordant links transmit the disease
    
    # let me update the infection vector in three steps to make it easier to read:
    transmitter.links <- discordant.links[transmit==1]
    nodes.of.transmitter.links <- unique(as.vector(as.integer(links[transmitter.links,1:2]))) # gets both nodes of the transmitter links into a single vector; unique just filters out repetitions
    
    removed.indices = which(removed==TRUE)
    if(length(removed.indices)>0){
      nodes.of.transmitter.links.new <- which(!nodes.of.transmitter.links %in% removed.indices)
      nodes.of.transmitter.links <- nodes.of.transmitter.links.new
    }
    
    infected[nodes.of.transmitter.links] <- TRUE # here I simply set both nodes to TRUE (although the transmitter already had 'TRUE'). In more complex models, you might want to do a further check here and overwrite only the newly infected nodes.
    susceptible[nodes.of.transmitter.links] <- FALSE
    removed[nodes.of.transmitter.links] <- FALSE
    
    # At some point in this loop, you need to update the number infected, and for advanced lab, number removed and susceptible
    # Increment length of time infected, for those that are infected.
    if(length(all_time_infected_nodes)<1){
      all_time_infected_nodes = nodes.of.transmitter.links
    } else {
      all_time_infected_nodes = append(all_time_infected_nodes, nodes.of.transmitter.links)
    }
    remove <- all_time_infected_nodes[duplicated(all_time_infected_nodes)]
    nodes.of.transmitter.links = nodes.of.transmitter.links [! nodes.of.transmitter.links %in% remove]
    
    infected_time.links = infected_time[infected_time==0]
    infected_time[nodes.of.transmitter.links] <- i

    if(i>removeafter){
      r_time = i - removeafter
      cat("i --> ",i,"removeafter --> ",removeafter,"r_time --> ",r_time,"\n")
      recovered_list <- which(infected_time==r_time)
      cat("Recovered List: ",recovered_list,"\n")
      infected_time[recovered_list] <- 0
      infected[recovered_list] <- FALSE
      removed[recovered_list] <- TRUE
      susceptible_time[recovered_list] <- TRUE
    }
    
    if(i>susceptibleafter){
      s_time = i - susceptibleafter
      susceptible_again_list <- which(susceptible_time==s_time)
      susceptible[susceptible_again_list] <- TRUE
      removed[susceptible_again_list] <- FALSE
      infected[susceptible_again_list] <- FALSE
      
      ## Susceptible Again - different color logic
      if(length(all_time_susceptible_again_nodes)<1){
        all_time_susceptible_again_nodes = susceptible_again_list
      } else {
        all_time_susceptible_again_nodes = append(all_time_susceptible_again_nodes, susceptible_again_list)
      }
    }
    
    time_stats$infected_t[i] = sum(infected, na.rm = TRUE)
    time_stats$removed_t[i] = sum(removed, na.rm = TRUE)
    time_stats$susceptible_t[i] = sum(susceptible, na.rm = TRUE)
    cat("Infected count -->",time_stats$infected_t[i],"\n")
    cat("Recovered count -->",time_stats$removed_t[i],"\n")
    cat("Susceptible count -->",time_stats$susceptible_t[i],"\n")
    cat("S + I + R count -->",sum(time_stats$infected_t[i],time_stats$removed_t[i],time_stats$susceptible_t[i]),"\n")
    
    
    if (display_net) {
      node.colour[infected] <- "red"
      # Make the removed points as yellow, susceptible points as skyblue. Uncomment these for advanced lab.
      node.colour[removed] <- "yellow"
      
      
      ## Susceptible Again - different color logic
      susceptible_again_temp <-logical(N)
      susceptible_again1 = which(susceptible==TRUE)
      susceptible_again.links = which(susceptible_again1 %in% all_time_susceptible_again_nodes)
      susceptible_again_temp[susceptible_again.links] <- TRUE
      node.colour[susceptible_again_temp] <-"black"
      #  cat("susceptible_again_temp", susceptible_again_temp,"\n")
      
      node.colour[susceptible] <-"SkyBlue2"
      input<-readline() # waits for the user to press <ENTER> before proceeding; you need to switch to the console to do this
      
      # TIP: Hit q to break the simulation between renderings of the network
      if (input=="q") {break}
      plot(network.i,layout=fixlayout, main=paste("Time =", i, " out of ", simlength), vertex.color=node.colour)
    }
    
  }  
  # time_stats is a list of three vectors that keeps track of number of infected, removed, and susceptible nodes over each round of time.
  return(time_stats)
}

library(igraph)

#creating graph object using the edge_data
g_hospital<-graph.data.frame(edge_data, directed = TRUE)

#plot(g_hospital)

V(g_hospital)$name<-as.integer(1:vcount(g_hospital))

par(mfrow=c(1,1))

#simulate to display the network simulation
infected_school_tm<-simulate_sir(g_hospital, simlength= 20,  p.t=0.05)

#simulate function with display=FALSE
infected_school_tm<-simulate_sir(g_hospital, simlength= 20,  p.t=0.05, display_net=FALSE)

par(mfrow=c(3,1))
plot(infected_school_tm$infected_t, type="l", col="red", ylab="Infected over time", xlab = "time index")
plot(infected_school_tm$removed_t, type="l", col="black", ylab="Removed over time", xlab = "time index")
plot(infected_school_tm$susceptible_t, type="l", col="blue", ylab="Susceptible over time", xlab = "time index")
