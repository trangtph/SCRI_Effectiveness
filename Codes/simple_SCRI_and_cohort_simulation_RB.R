library(reshape2)

### Function for creating data based on Bernouilli per day approach

bin.days <- function (n, days, weekly.prob, vaccine.effectiveness) {
  prob <- weekly.prob
  ve <- vaccine.effectiveness
  id <- c(1:n)
  exp <- c(rep(1,n/2), rep(0,n/2))
  simul.bin <- matrix(nrow = n, ncol = days + 2)
  simul.bin[,1] <- id
  simul.bin[,2] <- exp
  
  ## Create data for first week (including 3-7 days baseline interval)
  # vector of events
  events <- rbinom(n,1, prob*7/7)
  # vector of random day(column numbers)
  rand.col <- 2 + sample(1:7,n, replace = TRUE)
  # references day (& therefore in which column), the random event should appear
  id.cols <- rep(0,n)
  id.cols[which(events==1)] <- rand.col[which(events==1)]
  # create matrix column names
  colnames <- c("ID", "exposure")
  # populate each column representing day with 0 or 1 events based on column random ref number
  for (col in 3:9) {
    simul.bin[,col] <- ifelse(id.cols==col, 1, 0)
    colnames <- c(colnames, paste0("day", col-2))
  }
  
  # Create data for days 8 to 16 
  exp.prob <- prob*(1-ve/2)
  events <- c(rbinom(n/2,1, exp.prob*(16-7)/7), rbinom(n/2,1,prob*(16-7)/7))
  rand.col <- 2 + sample(10:18,n, replace = TRUE)
  id.cols <- rep(0,n)
  id.cols[which(events==1)] <- rand.col[which(events==1)]
  for (col in 10:18) {
    simul.bin[,col] <- ifelse(id.cols==col  & rowSums(simul.bin[,3:9])==0, 1, 0)
    colnames <- c(colnames, paste0("day", col-2))
  }
  
  # Create data for days 17 to 35 
  exp.prob <- prob*(1-ve)
  events <- c(rbinom(n/2,1, exp.prob*(days-16)/7), rbinom(n/2,1,prob*(days-16)/7))
  rand.col <- 2 + sample(17:days, n, replace = TRUE)
  id.cols <- rep(0,n)
  id.cols[which(events==1)] <- rand.col[which(events==1)]
  for (col in 19:(days+2)) {
    simul.bin[,col] <- ifelse(id.cols==col & rowSums(simul.bin[,3:18])==0, 1, 0)
    colnames <- c(colnames, paste0("day", col-2))
  }
  
  colnames(simul.bin) <- colnames
  return(simul.bin)
  
  # Create data for days 36 to 150 
  exp.prob <- prob*(1-ve)
  events <- c(rbinom(n/2,1, exp.prob*(days-35)/7), rbinom(n/2,1,prob*(days-35)/7))
  rand.col <- 2 + sample(36:days, n, replace = TRUE)
  id.cols <- rep(0,n)
  id.cols[which(events==1)] <- rand.col[which(events==1)]
  for (col in 38:(days+2)) {
    simul.bin[,col] <- ifelse(id.cols==col & rowSums(simul.bin[,18:36])==0, 1, 0)
    colnames <- c(colnames, paste0("day", col-2))
  }
  
  colnames(simul.bin) <- colnames
  return(simul.bin)
  
}

# Create group of 2000 individuals (1:1 exp-to-ctrls ratio) with 150 days of observations,
# with a weekly probability of 0.05 of the outcome & vaccine effectiveness of 66%
s1  <- bin.days(4000, 150, 0.05, 0.66)

##########################################################
########### Reshape and analyse
##  Create a SCRI dataset subsetted on events and exposures
scri.s1 <- as.data.frame(s1[which(s1[,2]==1 & rowSums(s1[,c(5:9, 19:37)])>0), c(1:37)])

scri.s1$base <- ifelse(rowSums(scri.s1[5:9])>0, 1, 0)
scri.s1$risk <- ifelse(rowSums(scri.s1[19:37])>0, 1, 0)
scri.s1 <- scri.s1[,c(1:2, 38:39)]

scri.s1.long <- melt(scri.s1, id.vars = c("ID", "exposure"))
names(scri.s1.long)[3] <- "Period"
names(scri.s1.long)[4] <- "Event"
scri.s1.long$loginterval <- ifelse(scri.s1.long$Period=="risk", log(35-16), log(7-2))

library(survival)
summary(clogit(Event ~ Period + strata(ID) + offset(loginterval), data = scri.s1.long))