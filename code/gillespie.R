library(tidyverse)
library(stringr)

gillespie <- function(states,parms,rates,maxits=1000){

	# Define function for pulling out the from/to states from 'rates'
	get_from_to <- function(rates, i){
		unlist(lapply(unlist(str_split(names(rates)[i],"->")),str_trim))	
	}

	# Make sure all compartments in 'rates' are defined in 'states'
	if(!setequal(
		names(states),
		unique(unlist(lapply(1:length(rates), function(x){get_from_to(rates,x)}))))){
		stop("All compartments in rates must have initial conditions in states")
	}

	# Initialize looping index 
	indexA <- 1
	# Initialize output table
	out <- list(c(t=0, states))

	while(indexA < maxits){

		# Evaluate the transition rates
		rates_eval <- with(as.list(c(states,parms)),{
			return(unlist(lapply(rates, function(x){eval(parse(text=x))})))
			})

		# Calculate the total rate
		totrate <- sum(rates_eval)

		# If the rate is basically 0, stop looping
		if(totrate<1e-6){break()}

		# Draw the time of the next event
		nexttime <- rexp(1,totrate)
		# Draw which event occurred
		nexttrans <- sample(length(rates_eval), 1, prob=rates_eval)

		# Extract the state-from and the state-to
		from_to <- get_from_to(rates, nexttrans)

		# Update the state vector
		states[from_to[1]] <- states[from_to[1]]-1
		states[from_to[2]] <- states[from_to[2]]+1

		# Append the new state vector to the output table
		out[[length(out)+1]] <- c(t=nexttime, states)
		# Increment the looping index
		indexA <- indexA + 1

	}

	# Bind and format the output table 
	out <- out %>% 
	 	map(~ as_tibble(t(.))) %>% 
		bind_rows() %>% 
		mutate(t=cumsum(t))

	return(out)
}

# =============================================================================
# Testing
# =============================================================================

states <- c(S=999, I=1, R=0)
parms <- c(beta=2, gamma=1, N=sum(states))
rates <- c(
	"S -> I"="beta*S*I/N", 
	"I -> R"="gamma*I")
temp <- gillespie(states,parms,rates,maxits=10000)

tempfig <- temp %>% 
	ggplot(aes(x=t, y=I)) + 
		geom_line() + 
		theme_classic()




# states <- c(S=998, E=1, I=1, R=0)
# parms <- c(beta=2, nu=1, gamma=1, sigma=1/10, N=sum(states))
# rates <- c(
# 	"S -> E"="beta*S*I/N", 
# 	"E -> I"="nu*E",
# 	"I -> R"="gamma*I",
# 	"R -> S"="sigma*R")
# temp2 <- gillespie(states,parms,rates,maxits=10000)

# temp2fig <- temp2 %>% 
# 	ggplot(aes(x=t, y=I)) + 
# 		geom_line() + 
# 		theme_classic()

