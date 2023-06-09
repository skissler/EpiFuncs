library(tidyverse) 
library(purrr)
library(deSolve)

# Bisection algorithm for a monotonically increasing function
bisect <- function(target, lwr, upr, func, tol=1e-6){

	x <- lwr + (upr-lwr)/2
	jump <- (upr-lwr)/4
	testval <- func(x)

	while(abs(testval-target)>tol){

		if(testval<target){
			x <- x + jump
		} else {
			x <- x - jump
		}

		jump <- jump/2
		testval <- func(x)

	}

	return(x)

}

# Time series of cumulative new infections: 
plotcuminf <- function(simout){

	N <- nrow(simout) 

	simoutshort <- filter(simout, tinf<Inf)
	tvals <- seq(from=0, to=max(simoutshort$tinf)+1, by=0.01)
	tinf <- simoutshort$tinf
	yvals <- unlist(lapply(tvals, function(x){sum(tinf<x)}))

	out <- tibble(t=tvals, y=yvals) %>% 
		ggplot(aes(x=t, y=y)) + 
			geom_line(size=1) + 
			geom_hline(yintercept=N, lty="dashed", col="grey") + 
			theme_classic() + 
			labs(x="Time", y="Cumulative infections")

	return(out)

}

# Time series of the current force of infection: 
plotfoi <- function(simout, pkernel, parslist){

	simoutshort <- filter(simout, tinf<Inf)
	tvals <- seq(from=0, to=max(simoutshort$tinf)+1, by=0.01)
	tinf <- simoutshort$tinf
	parslistshort <- parslist[simoutshort$id]

	yvals <- as.list(tvals) %>% 
		map(~ sum(unlist(pmap(list(.-tinf, parslistshort), pkernel)))) %>% 
		unlist()

	out <- tibble(t=tvals, y=yvals) %>% 
		ggplot(aes(x=t, y=y)) + 
			geom_line(size=1) + 
			theme_classic() + 
			labs(x="Time", y="Force of infection")

	return(out)

}

# Plot a WAIFW tree: 


# Plot individual infection kernels: 
plotai <- function(pkernel,parslist,tmax,tmin=0,yjittersd=0.001){

	tau <- seq(from=tmin, to=tmax, by=0.01)
	out <- as.list(tau) %>% 
		map(~ unlist(pmap(
			list(rep(., length(parslist)), parslist),
			pkernel
			))) %>% 
		map(~ tibble(id=1:length(.), y=.)) %>% 
		imap(~ mutate(., tindex=.y)) %>% 
		bind_rows() %>% 
		left_join(tibble(tau=tau, tindex=1:length(tau)), by="tindex") %>% 
		group_by(id) %>% 
		mutate(yjitter=rnorm(1,0,yjittersd)) %>% 
		ggplot(aes(x=tau, y=y+yjitter, group=factor(id))) + 
			geom_line(alpha=0.2) + 
			theme_classic() + 
			labs(x="Time since infection", y="Infectiousness")

	return(out)

}


# Plot the summed infection kernel across all individuals: 
plotA <- function(pkernel,parslist,tmax,tmin=0,compfun=NA){

	N <- length(parslist)
	tau <- seq(from=tmin, to=tmax, by=0.01)
	yvals <- as.list(tau) %>% 
		map(~ sum(unlist(pmap(
			list(rep(., length(parslist)), parslist),
			pkernel
			)))) %>% 
		unlist()
	out <- tibble(tau=tau, y=yvals/N) %>% 
		ggplot(aes(x=tau, y=y)) + 
			geom_line()  + 
			theme_classic() + 
			labs(x="Time since infection", y="Mean infectiousness")

	if(is.function(compfun)){
		ycomp <- unlist(lapply(tau, compfun))
		out <- out + geom_line(
			data=tibble(x=tau, y=ycomp), 
			aes(x=x, y=y), size=1.5, col="dodgerblue", alpha=0.4)
	}

	return(out)

}

# For evaluating the total kernel CDF: 
cumckernel <- function(t, tinf, parslist, ckernel){

	sum(unlist(pmap(
		list(t-tinf[tinf<Inf], parslist[tinf<Inf]),
		ckernel
		)))

}

# Main simulation algorithm
infdurabm <- function(pkernel, ckernel, parslist, invckernel=NA, initinf=1, initstep=1, maxits=1000, quiet=FALSE){

	# Initialize output vectors: 
	N <- length(parslist)
	tinf <- rep(Inf, N)
	tinf[initinf] <- 0
	whoinf <- rep(-1,N)
	whoinf[initinf] <- 0

	# Initialize looping variables: 
	t <- 0
	tstep <- initstep

	# Simulate infections:
	for(iteration in 1:maxits){

		# Propose the next time: 
		tprop <- t+tstep

		# Calculate the total force of infection over this time span: 
		cumforcestart <- cumckernel(t, tinf, parslist, ckernel)
		cumforceend <- cumckernel(tprop, tinf, parslist, ckernel)
		cumforce <- cumforceend - cumforcestart
			
		# Calculate the number of "events" that would occur in the current system given this rate: 

		pinf <- 1-exp(-cumforce/N)
		if(sum(tinf==Inf)*pinf < 1e-6){break()}
		nevents <- rbinom(1, sum(tinf==Inf), pinf)

		# --------
		# cumrate <- sum(tinf==Inf)*cumforce/N
		# if(cumrate < 1e-6){break()}
		# nevents <- rpois(1,cumrate)
		# --------

		# Get the timing of those events: 
		if(nevents==0){
			t <- tprop
			next()
		}

		# Get uniform draw for the inverse CDF
		draw <- min(runif(nevents))

		# Calculate "event times" using inverse CDF
		func <- function(x){
			(cumckernel(x, tinf, parslist, ckernel) - cumforcestart)/cumforce
		}

		# Update tprop to be the time of the new infection (the first)
		tprop <- bisect(draw, lwr=t, upr=tprop, func=func)		

		# Randomly draw a new person to infect
		newinf <- sample((1:N)[tinf==Inf], 1)

		# Who did the infecting? 
		pvals <- unlist(pmap(
			list(tprop-tinf[tinf<Inf], parslist[tinf<Inf]),
			pkernel
			))
		newwhoinf <- ((1:N)[tinf<Inf])[which(as.vector(rmultinom(1,1,pvals))==1)]

		# Update the output vectors: 
		tinf[newinf] <- tprop
		whoinf[newinf] <- newwhoinf

		# Update variables for new loop: 
		t <- tprop

	}

	if(!quiet){

		if(iteration >= maxits){
			print(paste0("Reached maxits at time ",round(t,2)))
		} else {
			print(paste0("Epidemic ended at time ",round(t,2)," after ",iteration," iterations"))
		}

	}

	out <- tibble(id=1:N, tinf=tinf, whoinf=whoinf)

	return(out)

}

# =============================================================================
# Example
# =============================================================================

pkernel_sir <- function(tau,pars){
	with(as.list(pars), {
		if(tau>0 & tau<=tstar){
			return(beta)
		} else {
			return(0)
		}
	})
}

ckernel_sir <- function(tau,pars){
	with(as.list(pars),{
		if(tau<=0){
			return(0)
		} else if(tau>0 & tau<=tstar){
			return(beta*tau)
		} else {
			return(beta*tstar)
		}
		})
}

beta <- 1/2
gamma <- 1/5
parslist <- as.list(rexp(200,gamma)) %>% 
	map(~ c(tstar=., beta=beta))

simout <- infdurabm(pkernel_sir, ckernel_sir, parslist, invckernel=NA, initinf=1, initstep=1, maxits=1000)

fig_cuminf <- plotcuminf(simout)
fig_foi <- plotfoi(simout, pkernel_sir, parslist)

fig_ai <- plotai(pkernel_sir, parslist, tmin=0.01, tmax=30)

cf <- function(x){beta*exp(-gamma*x)}
fig_A <- plotA(pkernel_sir, parslist, tmax=30, compfun=cf)

# =============================================================================
# Repeated simulations 
# =============================================================================

beta <- 1/2
gamma <- 1/5

simoutlist <- list() 
print(paste0("Loop started: ",Sys.time()))
for(simnum in 1:200){
	parslist <- as.list(rexp(200,gamma)) %>% 
		map(~ c(tstar=., beta=beta))
	simout <- infdurabm(pkernel_sir, ckernel_sir, parslist, invckernel=NA, initinf=1, initstep=1, maxits=1000, quiet=TRUE)
	tstardf <- tibble(id=1:length(parslist), tstar=unlist(map(parslist,~.["tstar"])))
	simout <- simout %>% 
		left_join(tstardf, by="id") %>% 
		mutate(trec=tinf+tstar)
	simoutlist[[simnum]] <- simout
	if(simnum%%10==0){
		print(paste0("Simulation number ",simnum," completed"))
	}
}
print(paste0("Loop ended: ",Sys.time()))

# save(simoutlist, file="~/Desktop/simoutlist.RData")

simoutdf <- simoutlist %>% 
	imap(~ mutate(.x, sim=.y)) %>% 
	bind_rows() %>% 
	filter(tinf<Inf) %>% 
	mutate(counter=1) %>% 
	group_by(sim) %>% 
	arrange(sim, trec) %>% 
	group_by(sim) %>% 
	mutate(cuminf=cumsum(counter)) %>% 
	select(-counter) %>% 
	split(.$sim) %>% 
	map(~ bind_rows(., tibble(id=-1, tinf=Inf, whoinf=0, tstar=0, trec=Inf, sim=min(.$sim), cuminf=last(.$cuminf)))) %>% 
	bind_rows()

sir <- function(t,state,parameters){
	with(as.list(c(state,parameters)), {

		dS <- -beta*I*S
		dI <- beta*I*S - gamma*I
		dR <- gamma*I

		list(c(dS, dI, dR))

		})
}

N <- length(parslist)
parameters <- c(beta=beta, gamma=gamma) 
state <- c(S=1-1/N, I=1/N, R=0)
times <- seq(from=0, to=100, by=0.01)

sirout <- ode(y = state, times = times, func = sir, parms = parameters) %>% 
	data.frame() %>% 
	as_tibble()

figsimcurves <- simoutdf %>% 
	ggplot() + 
		geom_step(aes(x=trec, y=cuminf, group=sim), col="grey", alpha=0.2) + 
		geom_line(data=sirout, aes(x=time, y=R*N), col="blue", size=1) + 
		theme_classic() 

# Compare with a gillespie algorithm: 
source('code/gillespie.R')

states <- c(S=199, I=1, R=0)
parms <- c(beta=beta, gamma=gamma, N=sum(states))
rates <- c(
	"S -> I"="beta*S*I/N", 
	"I -> R"="gamma*I")

gilloutlist <- list()
for(simnum in 1:200){
	gillout <- gillespie(states,parms,rates,maxits=10000)
	gilloutlist[[simnum]] <- gillout 
	if(simnum%%10==0){
		print(paste0("Simulation number ",simnum," completed"))
	}	
}

gilloutdf <- gilloutlist %>% 
	imap(~ mutate(.x, sim=.y)) %>% 
	map(~ bind_rows(., tibble(t=Inf, sim=min(.$sim), S=last(.$S), I=last(.$I), R=last(.$R)))) %>% 
	bind_rows()

figgillcurves <- gilloutdf %>% 
	ggplot() + 
		geom_step(aes(x=t, y=R, group=sim), col="grey", alpha=0.2) + 
		geom_line(data=sirout, aes(x=time, y=R*N), col="blue", size=1) + 
		theme_classic() 


