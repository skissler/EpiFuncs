library(tidyverse) 

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

# For evaluating the total kernel CDF: 
cumckernel <- function(t, tinf, parslist, ckernel){

	sum(unlist(pmap(
		list(t-tinf[tinf<Inf], parslist[tinf<Inf]),
		ckernel
		)))

}

# Main simulation algorithm
infdurabm <- function(pkernel, ckernel, parslist, invckernel=NA, initinf=1, initstep=1, maxits=1000){

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
		cumrate <- sum(tinf==Inf)*cumforce/N
		if(cumrate < 1e-4){break()}
		nevents <- rpois(1,cumrate)

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
		tstep <- tprop - t
		t <- tprop

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

beta <- 2
gamma <- 1/5
parslist <- as.list(rexp(20,gamma)) %>% 
	map(~ c(tstar=., beta=beta))




# simout <- infdurabm(pkernel_sir, ckernel_sir, parslist, invckernel=NA, initinf=1, initstep=1, maxits=1000)
