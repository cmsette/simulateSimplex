####################Defines function to analyze simplex (handles 2- through 5-players, discrete & continuous). Default 0 is 1e-14.
testSimplex <- function(W, gen_time, zero=1e-14){
	require(foreach); require(deSolve)
	nCores = parallel::detectCores(); doParallel::registerDoParallel(nCores)
	options(warn=-1)
	if(length(W) == 4){
		W_out <- matrix(W, ncol=2, byrow=T); rownames(W_out)<-colnames(W_out)<-c("R", "P")
		eq <- test2(W_out)
		if(eq[3]==1){prop2 <- matrix(c(eq[1]+0.01, 1-(eq[1]+0.01), eq[1]-0.01, 1-(eq[1]-0.01)), ncol=2, byrow=T)}else 									{prop2 <- matrix(c(1/3, 1-1/3, 2/3, 1-2/3), ncol=2, byrow=T)}
		
		simulations <- foreach::foreach(i = length(prop2[,1]), .combine="rbind") %dopar% {  
			sim_out <- step_2Cont(W=W, Time=c(0.1, 0.2), State=list(r=prop2[i,1], p=prop2[i,2]), Pars=NA)[1]; sim_out
		}		
		
		mapply(step_2Cont, W=W, Time=c(0.1, 0.2), State=list(r=prop2[1], p=prop2[2]), Pars=NA)

		return(prop2)
	}
	else if(length(W) == 9){
		W_out<-matrix(W, ncol=3, byrow=T); rownames(W_out)<-colnames(W_out)<-c("R", "P", "S")
		prop3 <- matrix(c(0.01, 0.01, 0.98, 0.01, 0.98, 0.01, 0.98, 0.01, 0.01, 0.30, 0.30, 0.40, 0.30, 0.40, 0.30, 0.40, 0.30, 0.30), ncol=3, 		byrow=T)
		#if(test3(W_out, zero)[4]==1){}


		return(prop3)
	}
	else if(length(W) == 16){
		W_out<-matrix(W, ncol=4, byrow=T); rownames(W_out)<-colnames(W_out)<-c("R", "P", "S", "L")
		prop4 <- matrix(c(0.01, 0.01, 0.01, 0.97, 0.01, 0.01, 0.97, 0.01, 0.01, 0.97, 0.01, 0.01, 0.97, 0.01, 0.01, 0.01, 0.22, 0.22, 0.22, 		0.34, 0.22, 0.22, 0.34, 0.22, 0.22, 0.34, 0.22, 0.22, 0.34, 0.22, 0.22, 0.22), ncol=4, byrow=T)
		#if(test4(W_out, zero)[5]==1){}


		return(prop4)
	}
	else if(length(W) == 25){
		W_out<-matrix(as.numeric(W), ncol=5, byrow=T); rownames(W_out)<-colnames(W_out)<-c("R", "P", "S", "L", "K")
		prop5 <- matrix(c(0.01, 0.01, 0.01, 0.01, 0.96, 0.01, 0.01, 0.01, 0.96, 0.01, 0.01, 0.01, 0.96, 0.01, 0.01, 0.01, 0.96, 0.01, 0.01, 		0.01, 0.96, 0.01, 0.01, 0.01, 0.01, 0.18, 0.18, 0.18, 0.18, 0.28, 0.18, 0.18, 0.18, 0.28, 0.18, 0.18, 0.18, 0.28, 0.18, 0.18, 0.18, 		0.28, 0.18, 0.18, 0.18, 0.28, 0.18, 0.18, 0.18, 0.18), ncol=5, byrow=T)
		#if(test5(W_out, zero)[6]==1){}
		
		
			times <- c(0,0.1,0.2); stop <- F
			yini <- c(r = prop5[i,1], p = prop5[i,2], s = prop5[i,3], l = prop5[i,4], ck = prop5[i,5])
			
		if(gen_time = "continuous"){}else
		if(gen_time = "discrete"){}
		else{print("Invalid generation time")}
	
	
	
		
		

		return(prop5)
	}
	else{print("Error in simplex dimension")}
	options(warn=0)
}


# Solve for eq, then just rotate around it, 360˚
# 1 time step, calculate the direction (inward vs. outward)


### testSimplex
# 1) read in W
# 2) step from values list
# 3) define stop points -> delta is really small or any value approaches 0
# 4) returns final value for eq (or zero)

####################Simulates flow from starting shares
##########Runs deSolve ode solver for incresing time steps, as long as the shares difference between time steps > zero
steps <- function(FUN, yini, times){
	while(stop == F){
		ode(func = FUN, W = W, y = yini, parms = NA, times = times)
		if(all(out[3,-1]-out[2,-1] <= zero)){stop <- T} else{times[2:3] <- times[2:3] + 0.2}
	}		
	ifelse(any(out[3,-1] <= zero), return(out), return(out))
}
	


####################Functions to calculate equilibrium values (2- through 5-players)
##########Calculates delta matrix, finds zero eigenvector (equilibrium), scales eigenvector to unity (simplex solution)		
###Equilibrium function for 2-side games
test2<-function(W){ 
	eqn<-c(W[1,1]-W[2,1], W[1,2]-W[2,2])												#creates delta vector
	soln<-c(-eqn[2]/(eqn[1]-eqn[2]), 1-(-eqn[2]/(eqn[1]-eqn[2])))					#algebraically solves for eq
	if(is.na(all(soln))){return(c(soln,NA))} else{
		if(all(as.matrix(soln) >= 0)){return(c(soln,1))} 							#returns eq solution	
		else if(-(eqn[2]-eqn[1]) < 1){soln<-c(0,1); return(c(soln,0))} 				#dominant P has negative slope
		else if(-(eqn[2]-eqn[1]) > 1){soln<-c(1,0); return(c(soln,0))}				#dominant R has positive slope
		else if(eqn[1] > 0 && eqn[2] > 0){return(c(1,0,0))}else{return(c(0,1,0))} 	# dominant has > eqn values
	}
}
###Equilibrium function for 3-side games
test3<-function(W, zero){ 
	eqn<-matrix(c(W[1,1]-W[2,1], W[2,1]-W[3,1], W[3,1]-W[1,1],  W[1,2]-W[2,2], W[2,2]-W[3,2], W[3,2]-W[1,2], 								W[1,3]-W[2,3], W[2,3]-W[3,3], W[3,3]-W[1,3]), ncol=3, byrow=F)								#creates delta matrix
	soln<-eigen(eqn)$vector; soln<-soln[,3]/sum(soln[,3])						#calculates & scales eigenvector for null set, eq
	test1<-eigen(eqn)$value[3]													#saves corresponding eigenvalue, should = 0
	test2<-abs(eqn %*% soln)													#multiplies delta matrix by solution, should = 0
	if(abs(Re(test1)) > zero | any(abs(test2)  > zero)){return(c(rep(NA,3),NA))} else 				#validation
	if(all(as.matrix(Re(soln)) >= 0)){return(c(Re(soln),1))} else{return(c(rep(NA,3),0))}			#returns eq soln
}
###Equilibrium function for 4-side games
test4<-function(W, zero){ 
	eqn<-matrix(c(W[1,1]-W[2,1], W[2,1]-W[3,1], W[3,1]-W[4,1], W[4,1]-W[1,1], 																W[1,2]-W[2,2], W[2,2]-W[3,2], W[3,2]-W[4,2], W[4,2]-W[1,2],																			W[1,3]-W[2,3], W[2,3]-W[3,3], W[3,3]-W[4,3], W[4,3]-W[1,3],																			W[1,4]-W[2,4], W[2,4]-W[3,4], W[3,4]-W[4,4], W[4,4]-W[1,4]), ncol=4, byrow=F)				#creates delta matrix
	soln<-eigen(eqn)$vector; soln<-soln[,4]/sum(soln[,4])						#calculates & scales eigenvector for null set,
	test1<-eigen(eqn)$value[4]													#saves corresponding eigenvalue, should = 0
	test2<-abs(eqn %*% soln)													#multiplies delta matrix by solution, should = 0
	if(abs(Re(test1)) > zero | any(abs(test2) > zero)){return(c(rep(NA,4),NA))} else 				#validation 
	if(all(as.matrix(Re(soln)) >= 0)){return(c(Re(soln),1))} else{return(c(rep(NA,4),0))}			#returns eq soln
}
###Equilibrium function for 5-side games
test5<-function(W, zero){ 
	eqn<-matrix(c(W[1,1]-W[2,1], W[2,1]-W[3,1], W[3,1]-W[4,1], W[4,1]-W[5,1], W[5,1]-W[1,1],												W[1,2]-W[2,2], W[2,2]-W[3,2], W[3,2]-W[4,2], W[4,2]-W[5,2], W[5,2]-W[1,2],															W[1,3]-W[2,3], W[2,3]-W[3,3], W[3,3]-W[4,3], W[4,3]-W[5,3], W[5,3]-W[1,3],															W[1,4]-W[2,4], W[2,4]-W[3,4], W[3,4]-W[4,4], W[4,4]-W[5,4], W[5,4]-W[1,4],															W[1,5]-W[2,5], W[2,5]-W[3,5], W[3,5]-W[4,5], W[4,5]-W[5,5], W[5,5]-W[1,5]), ncol=5, byrow=F)		#creates delta
	soln<-eigen(eqn)$vector; soln<-soln[,5]/sum(soln[,5])						#calculates & scales eigenvector for null set
	test1<-eigen(eqn)$value[5]													#saves corresponding eigenvalue, should = 0
	test2<-abs(eqn %*% soln)													#multiplies delta matrix by solution, should = 0
	if(abs(Re(test1)) > zero | any(abs(test2) > zero)){return(c(rep(NA,5),NA))} else 			#validation
	if(all(as.matrix(Re(soln)) >= 0)){return(c(Re(soln),1))} else{return(c(rep(NA,5),0))}		#returns eq soln
}

#################### Functions to calculate delta shares (2- through 5-players) based on payoff matrix
########## Calculates deltas for all continuous time steps passed to function		
### 2-side continuous games
step_2Cont<-function(W, Time, State, Pars){
	with(as.list(c(State,Pars)),{
		w_r <- (W[1]*r) + (W[2]*p)
		w_p <- (W[3]*r) + (W[4]*p)
		w_bar <- r*w_r + p*w_p
		dr <- r * (w_r - w_bar)
		dp <- p * (w_p - w_bar)
    return(list(c(dr, dp)))}
)}
### 3-side continuous games
step_3Cont<-function(W, Time, State, Pars){
	with(as.list(c(State,Pars)),{
		w_r <- (W[1]*r) + (W[2]*p) + (W[3]*s)
		w_p <- (W[4]*r) + (W[5]*p) + (W[6]*s)
		w_s <- (W[7]*r) + (W[8]*p) + (W[9]*s)
		w_bar <- r*w_r + p*w_p + s*w_s
		dr <- r * (w_r - w_bar)
		dp <- p * (w_p - w_bar)
		ds <- s * (w_s - w_bar)
	return(list(c(dr, dp, ds)))}
)}
### 4-side continuous games
step_4Cont<-function(W, Time, State, Pars){
	with(as.list(c(State,Pars)),{
		w_r <- (W[1]*r) + (W[2]*p) + (W[3]*s) + (W[4]*d)
		w_p <- (W[5]*r) + (W[6]*p) + (W[7]*s) + (W[8]*d)
		w_s <- (W[9]*r) + (W[10]*p) + (W[11]*s) + (W[12]*d)
		w_d <- (W[13]*r) + (W[14]*p) + (W[15]*s) + (W[16]*d)
		w_bar <- r*w_r + p*w_p + s*w_s + d*w_d
		dr <- r * (w_r - w_bar)
		dp <- p * (w_p - w_bar)
		ds <- s * (w_s - w_bar)
		dd <- d * (w_d - w_bar)
	return(list(c(dr, dp, ds, dd)))}
)}
### 5-side continuous games
step_5Cont<-function(W, Time, State, Pars){
	with(as.list(c(State,Pars)),{
		w_r <- (W[1]*r) + (W[2]*p) + (W[3]*s) + (W[4]*l) + (W[5]*ck)
		w_p <- (W[6]*r) + (W[7]*p) + (W[8]*s) + (W[9]*l) + (W[10]*ck)
		w_s <- (W[11]*r) + (W[12]*p) + (W[13]*s) + (W[14]*l) + (W[15]*ck)
		w_l <- (W[16]*r) + (W[17]*p) + (W[18]*s) + (W[19]*l) + (W[20]*ck)
		w_ck <- (W[21]*r) + (W[22]*p) + (W[23]*s) + (W[24]*l) + (W[25]*ck)
		w_bar <- r*w_r + p*w_p + s*w_s + l*w_l + ck*w_ck
		dr <- r * (w_r - w_bar)
		dp <- p * (w_p - w_bar)
		ds <- s * (w_s - w_bar)
		dl <- l * (w_l - w_bar)
		dck <- ck * (w_ck - w_bar)
	return(list(c(dr, dp, ds, dl, dck)))}
)}
########## Calculates deltas for all discrete time steps passed to function	
### 2-side discrete games
step_2Disc<-function(W, Time, State, Pars){
	with(as.list(c(State,Pars)),{
		w_r <- (W[1]*r) + (W[2]*p)
		w_p <- (W[3]*r) + (W[4]*p)
		w_bar <- r*w_r + p*w_p
		dr <- (r * (w_r / w_bar)) - r
		dp <- (p * (w_p / w_bar)) - p
	return(list(c(dr, dp)))}
)}
### 3-side discrete games
step_3Disc<-function(W, Time, State, Pars){
	with(as.list(c(State,Pars)),{
		w_r <- (W[1]*r) + (W[2]*p) + (W[3]*s)
		w_p <- (W[4]*r) + (W[5]*p) + (W[6]*s)
		w_s <- (W[7]*r) + (W[8]*p) + (W[9]*s)
		w_bar <- r*w_r + p*w_p + s*w_s
		dr <- (r * (w_r / w_bar)) - r
		dp <- (p * (w_p / w_bar)) - p
		ds <- (s * (w_s / w_bar)) - s
	return(list(c(dr, dp, ds)))}
)}
### 4-side discrete games
step_4Disc<-function(W, Time, State, Pars){
	with(as.list(c(State,Pars)),{
		w_r <- (W[1]*r) + (W[2]*p) + (W[3]*s) + (W[4]*d)
		w_p <- (W[5]*r) + (W[6]*p) + (W[7]*s) + (W[8]*d)
		w_s <- (W[9]*r) + (W[10]*p) + (W[11]*s) + (W[12]*d)
		w_d <- (W[13]*r) + (W[14]*p) + (W[15]*s) + (W[16]*d)
		w_bar <- r*w_r + p*w_p + s*w_s + d*w_d
		dr <- (r * (w_r / w_bar)) - r
		dp <- (p * (w_p / w_bar)) - p
		ds <- (s * (w_s / w_bar)) - s
		dd <- (d * (w_d / w_bar)) - d
	return(list(c(dr, dp, ds, dd)))}
)}
### 5-side discrete games
step_5Disc<-function(W, Time, State, Pars){
	with(as.list(c(State,Pars)),{
		w_r <- (W[1]*r) + (W[2]*p) + (W[3]*s) + (W[4]*l) + (W[5]*ck)
		w_p <- (W[6]*r) + (W[7]*p) + (W[8]*s) + (W[9]*l) + (W[10]*ck)
		w_s <- (W[11]*r) + (W[12]*p) + (W[13]*s) + (W[14]*l) + (W[15]*ck)
		w_l <- (W[16]*r) + (W[17]*p) + (W[18]*s) + (W[19]*l) + (W[20]*ck)
		w_ck <- (W[21]*r) + (W[22]*p) + (W[23]*s) + (W[24]*l) + (W[25]*ck)
		w_bar <- r*w_r + p*w_p + s*w_s + l*w_l + ck*w_ck
		dr <- (r * (w_r / w_bar)) - r
		dp <- (p * (w_p / w_bar)) - p
		ds <- (s * (w_s / w_bar)) - s
		dl <- (l * (w_l / w_bar)) - l
		dck <- (ck * (w_ck / w_bar)) - ck
	return(list(c(dr, dp, ds, dl, dck)))}
)}




