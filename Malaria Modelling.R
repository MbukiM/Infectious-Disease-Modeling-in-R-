library(deSolve)
library(ggplot2)

lambda=0.05 # Transmission rate 
r=1/200 # recovery rate in human
N=100 # Total Population 


theta_init =c("lambda"=lambda,"r"=r, "N"=N)


##SIR model Function 
# has three parameter time state parameter 


# First function 
SIS_deterministic<-function(t,x,parms){
  # version from Smith and McKenzie, 2004
  
  ## state variables 
  # the number of infected 
  # X is a vector I 
  Hi = x[1]
  
  ## parameter values
  lambda = parms["lambda"]
  r = parms["r"]
  N= parms["N"]
  
  ##variations
  dHi=lambda*(N-Hi)*Hi/N - r*Hi # Model defination 
  res  = c(dHi) 
  list(res)
}


# Second Function to simulate 
simulate_SIS=function(parameters, maxtime = 500){
  
  #parameters
  # Renaming 
  #parms = c(parameters["lambda"],parameters["r"], parameters["N"])
  
  #initial conditions
  init <- c(0.005*parameters["N"])  
  
  #simulation of the Ordinary Differential Equation
  # Initial variable time step for simulation the model 
  # Unit of time is in days 
  
  temps <- seq(0, maxtime) 
  solveSIS <- lsoda(y =init, times=temps, func = SIS_deterministic, 
                    parms = parameters)
  solutionSIS=as.data.frame(solveSIS)
  names(solutionSIS)=c("time","Hi")
  
  return(solutionSIS)
}



simul=simulate_SIS(parameters = theta_init)

ggplot(simul)+
  geom_line(aes(x=time, y=Hi))+ 
  labs(y="Proportion of\ninfected individuals", x="Time (days)")

##########################################################################################
# Stochastic Model 

simulate_SIS_stochastic=function(lambda, r, N, init, steps=1000){
  
  #initial conditions
  # 
  solution=data.frame(I=rep(NA,steps ), time=rep(NA,steps))
  solution[1,]=c(0.005*N,1)
  
  U1=runif(steps)
  U2=runif(steps)
  
# For loop 
  
  for(t in 2:steps){
    current_i=solution$I[t-1]
    rate_infect=lambda*current_i*(N-current_i)/N
    rate_recovery=r*current_i
    alpha0=rate_infect+rate_recovery
    if(alpha0>0){
      solution$time[t]=solution$time[t-1]+(1/alpha0)*log(1/U1[t])
      if(rate_infect/alpha0>U2[t]){
        solution$I[t]=solution$I[t-1]+1
      } else {
        solution$I[t]=solution$I[t-1]-1
        
      }
    } else {
      solution$I[t]=0
      solution$time[t]=solution$time[t-1]+1
    }
    
  }
  
  return(solution)
}


simul_sto=simulate_SIS_stochastic(lambda, r, N=N, init=50, steps=1000)
simul_sto2=simulate_SIS_stochastic(lambda, r, N=N, init=50, steps=1000)
simul_sto3=simulate_SIS_stochastic(lambda, r, N=N, init=50, steps=1000)


ggplot(simul)+
  geom_line(aes(x=time, y=Hi))+ 
  geom_line(data=simul_sto,aes(x=time, y=I, col="sto"))+ 
  geom_line(data=simul_sto2,aes(x=time, y=I, col="sto2"))+ 
  geom_line(data=simul_sto3,aes(x=time, y=I, col="sto3"))+ 
  labs(y="Proportion of\ninfected individuals", x="Time (days)")





library(deSolve)
library(ggplot2)

theta_before=c("m"=5,"a"=0.5,"b"=0.09, "c"=1,"p"=0.95,"v"=12, "r"=1/80)
theta_control=c("m"=0.05,"a"=0.5,"b"=0.09, "c"=1,"p"=0.95,"v"=12, "r"=1/80)


##SIR model
RossMacDonald<-function(t,x,parms){
  # version from Smith and McKenzie, 2004
  
  ## state variables
  Hi = x[1]
  Vi = x[2]
  
  ## parameter values
  m = parms["m"]
  a = parms["a"]
  b = parms["b"]
  c = parms["c"]
  p = parms["p"]
  v = parms["v"]
  r = parms["r"]
  g= -log(p)
  
  ##variations
  dHi=m*a*b*(1-Hi)*Vi - r*Hi
  dVi= a*c*Hi*(exp(-g*v)-Vi) - g*Vi
  res  = c(dHi, dVi)
  list(res)
}

simulate_RM=function(parameters, init,max_time){
  
  #parameters
  parms = c(parameters["m"],parameters["a"],
            parameters["b"],parameters["c"],
            parameters["p"],parameters["v"],
            parameters["r"])
  
  #simulation
  temps <- seq(0,max_time) 
  solveRM <- lsoda(y =init, times=temps, func = RossMacDonald, 
                   parms = parms)
  solutionRM=as.data.frame(solveRM)
  names(solutionRM)=c("time","Hi","Vi")
  
  return(solutionRM)
}

calculate_R0=function(parameters){
  return(c("R0"=as.numeric(parameters["m"]*(parameters["a"]^2)*parameters["b"]*parameters["c"]*(parameters["p"]^parameters["v"])/parameters["r"]/(-log(parameters["p"])))))
}

simul_before=simulate_RM(parameters = theta_before, init=c(0.5, 0.5),max_time = 2000)
simul_control=simulate_RM(parameters = theta_control, init=as.numeric(simul_before[2000,c(2,3)]), max_time = 5000)
simul_after=simulate_RM(parameters = theta_before, init=as.numeric(simul_control[2000,c(2,3)]), max_time = 5000)
simul_control$time=simul_control$time+2000
simul_after$time=simul_after$time+7000


calculate_R0(theta_before)
calculate_R0(theta_control)

ggplot(rbind(simul_before,simul_control, simul_after))+
  geom_line(aes(x=time, y=Hi))+ ylim(0,1)+xlim(100,NA)+
  labs(y="Proportion of\ninfected individuals", x="Time (days)")






















