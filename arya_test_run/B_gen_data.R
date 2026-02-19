gen_outcome_raw=function(Nat,Nnt,Nc){

    # This function generates a table of potential outcomes and complier status
    # Nat: number of always-takers 
    # Nnt: number of never-takers
    # Nc: number of compliers

    # By saying "raw", we mean that we have not attempted to add any effects to the outcomes yet

    # N: total number of units
    N=Nat+Nnt+Nc

    # Generate potential outcomes
    
    y1=rnorm(N,0,.01) #note what y1 equals here does not matter. They are going to be overwritten later in the function add_effects_to_raw_outcomes
    #y0=rbeta(N,2,1)
    #y0=c(rnorm(Nat,1,1),rnorm(Nnt,0,1),rnorm(Nc,1,1))
    #hist(y0)  
    #sort y0 


    d1=c(rep(1,Nat),rep(0,Nnt),rep(1,Nc))
    d0=c(rep(1,Nat),rep(0,Nnt),rep(0,Nc))
    y1=-(d1+d0)+rnorm(N,0,0.2)
    y0=y1


    # Genereate treatment takeup decision
    #u=runif(N)
    #d1=u <= pnorm(0.5+epsilon)
    #d1=u <= pbeta(0.2+y0,2,1)
    #d0=u<=pbeta(y0,2,1)
    
    #d0[order(y0)[1:3]]=1
    #d0=which.min(y0)
    # indices of units
    indices=1:N

    # Combine into a table
    output = cbind(indices,y1,y0,d1,d0)

    colnames(output)=c('indices',"y1_raw","y0_raw","d1","d0")
    
    return(output) 

}

add_effects_to_raw_outcomes_additive=function(outcome_raw,tau){

    # This function adds treatment effects to the raw outcomes
    # outcome_raw: table of raw outcomes and complier status
    # tau1: treatment effect for treated units
    # tau0: treatment effect for control units

    # Not we are assuming that in the row outcomes the order is always-takers, never-takers, compliers


    outcome_raw[,2]=outcome_raw[,3] + tau

    #calibrate the effects to match exactly 
    colnames(outcome_raw)=c('indices',"y1","y0","d1","d0")

    return(outcome_raw) 

}

add_effects_to_raw_outcomes=function(outcome_raw,tau_at,tau_nt,tau_c,N_at,N_nt,N_c){

    # This function adds treatment effects to the raw outcomes
    # outcome_raw: table of raw outcomes and complier status
    # tau1: treatment effect for treated units
    # tau0: treatment effect for control units

    # Not we are assuming that in the row outcomes the order is always-takers, never-takers, compliers


    # Extract raw outcomes
    if (tau_at!=0){
        eff_at =  tau_at*rnorm(N_at,1,0.1)  #
        #eff_at =  rep(tau_at,N_at)

        
    }else {
       #eff_at = rnorm(N_at,0,0.1)
       eff_at =  rep(tau_at,N_at)
       
    }

    if (tau_nt!=0){
        eff_nt =  tau_nt*rnorm(N_nt,1,0.1)
        #eff_nt = rep(tau_nt,N_nt)

    }else {
        #eff_nt = rnorm(N_nt,0.1)
        eff_nt =rep(tau_nt,N_nt)
        
    }

    if (tau_c!=0){
        eff_c =  tau_c*rnorm(N_c,1,0.1)
        #eff_c = rep(tau_c,N_c)
        #eff_c =  tau_c
    }else {
        eff_c = rnorm(N_c,0,0.1)
        #eff_c = rep(tau_c,N_c)
    }
    
    outcome_raw[,2]=outcome_raw[,3] + c(eff_at,eff_nt,eff_c)
    #browser()
    #calibrate the effects to match exactly 
    #For always-takers
    factor_at =  mean(outcome_raw[1:N_at,2]-outcome_raw[1:N_at,3]) -tau_at 
    outcome_raw[1:N_at,2] = outcome_raw[1:N_at,2] - factor_at 

    #For never-takers
    factor_nt = mean(outcome_raw[(N_at+1):(N_at+N_nt),2]-outcome_raw[(N_at+1):(N_at+N_nt),3]) -tau_nt 
    outcome_raw[(N_at+1):(N_at+N_nt),2] = outcome_raw[(N_at+1):(N_at+N_nt),2] -factor_nt 

    #For compliers
    factor_c = mean(outcome_raw[(N_at+N_nt+1):(N_at+N_nt+N_c),2]-outcome_raw[(N_at+N_nt+1):(N_at+N_nt+N_c),3]) - tau_c 
    outcome_raw[(N_at+N_nt+1):(N_at+N_nt+N_c),2] =  outcome_raw[(N_at+N_nt+1):(N_at+N_nt+N_c),2] - factor_c
    

    

    
    colnames(outcome_raw)=c('indices',"y1","y0","d1","d0")

    return(outcome_raw) 

}

gen_assignment_CR=function(N1,N0){ # nolint

    # Generate assignment to treatment and control groups: completly randomized
    # N1: number of units in treatment group
    # N0: number of units in control group

    # N: total number of units
    N=N1+N0

    treated=sample(1:N,N1,replace=FALSE)

    # indices of units
    indices=1:N

    # Generate assignment
    assignment=rep(0,N)
    assignment[treated]=1

    # Combine into a table
    output_assignment = cbind(indices,assignment)

    colnames(output_assignment)=c('indices',"assignment")

    return(output_assignment)

}


gen_assignment_CR_index=function(N1,N0,i){ # nolint
  
  # Generate assignment to treatment and control groups: completly randomized
  # N1: number of units in treatment group
  # N0: number of units in control group
  
  # N: total number of units
  N=N1+N0
  
  treated=sample(1:N,N1,replace=FALSE)
  
  # indices of units
  indices=1:N
  
  # Generate assignment
  assignment=rep(0,N)
  assignment[treated]=1
  

  
  return(assignment)
  
}
gen_data=function(N_at,N_nt,N_c,tau_at,tau_nt,tau_c){

    #generate an outcome table
    outcome_table=gen_outcome_raw(N_at,N_nt,N_c)

    #generate heterogeneous treatment effects 
    outcome_table=add_effects_to_raw_outcomes(outcome_table,tau_at,tau_nt,tau_c,N_at,N_nt,N_c)

    return(outcome_table)

}

gen_data_additive=function(N_at,N_nt,N_c,tau){

    #generate an outcome table
    outcome_table=gen_outcome_raw(N_at,N_nt,N_c)

    #generate heterogeneous treatment effects 
    outcome_table=add_effects_to_raw_outcomes_additive(outcome_table,tau)

    return(outcome_table)

}

gen_data_onesim = function(outcome_table,N1,N0){


    assignment=gen_assignment_CR(N1,N0)[,2]

    D = outcome_table[,'d1'] * assignment + outcome_table[,'d0']*(1-assignment)

    Y = outcome_table[,'y1'] * D + outcome_table[,'y0'] * (1-D)

    output = cbind(Y,D,assignment)

    colnames(output)=c('Y_observed',"D_observed","assignment")

    return(output)
}
