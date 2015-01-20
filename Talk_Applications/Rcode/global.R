#SCRIPT TO A FLOW-LIMITED ORGAN

#remove all current objects in the workspace
 #rm(list=ls(all=TRUE))
 #graphics.off()
 
#Set the working directory
  #master.dir <- "E:/a-jobs/PAGANZ15/Talk_Principles/Rcode"    #PRI Server
  #setwd(master.dir)  

#Load libraries
  library(ggplot2)
  library(doBy)
  library(plyr)
  library(grid)
  library(deSolve)
   

#--------------------------------------------------------------------------------------------------   

#Customize ggplot2 theme - R 2.15.3+
 theme_bw2 <- theme_set(theme_bw(base_size = 22))  
 theme_bw2 <- theme_update(plot.margin = unit(c(0.1,0.1,0.1,0.1), "npc"),
 axis.title.x=element_text(size = 18, vjust = 0),
 axis.title.y=element_text(size = 18, vjust = 1, angle = 90),
 strip.text.x=element_text(size = 16),
 strip.text.y=element_text(size = 16, angle = 90))
 
  
#---------------------------------------------------------------------------------------------------
#Define simulation function
#Define Parameter values - i.v. doses 

sim_flowlimited <- function(Rate,Tinf,lasttime,Q,V,CL)
{  #begin
  
  #Debug
  #Parameter values
  #Tinf <- 10
  #Rate <- 100
  #lasttime <- 20
  #Q <- 10                              
  #V <- 20
  #CL <- 5
    
  TIMEinf <- c(0,Tinf,lasttime*2)  #*2 needed to make the function work long after the last sim time
  RATEinf <- c(Rate,0,0)
    
  #Define and interpolation function that returns rate when given time - "const" give step interpolation 
  step.doseinf <- approxfun(TIMEinf, RATEinf, method = "const")
  
  #Define times
  #units are min, L, mg - mg/L = ug/ml
  TIME <- sort(unique(c(seq(from=0,to=lasttime,by=0.5),0.001,Tinf-0.001,Tinf+0.001))) 
  
               
    #Function containing differential equations for amounts in compartments (A) - see help for lsoda
          DES <- function(T, A, THETAin)
          {
              #Infusion specifications - by approxfun interpolation
              RateIn <- step.doseinf(T)
                            
              dA <- vector(len=1)
                 Cart <- RateIn/Q
              dA[1] <- (RateIn -Q*A[1] -CL*Cart)/V  
            
              list(dA,"Cart"=Cart)
          }

    #Set initial conditions - use names here to set names of output dataframe made by lsoda 
     A_0 <- c(A1=0) 

     paramlist <- c("RateIn"=RATEinf,"Q"=Q,"V"=V, "CL"=CL)

    #Run differential equation solver 
     sim.data <- lsoda(A_0, TIME, DES, paramlist)  
       
    #Process the simulated output 
      sim.data <- data.frame(sim.data)
      sim.data$Cven <- sim.data$A1
       sim.data$Cart[sim.data$time==0] <- 0
                       
    #Draw the plot
     plotobj <- NULL
     titletext <- paste("Flow-limited Organ Kinetics")
     plotobj <- ggplot(data=sim.data)
     plotobj <- plotobj + geom_line(aes(x=time, y=Cart), size=1, alpha=0.75, colour="red")
     plotobj <- plotobj + geom_line(aes(x=time, y=Cven), size=1, alpha=0.75, colour="blue")
     plotobj <- plotobj + annotate("text", x = max(sim.data$time)*0.8, y = max(sim.data$Cart)*0.8, label = paste("Volume =",V), colour = "black", size = 6)	
     plotobj <- plotobj + annotate("text", x = max(sim.data$time)*0.8, y = max(sim.data$Cart)*0.7, label = paste("Flow =",Q), colour = "black", size = 6)	
     if(CL!=0) plotobj <- plotobj + annotate("text", x = max(sim.data$time)*0.8, y = max(sim.data$Cart)*0.6, label = paste("Clearance =",CL), colour = "black", size = 6)	
     plotobj <- plotobj + scale_y_continuous("Concentration (ug/ml)")
     plotobj <- plotobj + scale_x_continuous("Time (min)")
     #plotobj <- plotobj + ggtitle(titletext) 
     print(plotobj)
         
} #end


#Test - Flow limited kinetics
  #sim_flowlimited(Rate=100,Tinf=10,lasttime=20,Q=10,V=20,CL=0)

