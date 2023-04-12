##############################################
##############################################
##########SimulateX############################
###############################################

#Credit is due, this was taken from 
#https://stackoverflow.com/questions/14967813/is-there-a-function-or-package-which-will-simulate-predictions-for-an-object-ret

simulateX <- function(object, seed = NULL, X, ...) {
  nsim = 1 
  object$fitted.values <- predict(object, X)
  simulate(object = object, nsim = nsim, seed = seed, ...)
}


############################################
###############################################
############MoveSimp##########################
###############################################

moveSimp<-function(pt1,pt2, RSS1, RSS2,YagiMod,WhipMod, mylev = 0.80, exp_predictions = TRUE,
                   ANT1, ANT2){
  ##Predictions for first point using mylev
  if(ANT1 == "Yagi"){
  r1<-predict(YagiMod, newdat = data.frame(Power = RSS1), interval = "prediction",
                  level = mylev)[3]} else {
                    r1<-predict(WhipMod, newdat = data.frame(Power = RSS1), interval = "prediction",
                                level = mylev)[3]    
                  }
  if(ANT2 == "Yagi"){
  r2<-predict(YagiMod, newdat = data.frame(Power = RSS2), interval = "prediction",
                  level = mylev)[3] } else { 
                    r2<-predict(WhipMod, newdat = data.frame(Power = RSS2), interval = "prediction",
                                level = mylev)[3]     }
  
  if(exp_predictions == TRUE){
    r1<-exp(r1)
    r2<-exp(r2)
  }
  buff1<-st_buffer(pt1, r1)##add buffer
  buff2<-st_buffer(pt2, r2)##add buffer

  #plot(rbind(buff1,buff2))
  #do they overlap
  moveSimpResult1<-TRUE #default to true
  
  
  if(length(st_intersection(st_geometry(buff1), st_geometry(buff2))) == 1 ){
    moveSimpResult1 <- FALSE}

  
  
  
  
  ######Do same thing but use the predicted dat Not upper cprediction inteval
  
  
  if(ANT1 == "Yagi"){
    r1<-predict(YagiMod, newdat = data.frame(Power = RSS1))} else {
                  r1<-predict(WhipMod, newdat = data.frame(Power = RSS1))   
                }
  if(ANT2 == "Yagi"){
    r2<-predict(YagiMod, newdat = data.frame(Power = RSS2)) } else { 
                  r2<-predict(WhipMod, newdat = data.frame(Power = RSS2))    }
  
  if(exp_predictions == TRUE){
    r1<-exp(r1)
    r2<-exp(r2)
  }
 
  buff1<-st_buffer(pt1, r1)##add buffer
  buff2<-st_buffer(pt2, r2)##add buffer
  
  #plot(rbind(buff1,buff2))
  #do they overlap
  moveSimpResult2<-TRUE #default to true
  
  
  if(length(st_intersection(st_geometry(buff1), st_geometry(buff2))) == 1 ){
    moveSimpResult2 <- FALSE}
  
  
  
  mydf<-data.frame(DIST = as.numeric(st_distance(pt1, pt2)), 
                   MOVE1 = moveSimpResult1, 
                   MOVE = moveSimpResult2)
  colnames(mydf)[2]<-as.character(paste("MOVE_", mylev, sep =""))
  return(mydf)
}




##############################################
###############################################
#############MAake New Point#################
##############################################

MakeNewPoint<-function(X, #The UTM X coordinate
                       Y, #The UTM Y coordinate
                       YagiMod, #The RSS Model for Yagi
                       WhipMod, #The RSS model for Whip
                       RSS, #The RSS observed during the location estimate
                       ANT, #The antenna type used
                       maxR,#Maximum distance allowed 
                       nsim,#The number of simulations to do
                       exp_predictions = TRUE){
  
  newX <- data.frame(Power = rep(RSS,nsim))
  
  if(ANT == "Whip"){mymod<-WhipMod} else { mymod<-YagiMod}
  
  
  newrs<-simulateX(mymod,  X = newX)#the new simulated responses
  
  #If they need to be exponentiated then do this
  if(exp_predictions == TRUE){
    newrs<-exp(newrs)
  }
  
  #check if they are above max allowed then change
  newrs[newrs$sim_1 > maxR, "sim_1"]<-maxR
  
  #generate an angle
  a<-runif(nsim,0, 2*pi)
  
  #Compute new point
  NewX<- X + newrs$sim_1 * cos(a)
  NewY<-Y + newrs$sim_1 * sin(a)
  return(data.frame(X = NewX,Y = NewY))
         
}









##############################################
###############################################
########MC method along the river#################
##############################################


moveMC<-function(X1,X2,Y1, Y2,#THe points
                 RSS1, RSS2, #The RSS values of the points
                 ANT1, ANT2, #The antennas used for the poitns
                 YagiMod, WhipMod, #RSS models
                 nsim, #Number of simulations to do
                 maxR, #max allowed value for radius
                 exp_predictions = TRUE){ #should we exponentiate model predictions (see above)

P1<-MakeNewPoint(X1,Y1,YagiMod,WhipMod, RSS = RSS1, ANT = ANT1, maxR = 1000, nsim = nsim)
P2<-MakeNewPoint(X2,Y2,YagiMod,WhipMod, RSS = RSS2, ANT = ANT2, maxR = 1000, nsim = nsim)    

P1_snap<-xy2segvert(x = P1$X, y = P1$Y, riv)
P1$seg<-P1_snap$seg
P1$vert<-P1_snap$vert

P2_snap<-xy2segvert(x = P2$X, y = P2$Y, riv)
P2$seg<-P2_snap$seg
P2$vert<-P2_snap$vert


#compute distances and put on P1
P1$Distance<-NA
for(i in 1:nsim){

P1$Distance[i]<-upstream(startseg = P1$seg[i],
                   endseg = P2$seg[i],
                  startvert = P1$vert[i],
                 endvert = P2$vert[i],
                rivers = riv)
}

#hist(P1$Distance)
P1$US_DS<-"DS"
P1[P1$Distance > 0,"US_DS"]<-"US"
#hist(P1$Distance, breaks = 50)
#table(P1$US_DS)

return(data.frame(nsim = nsim,
           US_p = 1 - nrow(P1[P1$US_DS == "US",])/nsim,
           DS_p = nrow(P1[P1$US_DS == "US",])/nsim,
           dist_ABS_mean = mean(abs(P1$Distance)),
           dist_ABS_sd = sd(abs(P1$Distance)),
           dist_mean = mean(P1$Distance),
           dist_sd = sd(P1$Distance)
           ))


}
                 
