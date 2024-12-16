# The  curv function is used to iteratively apply the curveball algorithm to the polygon by species matrix 
# This function takes an adjacency matrix and number of iterations as inputs and returns a dataframe 

#we have to make the constraint the networks to connected networks

curv <- function(mat, n){
  ran.poly<-list()
  for(i in 1:n){
    newmat <- curve_ball(mat)
    ran.poly[[i]]<- newmat
  }	
  ran.poly
}
