### Title:  Nimble function for if else statement in area correction
### Author: Abbey Feuka
### Date: 11JUL22
### Notes:
###################
library(nimble)
#double(0) <- scalar input/output
#double(1) <- vector input/output
avail_fun <- nimbleFunction(
  run = function(removal = double(0), gamma = double(0), 
                 gamma.past = double(1), theta = double(0),
                 theta.past = double(1)) {
    returnType(double(0))
    if(removal!=1) {
      return(gamma*theta*prod((1-gamma.past)+gamma.past*(1-theta.past)))
    } else {
      return(gamma*theta)
    }
  })
