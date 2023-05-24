### Function for storing log probabilities in nimble
### Author: Perry de Valpine (adapted by Abbey Feuka)
### Date: 21JUN2020
### Notes: found in nimble's google group 
###########################
require(nimble)

# A proxy "sampler" that isn't really an MCMC sampler but will be called
# in the list of samplers.  It will store the sum of nodes give in control$logProbNodes
# in the target node.
storeLogProb <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    logProbNodes <- model$expandNodeNames(control[['logProbNodes']])    
  },
  run = function() {
    sumLogProb <- model$getLogProb(logProbNodes)
    model[[target]] <<- sumLogProb
  },
  methods = list(reset = function() {})
)

# A utility function to modify an mcmc configuration to use storeLogProb.
# Arguments are:
# mcmcConf: an mcmc configuration such as returned by configureMCMC
# model:    an uncompiled model object
# target:   the node name to be used for storing the sum of log probabilities
# nodes:    the nodes whose sum of log probabilities will be stored in the target node.
#           default value for nodes is all stochastic nodes in the model (except target).
configureStoreLogProb <- function(mcmcConf, model, target, nodes) {
  if(missing(nodes))
    nodes <- model$getNodeNames(stochOnly = TRUE)
  nodes <- nodes[ nodes != target ]
  mcmcConf$removeSamplers(target)
  mcmcConf$addSampler(type = 'storeLogProb', target = target, control = list(logProbNodes = nodes))
  mcmcConf
}
