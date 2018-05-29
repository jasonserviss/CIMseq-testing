
# Function to create simulated derivative new cell types based on a rofile,
# by switching genes and adding some random noise.
createCelltypeFromProfile <- function(origProfile, nSwitched, mulNoise = 0, linNoise = 0) {
  newProfile <- origProfile
  if(nSwitched > 0) {
    switchIdx <- sample(1:(length(origProfile)), nSwitched * 2, replace = FALSE)
    from <- switchIdx[1:(length(switchIdx) / 2)]
    to <- switchIdx[-1:-(length(switchIdx) / 2)]
    newProfile[from] <- origProfile[to]
    newProfile[to] <- origProfile[from]
  }
  # Add some noise - multiplicative and additive.
  multNoise <- runif(n = length(newProfile), min = 1 - mulNoise, max = 1 + mulNoise)
  addNoise <- rbinom(n = length(newProfile), size = linNoise, prob = 0.5) - linNoise / 2
  newProfile <- newProfile * (multNoise) + addNoise
  
  # Make sure we have no negative values.
  newProfile[newProfile < 0] <- 0
  newProfile
}


