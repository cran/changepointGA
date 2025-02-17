## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(changepointGA)

## -----------------------------------------------------------------------------
Ts = 200
betaT = c(0.5) # intercept
XMatT = matrix(rep(1, Ts), ncol=1)
colnames(XMatT) = c("intercept")
sigmaT = 1
phiT = c(0.5)
thetaT = c(0.8)
DeltaT = c(2, -2)
CpLocT = c(50, 150)

myts = ts.sim(beta=betaT, XMat=XMatT, sigma=sigmaT, phi=phiT, theta=thetaT, Delta=DeltaT, CpLoc=CpLocT, seed=1234)
str(myts)

## ----fig.align = "center", fig.height=4, fig.width=6--------------------------
TsPlotCheck(Y=myts, tau=CpLocT)

## -----------------------------------------------------------------------------
ARIMA.BIC.Order(chromosome=c(2, 1, 1, 50, 150, Ts+1), plen=2, XMat=XMatT, Xt=myts)

## -----------------------------------------------------------------------------
N = Ts
p.range = list(ar=c(0,2), ma=c(0,2))

GA_param = list(
  popsize      = 80,
  Pcrossover   = 0.95,
  Pmutation    = 0.3,
  Pchangepoint = 10/N,
  minDist      = 1,
  mmax         = N/2 - 1,
  lmax         = 2 + N/2 - 1,
  maxgen       = 5000,
  maxconv      = 100,
  option       = "both",
  monitoring   = FALSE,
  parallel     = FALSE,
  nCore        = NULL,
  tol          = 1e-5,
  seed         = 1234
)

## -----------------------------------------------------------------------------
pop = matrix(0, nrow=2 + N/2 - 1, ncol=80)
# put the two desired candidate solutions into pop matrix
pop[1:6,1] = c(2, 1, 1, 50, 150, N+1)
pop[1:7,2] = c(3, 1, 1, 50, 100, 150, N+1)

# fill the remain 38 individuals slot from default generation
tmppop = random_population(popsize=80, prange=p.range, N=N, minDist=1, Pb=10/N, 
                           mmax=N/2 - 1, 2 + N/2 - 1)
pop[,3:80] = tmppop[,3:80]
pop[1:10, 1:10]

operators.self = list(population = pop,
                      selection  = "selection_linearrank",
                      crossover  = "uniformcrossover",
                      mutation   = "mutation")

## -----------------------------------------------------------------------------
XMatEst = matrix(1, nrow=N, ncol=1)

## -----------------------------------------------------------------------------
res.changepointGA = suppressWarnings(GA(ObjFunc=ARIMA.BIC.Order, N=N, GA_param, GA_operators = operators.self, 
                                        p.range=p.range, XMat=XMatEst, Xt=myts))
fit.changepointGA = res.changepointGA$overbestfit
fit.changepointGA
tau.changepointGA = res.changepointGA$overbestchrom
tau.changepointGA

## -----------------------------------------------------------------------------
IslandGA_param = list(
  subpopsize   = 80,
  Islandsize   = 2,
  Pcrossover   = 0.95,
  Pmutation    = 0.15,
  Pchangepoint = 10/N,
  minDist      = 1,
  mmax         = N/2 - 1,
  lmax         = 2 + N/2 - 1,
  maxMig       = 1000,
  maxgen       = 50,
  maxconv      = 20,
  option       = "both",
  monitoring   = FALSE,
  parallel     = FALSE, ###
  nCore        = NULL,
  tol          = 1e-5,
  seed         = 1234
)
tim1 = Sys.time()
res.Island.changepointGA = suppressWarnings(IslandGA(ObjFunc=ARIMA.BIC.Order, N=N, IslandGA_param, 
                                                     p.range=p.range, XMat=XMatEst, Xt=myts))
tim2 = Sys.time()
tim2 - tim1
fit.Island.changepointGA = res.Island.changepointGA$overbestfit
fit.Island.changepointGA
tau.Island.changepointGA = res.Island.changepointGA$overbestchrom
tau.Island.changepointGA

## -----------------------------------------------------------------------------
IslandGA_param = list(
  subpopsize   = 80,
  Islandsize   = 2,
  Pcrossover   = 0.95,
  Pmutation    = 0.15,
  Pchangepoint = 10/N,
  minDist      = 1,
  mmax         = N/2 - 1,
  lmax         = 2 + N/2 - 1,
  maxMig       = 1000,
  maxgen       = 50,
  maxconv      = 20,
  option       = "both",
  monitoring   = FALSE,
  parallel     = TRUE, ###
  nCore        = 2,
  tol          = 1e-5,
  seed         = 1234
)

tim3 = Sys.time()
res.Island.changepointGA = suppressWarnings(IslandGA(ObjFunc=ARIMA.BIC.Order, N=N, IslandGA_param, 
                                                     p.range=p.range, XMat=XMatEst, Xt=myts))

tim4 = Sys.time()
tim4 - tim3
fit.Island.changepointGA = res.Island.changepointGA$overbestfit
fit.Island.changepointGA
tau.Island.changepointGA = res.Island.changepointGA$overbestchrom
tau.Island.changepointGA

## -----------------------------------------------------------------------------
true.tau = c(50, 150)
est.tau = c(tau.Island.changepointGA[4:(4+tau.Island.changepointGA[1]-1)])
cptDist(tau1=true.tau, tau2=est.tau, N=N)

## -----------------------------------------------------------------------------
sessionInfo()

