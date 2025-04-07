library(tidyverse)
library(pomp)
library(foreach)
library(future)
library(doFuture)
library(iterators)

plan(multisession)

set.seed(1350254336)

NP = 5000; NMIF = 100; NUM_GUESSES = 400
# NP = 200; NMIF = 10; NUM_GUESSES = 40

# The code for the SEIR model is developed from https://kingaa.github.io/sbied/pfilter/model.R

source("https://kingaa.github.io/sbied/pfilter/model.R")


seir_step <- Csnippet("
  double dN_SE = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_EI = rbinom(E,1-exp(-mu_EI*dt));
  double dN_IR = rbinom(I,1-exp(-mu_IR*dt));
  S -= dN_SE;
  E += dN_SE - dN_EI;
  I += dN_EI - dN_IR;
  R += dN_IR;
  H += dN_IR;
")


seir_init <- Csnippet("
  S = nearbyint(eta*N);
  E = 0;
  I = 1;
  R = nearbyint((1-eta)*N);
  H = 0;
")

measSIR |>
  pomp(
    rprocess=euler(seir_step,delta.t=1/7),
    rinit=seir_init,
    paramnames=c("N","Beta","mu_EI","mu_IR","rho","eta"),
    statenames=c("S","E","I","R","H")
  ) -> measSEIR


fixed_params <- c(N=38000, mu_EI = 0.8, mu_IR=2., k=10)
coef(measSIR,names(fixed_params)) <- fixed_params


## LOCAL SEARCH

## What is this 'bake' function?
## See https://kingaa.github.io/sbied/pfilter/bake.html
## for an explanation.

bake(file="local_search.rds",{
  foreach(i=1:20,.combine=c,
    .options.future=list(seed=482947940)
  ) %dofuture% {
    measSIR |>
      mif2(
        Np=2000, Nmif=50,
        cooling.fraction.50=0.5,
        rw.sd=rw_sd(Beta=0.02, rho=0.02, eta=ivp(0.02)),
        partrans=parameter_trans(log="Beta",logit=c("rho","eta")),
        paramnames=c("Beta","rho","eta")
      )
  } -> mifs_local
  attr(mifs_local,"ncpu") <- nbrOfWorkers()
  mifs_local
}) -> mifs_local


bake(file="lik_local.rds",{
  foreach(mf=mifs_local,.combine=rbind,
    .options.future=list(seed=900242057)
  ) %dofuture% {
    evals <- replicate(10, logLik(pfilter(mf,Np=5000)))
    ll <- logmeanexp(evals,se=TRUE)
    mf |> coef() |> bind_rows() |>
      bind_cols(loglik=ll[1],loglik.se=ll[2])
  } -> results
  attr(results,"ncpu") <- nbrOfWorkers()
  results
}) -> results

read_csv("measles_params.csv") |>
  bind_rows(results) |>
  arrange(-loglik) |>
  write_csv("measles_params.csv")

## GLOBAL SEARCH

runif_design(
  lower=c(Beta=5,rho=0.2,eta=0),
  upper=c(Beta=80,rho=0.9,eta=1),
  nseq=NUM_GUESSES
) -> guesses

mf1 <- mifs_local[[1]]

bake(file="global_search.rds",
  dependson=guesses,{
    foreach(guess=iter(guesses,"row"), .combine=rbind,
      .options.future=list(seed=1270401374)
    ) %dofuture% {
      mf1 |>
        mif2(params=c(guess,fixed_params)) |>
        mif2(Nmif = NMIF) -> mf
      replicate(
        10,
        mf |> pfilter(Np = NP) |> logLik()
      ) |>
        logmeanexp(se=TRUE) -> ll
      mf |> coef() |> bind_rows() |>
        bind_cols(loglik=ll[1],loglik.se=ll[2])
    } -> results
    attr(results,"ncpu") <- nbrOfWorkers()
    results
  }) |>
  filter(is.finite(loglik)) -> results

read_csv("measles_params.csv") |>
  bind_rows(results) |>
  filter(is.finite(loglik)) |>
  arrange(-loglik) |>
  write_csv("measles_params.csv")

read_csv("measles_params.csv") |>
  filter(loglik>max(loglik)-20,loglik.se<2) |>
  sapply(range) -> box


## PROFILE LIKELIHOOD FOR RHO

freeze(seed=1196696958,
  profile_design(
    eta=seq(0.01,0.95,length=NUM_GUESSES/10),
    lower=box[1,c("Beta","rho")],
    upper=box[2,c("Beta","rho")],
    nprof=10, type="runif"
  )) -> guesses


mf1 <- mifs_local[[1]]
bake(file="rho_profile.rds",
  dependson=guesses,{
    foreach(guess=iter(guesses,"row"), .combine=rbind,
      .options.future=list(seed=2105684752)
    ) %dofuture% {
      mf1 |>
        mif2(params=c(guess, fixed_params),
          rw.sd=rw_sd(Beta=0.02,eta=ivp(0.02))) |>
        mif2(Nmif=NMIF,cooling.fraction.50=0.3) |>
        mif2() -> mf
      replicate(
        10,
        mf |> pfilter(Np=NP) |> logLik()) |>
        logmeanexp(se=TRUE) -> ll
      mf |> coef() |> bind_rows() |>
        bind_cols(loglik=ll[1],loglik.se=ll[2])
    } -> results
    attr(results,"ncpu") <- nbrOfWorkers()
    results
  }) -> results
