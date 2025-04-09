library(tidyverse)
library(pomp)
library(foreach)
library(future)
library(doFuture)
library(iterators)

registerDoFuture()

cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(cores) || cores < 1) cores <- parallel::detectCores(logical=FALSE)

plan(multisession, workers = cores)

set.seed(1350254336)

## THE FOLLOWING CODE HAS BEEN MODIFIED TO INCORPORATE TIME-VARYING BETA

ETA = .2; MU_EI = 1.; MU_IR = 1.5; K = 10

KERALA_POP = 34530000



# NP = 5000; NMIF = 100; NUM_GUESSES = 400
NP = 200; NMIF = 10; NUM_GUESSES = 40

# The code for the SEIR model is developed from https://kingaa.github.io/sbied/pfilter/model.R

covid_data = read.csv("./data/weekly_df.csv")
# covid_data = covid_data[6:119,]
# covid_data$Week_Number = seq(1, 114)
# rownames(covid_data) = NULL

seir_step <- Csnippet("

  double Beta;
  if (interval == 1) Beta = b1;
  else if (interval == 2) Beta = b2;
  else Beta = b3;

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
  I = 2000;
  R = nearbyint((1-eta)*N);
  H = 0;
")

# dmeas <- Csnippet("
#   lik = dnbinom_mu(reports,k,rho*H,give_log)+1.0e-25;"
# )

dmeas <- Csnippet("
  if (rho*H <= 0.0) {
    lik = (give_log) ? -1.0e10 : 0.0;
  } else {
    lik = dnbinom_mu(reports, k, rho*H, give_log);
  }
")

rmeas <- Csnippet("
  reports = rnbinom_mu(k,rho*H);"
)

emeas <- Csnippet("
  E_reports = rho*H;"
)

time_indicators = covariate_table(
  t = covid_data$Week_Number,
  interval = c(rep(1, 61),
                   rep(2, 40),
                   rep(3, 18)),
  times = "t")

## MODEL INIT

init_params = c(b1=6,b2=12,b3=40,rho=.5, mu_EI=1,mu_IR=.5,k=K,eta=.3,N=KERALA_POP)

covid_data |>
  select(Week_Number,reports=Confirmed) |>
  filter(Week_Number<=119) |>
  pomp(
    times="Week_Number",t0=1,
    rprocess=euler(seir_step,delta.t=1/7),
    rinit=seir_init,
    rmeasure=rmeas,
    dmeasure=dmeas,
    emeasure=emeas,
    accumvars="H",
    statenames=c("S", "E","I","R","H"),
    paramnames=c("b1","b2","b3","mu_EI","mu_IR","eta","rho","k","N"),
    params=init_params,
    covar = time_indicators
  ) -> COVID_SEIR


### Simulation

sim_df <- simulate(COVID_SEIR, nsim = 1, format = "data.frame") |>
  select(Week_Number, .id, reports) |>
  mutate(source = "Simulated")

# Real data
real_df <- covid_data |>
  select(Week_Number, Confirmed) |>
  rename(reports = Confirmed) |>
  mutate(source = "Observed")

# Combine for plotting
bind_rows(sim_df, real_df) |>
  ggplot(aes(x = Week_Number, y = reports, color = source)) +
  geom_line(linewidth = 1) +
  labs(
    title = "Simulated vs. Observed Weekly COVID-19 Cases",
    x = "Week",
    y = "Reported Cases",
    color = "Data Source"
  ) +
  theme_minimal()


### A quick sanity check
ll <- replicate(10, logLik(pfilter(COVID_SEIR, Np = 5000))) |>
  logmeanexp(se = TRUE)
cat("[SAN CHECK] loglik =", round(ll[1], 2), " | SE =", round(ll[2], 4), "\n")

# fixed_params <- c(N=KERALA_POP, mu_EI = 1., mu_IR=1.5, k=10)
# coef(measSIR,names(fixed_params)) <- fixed_params


## LOCAL SEARCH

bake(file="local_search.rds",{
  foreach(i=1:cores,.combine=c,
    .options.future=list(seed=482947940)
  ) %dopar% {
    COVID_SEIR |>
      mif2(
        Np=NP, Nmif=NMIF,
        cooling.fraction.50=0.5,
        rw.sd=rw_sd(b1=.02,b2=.02,b3=.02, rho=0.02, eta=ivp(0.02)),
        partrans=parameter_trans(log=c("b1","b2","b3"),logit=c("rho","eta")),
        paramnames=c("b1","b2","b3","rho","eta")
      )
  } -> mifs_local
  attr(mifs_local,"ncpu") <- nbrOfWorkers()
  mifs_local
}) -> mifs_local


bake(file="lik_local.rds",{
  foreach(mf=mifs_local,.combine=rbind,
    .options.future=list(seed=900242057)
  ) %dopar% {
    evals <- replicate(10, logLik(pfilter(mf,Np=NP)))
    ll <- logmeanexp(evals,se=TRUE)
    mf |> coef() |> bind_rows() |>
      bind_cols(loglik=ll[1],loglik.se=ll[2])
  } -> results
  attr(results,"ncpu") <- nbrOfWorkers()
  results
}) -> results

# read_csv("measles_params.csv") |>
#   bind_rows(results) |>
#   arrange(-loglik) |>
#   write_csv("measles_params.csv")



## THE FOLLOWING CODE IS NOT MODIFIED AS OF APR 6 2025


## GLOBAL SEARCH

# runif_design(
#   lower=c(b1=5, b2=5, b3=5, rho=0.2,eta=0),
#   upper=c(b1=80, b2=80, b3=80, rho=0.9,eta=1),
#   nseq=NUM_GUESSES
# ) -> guesses

# guesses_full <- t(apply(guesses, 1, function(row) {
#   c(row, N = KERALA_POP, mu_EI = 1., mu_IR = 1.5, k = 10)
# })) |> as.data.frame()

fixed_params <- tibble(N=KERALA_POP, mu_EI = MU_EI, mu_IR=MU_IR, k=K)

guesses <- runif_design(
  lower = c(b1 = 5, b2 = 5, b3 = 5, rho = 0.2, eta = 0),
  upper = c(b1 = 80, b2 = 80, b3 = 80, rho = 0.9, eta = 1),
  nseq = NUM_GUESSES
) |>
  as_tibble() |>
  bind_cols(fixed_params |>
      slice(rep(1, NUM_GUESSES))  # repeat to match number of rows
  )

mf1 <- mifs_local[[1]]


batch_size <- ceiling(nrow(guesses) / cores) # automatically set batch size for higher parallel computing efficiency
batched_guesses <- split(guesses, ceiling(seq_len(nrow(guesses)) / batch_size)) # This line is generated by ChatGPT since I need to look up the split function.

# bake(file="global_search.rds",
#   dependson=guesses,{
#     foreach(guess = batched_guesses, .combine=rbind,
#       .options.future=list(seed=1270401374)
#     ) %dopar% {

#       mf1 |> mif2(params=guess, Nmif = NMIF) -> mf
      
#       replicate(10, mf |> pfilter(Np = NP) |> logLik()) |> logmeanexp(se=TRUE) -> ll
      
#       mf |> coef() |> bind_rows() |> bind_cols(loglik=ll[1], loglik.se=ll[2])
#     } -> results

bake(file="global_search.rds", dependson=guesses, {
  foreach(batch = batched_guesses, .combine = bind_rows,
          .options.future = list(seed = 1270401374)) %dopar% {
    
    map_dfr(seq_len(nrow(batch)), function(i) {
      guess <- batch[i, ]

      # Convert row to named numeric vector
      params <- guess |> unlist() |> as.numeric()
      names(params) <- names(guess)

      mf <- mf1 |> mif2(params = params, Nmif = NMIF)

      ll <- replicate(10, mf |> pfilter(Np = NP) |> logLik()) |>
        logmeanexp(se = TRUE)

      mf |> coef() |> bind_rows() |> 
        bind_cols(loglik = ll[1], loglik.se = ll[2])
    })
  } -> results
    attr(results,"ncpu") <- nbrOfWorkers()
    results
  }) |> 
  filter(is.finite(loglik)) -> results


read_csv("params.csv") |>
  bind_rows(results) |>
  filter(is.finite(loglik)) |>
  arrange(-loglik) |>
  filter(loglik>max(loglik)-20,loglik.se<2) |>
  sapply(range) -> box


# ## PROFILE LIKELIHOOD FOR RHO

# freeze(seed=1196696958,
#   profile_design(
#     eta=seq(0.01,0.95,length=NUM_GUESSES/10),
#     lower=box[1,c("Beta","rho")],
#     upper=box[2,c("Beta","rho")],
#     nprof=10, type="runif"
#   )) -> guesses


# mf1 <- mifs_local[[1]]
# bake(file="rho_profile.rds",
#   dependson=guesses,{
#     foreach(guess=iter(guesses,"row"), .combine=rbind,
#       .options.future=list(seed=2105684752)
#     ) %dofuture% {
#       mf1 |>
#         mif2(params=c(guess, fixed_params),
#           rw.sd=rw_sd(Beta=0.02,eta=ivp(0.02))) |>
#         mif2(Nmif=NMIF,cooling.fraction.50=0.3) |>
#         mif2() -> mf
#       replicate(
#         10,
#         mf |> pfilter(Np=NP) |> logLik()) |>
#         logmeanexp(se=TRUE) -> ll
#       mf |> coef() |> bind_rows() |>
#         bind_cols(loglik=ll[1],loglik.se=ll[2])
#     } -> results
#     attr(results,"ncpu") <- nbrOfWorkers()
#     results
#   }) -> results
