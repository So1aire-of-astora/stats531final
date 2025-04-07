library(tidyverse)
library(pomp)
library(foreach)
library(future)
library(doFuture)
library(iterators)

plan(multisession)

set.seed(1350254336)


source("https://kingaa.github.io/sbied/pfilter/model.R")

NP = 5000; NMIF = 100; NUM_GUESSES = 400

mifs_local = readRDS("local_search.rds")

t_loc <- attr(mifs_local,"system.time")
ncpu_loc <- attr(mifs_local,"ncpu")

mifs_local |>
  traces() |>
  melt() |>
  ggplot(aes(x=iteration,y=value,group=.L1,color=factor(.L1)))+
  geom_line()+
  guides(color="none")+
  facet_wrap(~name,scales="free_y")


results = readRDS("lik_local.rds")
t_local <- attr(results,"system.time")
ncpu_local <- attr(results,"ncpu")

pairs(~loglik+b1+b2+b3+eta+rho,data=results,pch=16)

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


results = readRDS("global_search.rds")

t_global <- attr(results,"system.time")
ncpu_global <- attr(results,"ncpu")



read_csv("measles_params.csv") |>
  bind_rows(results) |>
  filter(is.finite(loglik)) |>
  arrange(-loglik) |>
  write_csv("measles_params.csv")

read_csv("measles_params.csv") |>
  filter(loglik>max(loglik)-50) |>
  bind_rows(guesses) |>
  mutate(type=if_else(is.na(loglik),"guess","result")) |>
  arrange(type) -> all

pairs(~loglik+Beta+eta+rho, data=all, pch=16, cex=0.3,
  col=ifelse(all$type=="guess",grey(0.5),"red"))



all |>
  filter(type=="result") |>
  filter(loglik>max(loglik)-10) |>
  ggplot(aes(x=eta,y=loglik))+
  geom_point()+
  labs(
    x=expression(eta),
    title="poor man's profile likelihood"
  )

read_csv("measles_params.csv") |>
  filter(loglik>max(loglik)-20,loglik.se<2) |>
  sapply(range) -> box


## PROFILE LIKELIHOOD FOR RHO

freeze(seed=1196696958,
  profile_design(
    eta=seq(0.01,0.95,length=40),
    lower=box[1,c("Beta","rho")],
    upper=box[2,c("Beta","rho")],
    nprof=15, type="runif"
  )) -> guesses
plot(guesses)


results = readRDS("rho_profile.rds")

t_rho <- attr(results,"system.time")
ncpu_rho <- attr(results,"ncpu")

read_csv("measles_params.csv") |>
  bind_rows(results) |>
  filter(is.finite(loglik)) |>
  arrange(-loglik) |>
  write_csv("measles_params.csv")

results |>
  filter(is.finite(loglik)) -> results

pairs(~loglik+Beta+eta+rho,data=results,pch=16)

results |>
  filter(loglik>max(loglik)-10,loglik.se<1) |>
  group_by(round(rho,2)) |>
  filter(rank(-loglik)<3) |>
  ungroup() |>
  ggplot(aes(x=rho,y=loglik))+
  geom_point()+
  geom_hline(
    color="red",
    yintercept=max(results$loglik)-0.5*qchisq(df=1,p=0.95)
  )

results |>
  filter(loglik>max(loglik)-0.5*qchisq(df=1,p=0.95)) |>
  summarize(min=min(rho),max=max(rho)) -> rho_ci
