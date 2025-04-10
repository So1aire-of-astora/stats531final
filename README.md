# POMP-Based Time Series Analysis of COVID-19 in Kerala, India

## ARIMA Model

### Key Findings
+ Benchmark log likelihood: around -2,100
+ Outstanding fit
+ Diagnostics check out, except for significant heteroskedasticidy in the residuals. 

## SEIR(S) Model

### Key Findings - Local Search
+ Adding the transition from *R* to *S* improves the likelihood of the model.
+ The log likelihood reaches around -2,000, marginally outperforming the benchmark.
+ Simulations based on the local search are still shitty though. Check the graphs.

### Logging Each Run
Paul has created a Bash script that runs any given R script, collects info such as initial parameters, number of particles, etc., and writes
these info to a log file. You can run the bash script using the following command:

`$ ./rjob_runner.sh --script=your_script.R --log_file=your_logfile.log --clear_rds=true`

For example, to run the local search algorithm, simply use

`$ ./rjob_runner.sh --script=local.R --log_file=seirs_session.log --clear_rds=true`

If you see an error like

`bash: ./rjob_runner.sh: Permission denied`

or on MacOS,

`zsh: ./rjob_runner.sh: Permission denied`

you can make the script executable by running 

`chmod +x ./rjob_runner.sh`

Please do NOT edit the last flag `--clear_rds=true` for now, since cached *.rds* files may cause issues. (I will fix that later, since it's not a high priority task. --Paul) 

Note that the log will keep the records from all previous runs, unless you empty the log file manually, or just delete it.

### Acknowledgement
Some R code are borrowed from [STATS 531 course materials](https://ionides.github.io/531w25/). 
