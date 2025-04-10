### Login to Greatlake
`ssh uniqname@greatlakes.arc-ts.umich.edu`

### Moving Files
From local to remote:
`scp myfile uniqname@greatlakes-xfer.arc-ts.umich.edu:mydir`
`scp -r mydir uniqname@greatlakes-xfer.arc-ts.umich.edu:`

It also works the other way round.

### Batch Jobs
Create a new job:
`sbatch sample.sbat`

Query Job Status
`squeue -j jobid`
`squeue -u uniqname`

Delete a job:
`scancel jobid`

Check job state and est. start time:
`sbatch --test-only sample.sbat`

#### More Slurm commands:
Show job hist:
`sacct -u user`

Show CPU utilization for job:
`seff jobid`

List accounts with permission:
`my_accounts`

### R Modules
Load default R distro:
`module load R`

See other available distros:
`module avail R`


