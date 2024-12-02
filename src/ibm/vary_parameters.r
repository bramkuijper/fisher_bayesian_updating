#!/usr/bin/env Rscript

# R script to vary parameters and generate a batch file of things

date_str <- format(Sys.time(), "%d%m%Y_%H%M%S")

a <- 1
b <- c(0.1,0.01)
c <- 0.5

biast <- 0.99

mu_p <- 0.05
mu_t <- 0.05
sdmu <- 0.4
sdmu_prior <- 0.05

output_file_prefix <- paste0("sim_fisher_bayes_",date_str)

counter <- 0

exe <-  "./fisher_bayes.exe"

nrep <- 5

batch_file_contents <- ""

for (repi in 1:nrep)
{
    for (bi in b)
    {
        counter <- counter + 1

        file_name <- paste0(output_file_prefix,"_",counter)
        command_string <- paste(exe,
                a,
                bi,
                c,
                biast,
                mu_p,
                mu_t,
                sdmu,
                sdmu_prior,
                file_name)

        batch_file_contents <- paste0(
                batch_file_contents,
                "echo ",counter,
                "\n",
                command_string,
                "\n"
                )
    }
}

writeLines(text=batch_file_contents)

