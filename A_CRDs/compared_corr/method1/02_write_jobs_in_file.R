# compare the correlation based on independent groups  ------------------------------------
# create subfiles with sections to limit the ram

# cocor.indep.groups(r1.jk, r2.hm, n1, n2, alternative = "two.sided",
#                    test = "all", alpha = 0.05, conf.level = 0.95, null.value = 0,
#                    data.name = NULL, var.labels = NULL, return.htest = FALSE)

# Clean environment ------------------------------------

rm(list = ls())
gc()


# Packages ---------------------------------------------

library(cocor)
library(dplyr)
library(data.table) # call fread
library(R.utils) # to read .gz

args = commandArgs(trailingOnly=TRUE)

# Paths ---------------------------------------------
#chr_num=args[1]

# path mac
path = '/Users/dianaavalos/Programming/A_CRD_plots/1_CRD'
path_peak = '/Users/dianaavalos/Programming/A_CRD_plots/8_PEAKS'
outpath= '/Users/dianaavalos/Programming/A_CRD_plots/9_COMPARE_CORR_PEAKS'


# path baobab to be updated
path = '/home/users/a/avalosma/scratch/1_CRD'
path_peak = '/home/users/a/avalosma/scratch/8_PEAKS'
outpath= '/home/users/a/avalosma/scratch/9_COMPARE_CORR_PEAKS'
path_script = "/home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/compared_corr"
FOLDER="/home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/compared_corr"


# Load files ---------------------------------------------

type="hist"
Nneut <- 165
Nmono <- 160
Ntcell <- 94
sub <- 1
chr_num <- 5

system("touch /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/compared_corr/bash_commands.sh")
for (chr_num in seq(1,22)){
  print(chr_num)
  
  for (sub in seq(1,1000)){
    
    CORR_neut_sub_file <- file.path(outpath, paste0("neut_chr",chr_num,"_section_",sub,".RData"))
    CORR_mono_sub_file <- file.path(outpath, paste0("mono_chr",chr_num,"_section_",sub,".RData"))
    CORR_tcell_sub_file <- file.path(outpath, paste0("tcell_chr",chr_num,"_section_",sub,".RData"))
    
    
    # run jobs with system
    cat("\n",file =paste0(path_script,"/bash_commands.sh"),append=TRUE)   

 
    filename=file.path(outpath, paste0("COR_mono_vs_neut_chr",chr_num,"_section_",sub,".txt"))
    Rcommand <- paste0("Rscript ",path_script,"/worker_compare_coeff_per_chr.R ",CORR_mono_sub_file, " ",CORR_neut_sub_file,  " ",Nmono, " ", Nneut, " ", filename)
    job_name <- paste0("COR_mono_vs_neut_chr",chr_num,"_s_",sub)
    Bash_cmd <-  paste0( "wsbatch -J ", job_name, ".job --partition=shared-cpu --time=12:00:00 --mem=10000 -o ", FOLDER, "/log/", job_name, ".out -e ", FOLDER, "/log/", job_name, '.err --wrap="',Rcommand,'"')
#    cat(Bash_cmd)
    cat(Bash_cmd,file =paste0(path_script,"/bash_commands.sh"),append=TRUE)
    cat("\n",file =paste0(path_script,"/bash_commands.sh"),append=TRUE)
# system(Bash_cmd)
    
    filename=file.path(outpath, paste0("COR_mono_vs_tcell_chr",chr_num,"_section_",sub,".txt"))
    Rcommand <- paste0("Rscript ",path_script,"/worker_compare_coeff_per_chr.R ",CORR_mono_sub_file, " ", CORR_tcell_sub_file, " ", Nmono, " ", Ntcell, " ", filename)
    job_name <- paste0("COR_mono_vs_tcell_chr",chr_num,"_s_",sub)
    Bash_cmd <-  paste0( "wsbatch -J ", job_name, ".job --partition=shared-cpu --time=12:00:00 --mem=10000 -o ", FOLDER, "/log/", job_name, ".out -e ", FOLDER, "/log/", job_name, '.err --wrap="',Rcommand,'"')
 #   cat(Bash_cmd) 
    cat(Bash_cmd,file =paste0(path_script,"/bash_commands.sh"),append=TRUE)    
    cat("\n",file =paste0(path_script,"/bash_commands.sh"),append=TRUE)
    
    filename=file.path(outpath, paste0("COR_neut_vs_tcell_chr",chr_num,"_section_",sub,".txt"))
    Rcommand <- paste0("Rscript ",path_script,"/worker_compare_coeff_per_chr.R ",CORR_neut_sub_file, " ", CORR_tcell_sub_file,  " ",Nneut, " ", Ntcell, " ", filename)
    job_name <- paste0("COR_neut_vs_tcell_chr",chr_num,"_s_",sub)
    Bash_cmd <-  paste0( "wsbatch -J ", job_name, ".job --partition=shared-cpu --time=12:00:00 --mem=10000 -o ", FOLDER, "/log/", job_name, ".out -e ", FOLDER, "/log/", job_name, '.err --wrap="',Rcommand,'"')
  #  cat(Bash_cmd)    
    cat(Bash_cmd,file =paste0(path_script,"/bash_commands.sh"),append=TRUE)
    cat("\n",file =paste0(path_script,"/bash_commands.sh"),append=TRUE)
  }
}






