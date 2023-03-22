# genenv.R
# 
# Colects environment variables from R, required by the R *binary*.
# They are normally set by the R *script* which is circumvented if 
# we need to run the binary under the debugger as we want to debug R
# not bash...
# 
# This script is configured as a task in 'tasks.json' and referenced within
# 'launch.json' as the 'preLaunchTask'.

env <- Sys.getenv()
envnames <- names(env)
rnames <- envnames[startsWith(envnames, "R_")]
cached_names <- rnames
ld_lib_path <- Sys.getenv("LD_LIBRARY_PATH")
if (ld_lib_path != "") {
    cached_names <- c("LD_LIBRARY_PATH", rnames)
}
writeLines(paste0(cached_names, "=", env[cached_names]), ".vscode/.env")
