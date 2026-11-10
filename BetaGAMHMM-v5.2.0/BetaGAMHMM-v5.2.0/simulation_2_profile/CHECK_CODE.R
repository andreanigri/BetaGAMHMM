# Parse every R file without executing the analysis.
files <- sort(list.files(".", pattern="\\.R$", recursive=TRUE, full.names=TRUE))
bad <- list()
for (f in files) {
  ok <- tryCatch({parse(file=f); TRUE}, error=function(e){bad[[f]] <<- conditionMessage(e); FALSE})
  cat(if(ok) "OK   " else "FAIL ", f, "\n", sep="")
}
if(length(bad)){print(bad); stop("R parse check failed.")}
cat("PASS: ", length(files), " R files parsed.\n", sep="")
