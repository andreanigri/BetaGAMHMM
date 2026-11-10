# Parse every repository R source file without executing the analyses.
files <- sort(list.files(".", pattern = "\\.R$", recursive = TRUE, full.names = TRUE))
normalized <- gsub("\\\\", "/", files)
files <- files[!grepl("(^|/)(outputs|renv|\\.git)(/|$)", normalized)]
bad <- list()
for (f in files) {
  ok <- tryCatch({parse(file = f); TRUE}, error = function(e) {bad[[f]] <<- conditionMessage(e); FALSE})
  cat(if (ok) "OK   " else "FAIL ", f, "\n", sep = "")
}
if (length(bad)) {print(bad); stop("R parse check failed.")}
cat("PASS: ", length(files), " R files parsed.\n", sep = "")
