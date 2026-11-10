# Verify the local HMD Death Counts archive without fitting a model.
.check_hmd_data <- function() {
  this_ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  this_file <- if (!is.null(this_ofile) && length(this_ofile) && nzchar(this_ofile[1L])) {
    tryCatch(normalizePath(this_ofile[1L], winslash = "/", mustWork = TRUE), error = function(e) NA_character_)
  } else NA_character_
  root <- if (length(this_file) == 1L && !is.na(this_file)) dirname(this_file) else getwd()
  oldwd <- getwd(); on.exit(setwd(oldwd), add = TRUE)
  setwd(file.path(root, "real_data"))
  source(file.path("config", "real_data_config.R"), local = .GlobalEnv)
  source(file.path("R", "01_data.R"), local = .GlobalEnv)
  if (!file.exists(DATA_FILE)) stop("HMD archive not found: ", DATA_FILE, "\nSee docs/DATA.md.", call. = FALSE)
  meta <- read_hmd_metadata(DATA_FILE)
  cat("HMD archive: ", normalizePath(DATA_FILE, winslash = "/", mustWork = TRUE), "\n", sep = "")
  cat("Metadata: ", meta, "\n", sep = "")
  expected <- paste0("Last modified: ", HMD_EXPECTED_LAST_MODIFIED)
  if (!grepl(expected, meta, fixed = TRUE)) stop("Unexpected HMD snapshot; expected ", expected, call. = FALSE)
  tmp <- extract_hmd_member(DATA_FILE, "Deaths_1x1.txt")
  on.exit(unlink(tmp$exdir, recursive = TRUE, force = TRUE), add = TRUE)
  md5 <- unname(tools::md5sum(tmp$path))
  cat("Deaths_1x1 MD5: ", md5, "\n", sep = "")
  if (!identical(tolower(md5), tolower(HMD_EXPECTED_DEATHS_MD5))) {
    stop("Unexpected HMD Deaths_1x1 checksum; expected ", HMD_EXPECTED_DEATHS_MD5, call. = FALSE)
  }
  cat("PASS: expected HMD snapshot and Deaths_1x1 checksum detected.\n")
  invisible(TRUE)
}
.check_hmd_data()
rm(.check_hmd_data)
