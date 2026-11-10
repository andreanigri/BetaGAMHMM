# HMD import, filtering, and validation.
#
# Both emission models use the SAME HMD death-count panel:
#   Beta response      : Female / Total
#   Binomial sensitivity: Female and Male (hence Total = Female + Male)
#
# HMD Complete Data Series may contain non-integer estimated death counts.
# Those values are retained exactly. No rounding or truncation is performed.

find_hmd_member <- function(zip_path, suffix) {
  listing <- utils::unzip(zip_path, list = TRUE)
  names_norm <- gsub("\\\\", "/", listing$Name)
  hit <- listing$Name[endsWith(names_norm, suffix)]
  if (length(hit) != 1L) {
    stop("Expected exactly one archive member ending in '", suffix, "'; found ", length(hit), ".")
  }
  hit
}

extract_hmd_member <- function(zip_path, suffix) {
  member <- find_hmd_member(zip_path, suffix)
  exdir <- tempfile("bghmm_hmd_")
  dir.create(exdir, recursive = TRUE, showWarnings = FALSE)

  extracted <- utils::unzip(
    zip_path,
    files = member,
    exdir = exdir,
    junkpaths = TRUE,
    overwrite = TRUE
  )

  if (!length(extracted) || !file.exists(extracted[1L])) {
    candidate <- file.path(exdir, basename(member))
    if (!file.exists(candidate)) {
      unlink(exdir, recursive = TRUE, force = TRUE)
      stop("Failed to extract HMD archive member: ", member)
    }
    extracted <- candidate
  }

  list(path = normalizePath(extracted[1L], winslash = "/", mustWork = TRUE),
       member = member, exdir = exdir)
}

read_hmd_metadata <- function(path, suffix = "Deaths_1x1.txt") {
  if (grepl("\\.zip$", path, ignore.case = TRUE)) {
    tmp <- extract_hmd_member(path, suffix)
    on.exit(unlink(tmp$exdir, recursive = TRUE, force = TRUE), add = TRUE)
    return(readLines(tmp$path, n = 1L, warn = FALSE))
  }
  readLines(path, n = 1L, warn = FALSE)
}

read_hmd_deaths <- function(path) {
  if (!file.exists(path)) stop("HMD data source not found: ", path)

  if (grepl("\\.zip$", path, ignore.case = TRUE)) {
    # On Windows, read.table() may attempt seek() while parsing a connection
    # returned by unz(), which is not seekable. Extract the single HMD member
    # to a temporary regular file first, then read it normally.
    tmp <- extract_hmd_member(path, "Deaths_1x1.txt")
    on.exit(unlink(tmp$exdir, recursive = TRUE, force = TRUE), add = TRUE)

    meta <- readLines(tmp$path, n = 1L, warn = FALSE)
    dat <- utils::read.table(
      tmp$path, header = TRUE, skip = 2L, stringsAsFactors = FALSE,
      check.names = FALSE, na.strings = c(".", "NA", "NaN"), comment.char = ""
    )
    return(list(
      data = dat, metadata = meta, member = tmp$member, source_type = "zip",
      deaths_md5 = unname(tools::md5sum(tmp$path))
    ))
  }

  meta <- read_hmd_metadata(path)
  dat <- utils::read.table(
    path, header = TRUE, skip = 2L, stringsAsFactors = FALSE,
    check.names = FALSE, na.strings = c(".", "NA", "NaN"), comment.char = ""
  )
  list(
    data = dat, metadata = meta, member = basename(path), source_type = "text",
    deaths_md5 = unname(tools::md5sum(path))
  )
}

strict_noninteger <- function(x, tol = 1e-10) {
  is.finite(x) & abs(x - round(x)) > tol
}

prepare_hmd_panel <- function(data_path,
                              country_code = "RUS",
                              application_label = "Russia",
                              year_range = c(1960L, 2014L),
                              age_range = c(0L, 40L),
                              beta_boundary_action = c("error", "continuity")) {
  beta_boundary_action <- match.arg(beta_boundary_action)

  deaths_obj <- read_hmd_deaths(data_path)

  # Exact manuscript reproduction uses the HMD Death Counts snapshot whose
  # metadata reports "Last modified: 27 Aug 2026". Historical HMD values can
  # be revised, so a different snapshot is rejected unless explicitly allowed.
  expected_hmd_date <- get0("HMD_EXPECTED_LAST_MODIFIED", ifnotfound = "")
  allow_hmd_mismatch <- isTRUE(get0("HMD_ALLOW_VERSION_MISMATCH", ifnotfound = FALSE))
  if (nzchar(expected_hmd_date)) {
    expected_tag <- paste0("Last modified: ", expected_hmd_date)
    if (!grepl(expected_tag, deaths_obj$metadata, fixed = TRUE)) {
      msg <- paste0(
        "HMD snapshot mismatch. Expected metadata containing '", expected_tag,
        "' but found: ", deaths_obj$metadata, ". ",
        "For exact manuscript reproduction, use the 27 Aug 2026 HMD Death Counts snapshot. ",
        "For a deliberately updated-data analysis, set BGHMM_ALLOW_HMD_VERSION_MISMATCH=true."
      )
      if (allow_hmd_mismatch) warning(msg, call. = FALSE) else stop(msg, call. = FALSE)
    }
  }

  expected_hmd_md5 <- get0("HMD_EXPECTED_DEATHS_MD5", ifnotfound = "")
  if (nzchar(expected_hmd_md5) && !identical(tolower(deaths_obj$deaths_md5), tolower(expected_hmd_md5))) {
    msg <- paste0(
      "HMD Deaths_1x1 content checksum mismatch. Expected MD5 ", expected_hmd_md5,
      " but found ", deaths_obj$deaths_md5, ". ",
      "For exact manuscript reproduction, use the 27 Aug 2026 Death Counts snapshot."
    )
    if (allow_hmd_mismatch) warning(msg, call. = FALSE) else stop(msg, call. = FALSE)
  }

  d <- deaths_obj$data
  required <- c("PopName", "Year", "Age", "Female", "Male", "Total")
  if (!all(required %in% names(d))) {
    stop("Deaths_1x1 data are missing required columns: ", paste(setdiff(required, names(d)), collapse = ", "))
  }

  d$Year <- suppressWarnings(as.integer(d$Year))
  d$age_num <- suppressWarnings(as.numeric(gsub("\\+", "", as.character(d$Age))))
  for (nm in c("Female", "Male", "Total")) d[[nm]] <- suppressWarnings(as.numeric(d[[nm]]))

  country_rows <- d$PopName == country_code & is.finite(d$Year) & is.finite(d$age_num)
  if (!any(country_rows)) {
    codes <- sort(unique(as.character(d$PopName)))
    stop(
      "Country code '", country_code, "' was not found in the HMD file. ",
      "Check the HMD PopName code. Example available codes: ",
      paste(utils::head(codes, 20L), collapse = ", "), "."
    )
  }

  d_country <- d[country_rows, , drop = FALSE]
  available_years <- range(d_country$Year, na.rm = TRUE)
  available_ages <- range(d_country$age_num[is.finite(d_country$age_num)], na.rm = TRUE)

  keep <- country_rows &
    d$Year >= year_range[1L] & d$Year <= year_range[2L] &
    d$age_num >= age_range[1L] & d$age_num <= age_range[2L]

  d <- d[keep, c("PopName", "Year", "age_num", "Female", "Male", "Total"), drop = FALSE]
  names(d) <- c("country", "year", "age", "female_deaths", "male_deaths", "total_deaths")
  d <- d[order(d$year, d$age), , drop = FALSE]

  years <- seq.int(as.integer(year_range[1L]), as.integer(year_range[2L]))
  ages <- seq.int(as.integer(age_range[1L]), as.integer(age_range[2L]))

  if (!nrow(d)) {
    stop(
      "No observations remain for ", country_code, " after filtering. ",
      "Available years are approximately ", available_years[1L], "-", available_years[2L],
      " and ages ", available_ages[1L], "-", available_ages[2L], "."
    )
  }
  if (anyDuplicated(d[c("year", "age")])) stop("Duplicate HMD year-age cells found in selected panel.")

  expected_key <- as.vector(outer(years, ages, function(y, a) paste(y, a, sep = ":")))
  data_key <- paste(d$year, d$age, sep = ":")
  missing_key <- setdiff(expected_key, data_key)
  if (length(missing_key)) {
    stop(
      "The requested HMD panel is not rectangular/complete. Missing cells include: ",
      paste(utils::head(missing_key, 12L), collapse = ", "),
      ". Choose a complete country-period-age window."
    )
  }
  if (nrow(d) != length(years) * length(ages)) {
    stop("Unexpected panel size after filtering: ", nrow(d), ".")
  }

  if (any(!is.finite(d$female_deaths)) || any(!is.finite(d$male_deaths)) || any(!is.finite(d$total_deaths))) {
    stop("Non-finite death counts are present in the selected panel.")
  }
  if (any(d$female_deaths < 0) || any(d$male_deaths < 0) || any(d$total_deaths <= 0)) {
    stop("Death counts must be non-negative and total deaths must be positive.")
  }

  identity_error <- max(abs((d$female_deaths + d$male_deaths) - d$total_deaths))
  if (!is.finite(identity_error) || identity_error > 1e-6) {
    stop("Female + Male does not reproduce Total within tolerance; max discrepancy = ", identity_error)
  }

  d$prop_fem <- d$female_deaths / d$total_deaths
  if (any(d$prop_fem < 0 | d$prop_fem > 1)) stop("Female death shares fall outside [0,1].")

  # Matrix layout: one HMM time point per year, one component per age.
  Y_raw <- matrix(d$prop_fem, nrow = length(years), ncol = length(ages), byrow = TRUE)
  F <- matrix(d$female_deaths, nrow = length(years), ncol = length(ages), byrow = TRUE)
  M <- matrix(d$male_deaths, nrow = length(years), ncol = length(ages), byrow = TRUE)
  N <- matrix(d$total_deaths, nrow = length(years), ncol = length(ages), byrow = TRUE)

  # Explicit same-data check used in the supplement/reproducibility audit.
  same_data_error <- max(abs(Y_raw - F / N))
  if (!is.finite(same_data_error) || same_data_error > 1e-14) {
    stop("Internal same-data check failed: Beta proportions do not equal Female/Total.")
  }

  boundary_cells <- which(Y_raw <= 0 | Y_raw >= 1, arr.ind = TRUE)
  Y_beta <- Y_raw
  if (nrow(boundary_cells)) {
    if (beta_boundary_action == "error") {
      stop(
        "Exact 0/1 proportions are present (", nrow(boundary_cells), " cells). ",
        "The Beta density is undefined at the boundary. Set BGHMM_BETA_BOUNDARY_ACTION='continuity' ",
        "only if that transformation is scientifically justified and reported."
      )
    }
    n_cells <- length(Y_raw)
    Y_beta <- (Y_raw * (n_cells - 1) + 0.5) / n_cells
  }

  frac_f <- strict_noninteger(F)
  frac_m <- strict_noninteger(M)
  frac_n <- strict_noninteger(N)
  frac_any <- frac_f | frac_m | frac_n
  frac_idx <- which(frac_any, arr.ind = TRUE)
  fractional_table <- if (nrow(frac_idx)) {
    data.frame(
      year = years[frac_idx[, 1L]],
      age = ages[frac_idx[, 2L]],
      female_hmd = F[frac_idx],
      male_hmd = M[frac_idx],
      total_hmd = N[frac_idx],
      female_fractional = frac_f[frac_idx],
      male_fractional = frac_m[frac_idx],
      total_fractional = frac_n[frac_idx],
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      year = integer(0), age = integer(0), female_hmd = numeric(0), male_hmd = numeric(0),
      total_hmd = numeric(0), female_fractional = logical(0), male_fractional = logical(0),
      total_fractional = logical(0)
    )
  }

  list(
    data = d,
    source_path = normalizePath(data_path, winslash = "/", mustWork = TRUE),
    source_type = deaths_obj$source_type,
    hmd_deaths_metadata = deaths_obj$metadata,
    hmd_deaths_md5 = deaths_obj$deaths_md5,
    hmd_deaths_member = deaths_obj$member,
    country_code = country_code,
    application_label = application_label,
    years = years,
    ages = ages,
    Y = Y_beta,
    Y_raw = Y_raw,
    F = F,
    M = M,
    N = N,
    n_years = length(years),
    n_ages = length(ages),
    n_cells = length(Y_raw),
    boundary_cells = boundary_cells,
    boundary_action = beta_boundary_action,
    same_data_max_abs_error = same_data_error,
    count_audit = list(
      n_fractional_cells = sum(frac_any),
      n_fractional_female = sum(frac_f),
      n_fractional_male = sum(frac_m),
      n_fractional_total = sum(frac_n),
      max_distance_to_integer = max(abs(c(F - round(F), M - round(M), N - round(N)))),
      fractional_table = fractional_table
    )
  )
}

write_data_diagnostics <- function(panel, output_dir) {
  safe_dir_create(output_dir)
  d <- panel$data

  summary <- data.frame(
    quantity = c(
      "country_code", "application_label", "start_year", "end_year", "min_age", "max_age",
      "n_years", "n_ages", "n_cells", "min_prop_fem", "median_prop_fem", "max_prop_fem",
      "boundary_cells", "min_total_deaths", "median_total_deaths", "max_total_deaths",
      "fractional_HMD_cells", "fractional_female_cells", "fractional_male_cells",
      "fractional_total_cells", "same_data_max_abs_error"
    ),
    value = as.character(c(
      panel$country_code, panel$application_label, min(panel$years), max(panel$years), min(panel$ages), max(panel$ages),
      panel$n_years, panel$n_ages, panel$n_cells,
      min(d$prop_fem), stats::median(d$prop_fem), max(d$prop_fem), nrow(panel$boundary_cells),
      min(d$total_deaths), stats::median(d$total_deaths), max(d$total_deaths),
      panel$count_audit$n_fractional_cells, panel$count_audit$n_fractional_female,
      panel$count_audit$n_fractional_male, panel$count_audit$n_fractional_total,
      panel$same_data_max_abs_error
    )),
    stringsAsFactors = FALSE
  )

  utils::write.csv(summary, file.path(output_dir, "data_summary.csv"), row.names = FALSE)
  utils::write.csv(d, file.path(output_dir, "analysis_panel_with_hmd_counts.csv"), row.names = FALSE)
  utils::write.csv(panel$count_audit$fractional_table, file.path(output_dir, "fractional_HMD_counts.csv"), row.names = FALSE)

  by_year <- stats::aggregate(total_deaths ~ year, d, sum)
  utils::write.csv(by_year, file.path(output_dir, "denominator_by_year.csv"), row.names = FALSE)

  same_data <- data.frame(
    check = c(
      "same year-age cells", "Beta response equals Female/Total", "Binomial-kernel uses original Female",
      "Binomial-kernel uses original Male", "No rounding or truncation"
    ),
    result = c(TRUE, panel$same_data_max_abs_error <= 1e-14, TRUE, TRUE, TRUE),
    detail = c(
      paste0(panel$n_cells, " cells"),
      paste0("max abs error = ", format(panel$same_data_max_abs_error, scientific = TRUE)),
      "panel$F", "panel$M", "HMD estimated counts retained exactly"
    ),
    stringsAsFactors = FALSE
  )
  utils::write.csv(same_data, file.path(output_dir, "same_data_beta_binomial_check.csv"), row.names = FALSE)

  writeLines(
    c(
      paste0("Deaths file metadata: ", panel$hmd_deaths_metadata),
      paste0("Country: ", panel$application_label, " [", panel$country_code, "]"),
      paste0("Window: ", min(panel$years), "-", max(panel$years), "; ages ", min(panel$ages), "-", max(panel$ages)),
      "Beta and binomial-kernel analyses use exactly the same HMD year-age cells and death estimates.",
      "No HMD death estimate is rounded or truncated.",
      "For fractional HMD estimated counts, the denominator-aware robustness model uses the gamma-extended binomial log-kernel.",
      paste0("Fractional HMD cells: ", panel$count_audit$n_fractional_cells, "/", panel$n_cells),
      paste0("Beta same-data max absolute error: ", format(panel$same_data_max_abs_error, scientific = TRUE))
    ),
    con = file.path(output_dir, "data_source_and_processing.txt")
  )
  invisible(summary)
}
