# Conservative hardware diagnostics and automatic worker selection.
# No external package is required. The decision can be overridden with
# BGHMM_N_CORES before sourcing RUN_ALL.R.

read_total_ram_gb <- function() {
  sys <- Sys.info()[["sysname"]]
  out <- NA_real_
  if (identical(sys, "Windows")) {
    script <- shQuote("[math]::Round((Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory/1GB,2)")
    cmd <- c("-NoProfile", "-Command", script)
    ans <- tryCatch(system2("powershell", cmd, stdout = TRUE, stderr = FALSE), error = function(e) character(0))
    if (length(ans)) out <- suppressWarnings(as.numeric(ans[length(ans)]))
  } else if (identical(sys, "Linux") && file.exists("/proc/meminfo")) {
    x <- readLines("/proc/meminfo", warn = FALSE)
    hit <- grep("^MemTotal:", x, value = TRUE)
    if (length(hit)) {
      kb <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", hit[1L])))
      if (is.finite(kb)) out <- kb / 1024^2
    }
  } else if (identical(sys, "Darwin")) {
    ans <- tryCatch(system2("sysctl", c("-n", "hw.memsize"), stdout = TRUE, stderr = FALSE), error = function(e) character(0))
    if (length(ans)) {
      b <- suppressWarnings(as.numeric(ans[1L]))
      if (is.finite(b)) out <- b / 1024^3
    }
  }
  if (is.finite(out) && out > 0) out else NA_real_
}

read_cpu_name <- function() {
  sys <- Sys.info()[["sysname"]]
  if (identical(sys, "Windows")) {
    ans <- tryCatch(system2(
      "powershell",
      c("-NoProfile", "-Command", shQuote("(Get-CimInstance Win32_Processor | Select-Object -First 1 -ExpandProperty Name)")),
      stdout = TRUE, stderr = FALSE
    ), error = function(e) character(0))
    if (length(ans) && nzchar(trimws(ans[1L]))) return(trimws(ans[1L]))
  }
  if (identical(sys, "Linux") && file.exists("/proc/cpuinfo")) {
    x <- readLines("/proc/cpuinfo", warn = FALSE)
    hit <- grep("model name", x, value = TRUE)
    if (length(hit)) return(trimws(sub("^[^:]+:", "", hit[1L])))
  }
  paste(Sys.info()[c("machine", "sysname")], collapse = " / ")
}

choose_auto_workers <- function(reserve_physical = 1L, cap = 12L,
                                min_ram_for_os_gb = 4,
                                ram_per_worker_gb = 1.25) {
  physical <- suppressWarnings(parallel::detectCores(logical = FALSE))
  logical <- suppressWarnings(parallel::detectCores(logical = TRUE))
  if (!is.finite(logical) || logical < 1L) logical <- 2L
  # Some Windows/R combinations cannot report physical cores. In that case,
  # use half the logical-thread count rather than treating all threads as
  # independent physical cores. This is safer for long laptop workloads.
  if (!is.finite(physical) || physical < 1L) physical <- max(1L, ceiling(logical / 2))
  ram_gb <- read_total_ram_gb()

  cpu_cap <- max(1L, min(as.integer(cap), as.integer(physical) - as.integer(reserve_physical)))
  if (is.finite(ram_gb)) {
    ram_cap <- max(1L, floor((ram_gb - min_ram_for_os_gb) / ram_per_worker_gb))
    workers <- max(1L, min(cpu_cap, ram_cap))
  } else {
    ram_cap <- NA_integer_
    workers <- cpu_cap
  }

  override <- Sys.getenv("BGHMM_N_CORES", unset = "")
  overridden <- FALSE
  if (nzchar(override)) {
    z <- suppressWarnings(as.integer(override))
    if (!is.finite(z) || z < 1L) stop("BGHMM_N_CORES must be a positive integer.")
    workers <- min(z, as.integer(logical))
    overridden <- TRUE
  }

  list(
    cpu_name = read_cpu_name(),
    os = paste(Sys.info()[c("sysname", "release")], collapse = " "),
    physical_cores = as.integer(physical),
    logical_cores = as.integer(logical),
    ram_gb = ram_gb,
    cpu_cap = cpu_cap,
    ram_cap = ram_cap,
    selected_workers = as.integer(workers),
    user_override = overridden
  )
}

write_hardware_diagnostics <- function(hw, output_dir) {
  safe_dir_create(output_dir)
  tab <- data.frame(
    quantity = names(hw),
    value = vapply(hw, function(x) paste(x, collapse = ";"), character(1)),
    stringsAsFactors = FALSE
  )
  utils::write.csv(tab, file.path(output_dir, "hardware_diagnostics.csv"), row.names = FALSE)
  writeLines(c(
    paste0("CPU: ", hw$cpu_name),
    paste0("OS: ", hw$os),
    paste0("Physical cores: ", hw$physical_cores),
    paste0("Logical cores: ", hw$logical_cores),
    paste0("Detected RAM (GB): ", ifelse(is.finite(hw$ram_gb), format(hw$ram_gb, digits = 4), "not available")),
    paste0("Selected PSOCK workers: ", hw$selected_workers),
    paste0("User override: ", hw$user_override),
    "Worker selection is deliberately conservative for long multi-start optimisation on laptops."
  ), file.path(output_dir, "hardware_diagnostics.txt"))
  invisible(tab)
}

hardware_parallel_worker <- function(x) as.numeric(x)^2

test_parallel_backend <- function(n_workers, project_root = getwd(), seed = 123L) {
  if (n_workers <= 1L) return(TRUE)
  jobs <- as.list(seq_len(n_workers))
  ans <- tryCatch(
    parallel_lapply(jobs, hardware_parallel_worker, n_cores = n_workers,
                    project_root = project_root, seed = seed),
    error = function(e) NULL
  )
  !is.null(ans) && length(ans) == n_workers &&
    identical(as.numeric(unlist(ans)), as.numeric(seq_len(n_workers)^2))
}
