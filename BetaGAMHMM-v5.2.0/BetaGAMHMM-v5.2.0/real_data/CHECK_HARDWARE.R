# Optional hardware-only preview. Run from the project root.
if (!file.exists(file.path("R", "00_utils.R")) || !file.exists(file.path("R", "11_hardware.R"))) {
  stop("CHECK_HARDWARE.R must be run from the project root.")
}
source(file.path("R", "00_utils.R"), local = .GlobalEnv)
source(file.path("R", "11_hardware.R"), local = .GlobalEnv)
hw <- choose_auto_workers()
print(data.frame(
  CPU = hw$cpu_name,
  physical_cores = hw$physical_cores,
  logical_cores = hw$logical_cores,
  RAM_GB = hw$ram_gb,
  selected_workers = hw$selected_workers,
  user_override = hw$user_override,
  check.names = FALSE
))
invisible(hw)
