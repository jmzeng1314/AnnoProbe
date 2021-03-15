devtools::check()
rhub::check_for_cran()
rhub::check_for_cran(
  platform="windows-x86_64-devel",
  env_vars=c(R_COMPILE_AND_INSTALL_PACKAGES = "always")
)

rhub::platforms()
as.data.frame(platforms())
