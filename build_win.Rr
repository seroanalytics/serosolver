#!/usr/bin/env Rscript
dest <- "binaries"
dir.create(dest, FALSE)
temp <- tempfile()
dir.create(temp)
owd <- setwd(temp)
ok <- callr::rcmd_safe("build", c(owd, "--no-manual", "--no-build-vignettes"))
if (ok$status != 0) {
    stop(ok$stdout)
}
tgz <- file.path(getwd(), dir())
setwd(owd)

res <- buildr::build_binaries(tgz, "builderhv.dide.ic.ac.uk", 8732, timeout = 600, dest = dest)
message("built binary ", res)
