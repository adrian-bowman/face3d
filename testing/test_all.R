#     Invoke all the tests

library(face3d)
devtools::install("face3d")

library(rgl)
library(crayon)

test.prompt <- FALSE
reinstall   <- FALSE

test_label <- function(label, test.prompt) {
  cat("\n**", label, "...")
  if (test.prompt) readline(prompt = "   Press [enter] to continue ...") else cat("\n\n")
}

fls <- list.files("testing", full.names = TRUE)
fls <- fls[-match("testing/test_all.R", fls)]
for (fl in fls) {
   cat("\n******", fl, "******\n\n")
   source(fl)
}
cat("\nCompleted.\n")
