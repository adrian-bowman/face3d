#     Invoke all the tests

test.prompt <- FALSE
reinstall   <- FALSE

library(sm)
devtools::install("face3d")

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
