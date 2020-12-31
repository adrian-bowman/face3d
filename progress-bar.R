cat("Target  :", paste(rep("-", 50), collapse = ""), "\n")
cat("Progress: ")
del <- 50 / length(fls)
ipg <- 0
idl <- 0
for (i in 1:length(fls)) {
   idl <- idl + del
   inw <- round(ipg + idl) - ipg
   if (inw > 0) {
      cat(paste(rep(".", inw), collapse = ""))
      ipg <- ipg + inw
      idl <- idl - inw
   }
}
