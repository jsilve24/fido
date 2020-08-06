library(testthat)
library(fido)

#Sys.setenv(KMP_DUPLICATE_LIB_OK="TRUE")
test_check("fido")
#Sys.unsetenv("KMP_DUPLICATE_LIB_OK")
