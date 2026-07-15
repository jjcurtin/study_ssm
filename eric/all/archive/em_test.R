# EM script - based off of MARSSkem

# Set initial conditions from kalman filter
kf.x0 <- "x10"

# Import relevant quantities
y <- MODELobj[["data"]] # must have time going across columns
d <- MODELobj[["free"]] # D or free matrix
f <- MODELobj[["fixed"]] # f matrix
inits <- MLEobj[["start"]]
model.el <- attr(MODELobj, "par.names")
model.dims <- attr(MODELobj, "model.dims")
n <- model.dims[["data"]][1]
TT <- model.dims[["data"]][2]
m <- model.dims[["x"]][1]
Id <- list(m = diag(1, m), n = diag(1, n))
IIm <- diag(1, m) # identity matrices


