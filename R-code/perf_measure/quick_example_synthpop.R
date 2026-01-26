library(synthpop)

syn <- syn(iris)$syn

# with k samples
syn <- syn(iris, k = 10)$syn

# with smoothing for cont. variables
syn <- syn(iris, smoothing = "spline")$syn
syn <- syn(iris, smoothing = "density")$syn
