
library(dplyr)

w <- c(3.579174, 0.1751151, 0.9455746)
w <- w / sum(w)

containment <- tibble(
	"conquest" = w,
	"high" = c(0.9, 0.09, 0.01),
	"medium" = c(0.6, 0.3, 0.1),
	"low" = c(0.34, 0.33, 0.33),
	"scs1" = c(0.38, 0.47, 0.15),
	"scs2" = c(0.43, 0.24, 0.32)
)

save(containment, file="../data/containment.rdata")
