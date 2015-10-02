testdata <- data.frame(
  N=1:10,
  L=letters[1:10],
  TYPE=rep(c(1,2), each=5)
)


save(testdata, file=paste(getwd(), "data", "testdate.rda", sep="/"))
