## Example files of the vignette "Length-frequency data for TropFishR"
## (inst/extdata, issue #42): the import code of the vignette works

extfile <- function(f) system.file("extdata", f, package = "TropFishR")

test_that("data arrangement 1 (one row per fish) gives valid lfq data", {
    lfq1 <- read.csv(extfile("lfq1.csv"))
    expect_named(lfq1, c("date", "length"))
    lfq1$date <- as.Date(lfq1$date, format = "%d.%m.%Y")
    expect_false(anyNA(lfq1$date))
    res <- lfqCreate(data = lfq1, Lname = "length", Dname = "date", plot = FALSE)
    expect_s3_class(res, "lfq")
    expect_equal(sum(res$catch), nrow(lfq1))
})

test_that("data arrangement 2 (frequencies, csv2 and txt) gives valid lfq data", {
    lfq2 <- read.csv2(extfile("lfq2.csv"))
    expect_named(lfq2, c("DATE", "Length_CM", "Frequency"))
    expect_true(is.numeric(lfq2$Length_CM))
    lfq2$DATE <- as.Date(lfq2$DATE, format = "%d/%m/%y")
    expect_false(anyNA(lfq2$DATE))
    res <- lfqCreate(data = lfq2, Lname = "Length_CM", Dname = "DATE",
                     Fname = "Frequency", plot = FALSE)
    expect_s3_class(res, "lfq")
    expect_equal(sum(res$catch), sum(lfq2$Frequency))
    ## tab separated text file with the same data
    lfq2a <- read.table(extfile("lfq2.txt"), header = TRUE, sep = "\t", dec = ",")
    expect_equal(lfq2a$Length_CM, lfq2$Length_CM)
    expect_equal(lfq2a$Frequency, lfq2$Frequency)
})

test_that("data arrangement 3 (length classes x dates) gives valid lfq data", {
    lfq3 <- read.csv(extfile("lfq3.csv"), check.names = FALSE)
    dates <- as.Date(colnames(lfq3)[-1], format = "%Y-%m-%d")
    expect_false(anyNA(dates))
    expect_true(all(diff(dates) > 0))
    lfq3new <- list(dates = dates, midLengths = lfq3$lengthClass,
                    catch = as.matrix(lfq3[, -1]))
    class(lfq3new) <- "lfq"
    lfq3new$catch[is.na(lfq3new$catch)] <- 0
    expect_equal(dim(lfq3new$catch), c(length(lfq3new$midLengths), length(dates)))
    expect_true(all(colSums(lfq3new$catch) > 0))
    expect_equal(count_plot_pages(plot(lfq3new, Fname = "catch")), 1)
    ## lfq3.xlsx is written from lfq3.csv (make_data/makedata_extdata_lfqData.R)
    expect_true(file.exists(extfile("lfq3.xlsx")))
})
