test_that("date2yeardec gives the fraction of the year elapsed", {
    expect_equal(date2yeardec(as.Date("2021-01-01")), 2021)
    ## 2020 is a leap year: 1 March is day 60 (0-based) of 366
    expect_equal(date2yeardec(as.Date("2020-03-01")), 2020 + 60 / 366)
    expect_equal(date2yeardec(as.Date("2021-07-02")), 2021 + 182 / 365)
})

test_that("yeardec2date is the inverse of date2yeardec", {
    expect_equal(yeardec2date(date2yeardec(as.Date("2015-01-01"))),
                 as.Date("2015-01-01"))
    dates <- as.Date(c("2015-01-01", "2016-02-29", "2019-07-15", "2020-12-31"))
    expect_equal(yeardec2date(date2yeardec(dates)), dates)
})

test_that("date2yeardec rejects non-date input", {
    expect_error(date2yeardec("2020-01-01"), "not in POSIXt or Date")
})
