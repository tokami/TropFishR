## Evaluate `expr` without console output and with plots sent to a null
## device, so tests don't print progress text or create Rplots.pdf.
## Warnings and errors still propagate.
quiet <- function(expr) {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    res <- NULL
    utils::capture.output(res <- suppressMessages(expr))
    res
}

## Number of plot pages drawn while evaluating `expr`.
count_plot_pages <- function(expr) {
    dir <- tempfile()
    dir.create(dir)
    on.exit(unlink(dir, recursive = TRUE), add = TRUE)
    grDevices::png(file.path(dir, "page%03d.png"))
    utils::capture.output(suppressMessages(force(expr)))
    grDevices::dev.off()
    length(list.files(dir))
}

## Mark a test that documents a known bug. The test body is the correct
## expectation; remove the skip_known_bug() line once the bug is fixed.
## Set TROPFISHR_KNOWN_BUGS=true to run these tests instead of skipping them.
skip_known_bug <- function(description) {
    if (!identical(Sys.getenv("TROPFISHR_KNOWN_BUGS"), "true")) {
        skip(paste("Known bug:", description))
    }
}
