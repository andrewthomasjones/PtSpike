# Modified expresso function
#
# Qianqian Zhu, 8/2010
#
expresso = function (afbatch, bg.correct = TRUE, bgcorrect.method = NULL, 
    bgcorrect.param = list(), normalize = TRUE, normalize.method = NULL, 
    normalize.param = list(), pmcorrect.method = NULL, pmcorrect.param = list(), 
    summary.method = NULL, summary.param = list(), summary.subset = NULL, 
    verbose = TRUE, widget = FALSE) 
{
    setCorrections <- function() {
        bioc.opt <- getOption("BioC")
        if (bg.correct) {
            if (is.null(bgcorrect.method)) {
                BGMethods <- bgcorrect.methods
            }
            else {
                BGMethods <- bgcorrect.method
            }
        }
        else {
            BGMethods <- "None"
        }
        if (normalize) {
            if (is.null(normalize.method)) {
                normMethods <- normalize.methods(afbatch)
            }
            else {
                normMethods <- normalize.method
            }
        }
        else {
            normMethods <- "None"
        }
        if (is.null(pmcorrect.method)) {
            PMMethods <- pmcorrect.methods
        }
        else {
            PMMethods <- pmcorrect.method
        }
        if (is.null(summary.method)) {
            expMethods <- generateExprSet.methods
        }
        else {
            expMethods <- summary.method
        }
	cat(c(BGMethods, normMethods))
        corrections <- expressoWidget(BGMethods, normMethods, 
            PMMethods, expMethods, bioc.opt$affy$bgcorrect.method, 
            bioc.opt$affy$normalize.method, bioc.opt$affy$pmcorrect.method, 
            bioc.opt$affy$summary.method)
        if (!is.null(corrections)) {
            if (corrections[["BG"]] != "None") {
                bgcorrect.method <<- corrections[["BG"]]
            }
            if (corrections[["NORM"]] != "None") {
                normalize.method <<- corrections[["NORM"]]
            }
            if (corrections[["PM"]] != "None") {
                pmcorrect.method <<- corrections[["PM"]]
            }
            if (corrections[["EXP"]] != "None") {
                summary.method <<- corrections[["EXP"]]
            }
        }
        else {
            stop("Aborted by user")
        }
    }
    if (widget) {
        require(tkWidgets) || stop("library tkWidgets could not be found !")
    }
    nchips <- length(afbatch)
    if (widget) {
        setCorrections()
    }
    if (verbose) {
        if (bg.correct) {
            cat("background correction:", bgcorrect.method, "\n")
        }
        if (normalize) {
            cat("normalization:", normalize.method, "\n")
        }
        cat("PM/MM correction :", pmcorrect.method, "\n")
        cat("expression values:", summary.method, "\n")
    }
    if (bg.correct) {
        if (verbose) 
            cat("background correcting...")
        afbatch <- do.call("bg.correct", c(alist(afbatch, method = bgcorrect.method), 
            bgcorrect.param))
        if (verbose) 
            cat("done.\n")
    }
    if (normalize) {
        if (verbose) 
            cat("normalizing...")
       afbatch <- do.call("normalize", c(alist(afbatch, normalize.method), 
           normalize.param))
#        afbatch <- normalize(afbatch, normalize.method, group=normalize.param$group, nor.method=normalize.param$nor.method)
        if (verbose) 
            cat("done.\n")
    }
	if (is.null(pmcorrect.method ) | is.null(summary.method)) return(afbatch) # QQ 1/22/09
    eset <- computeExprSet(afbatch, summary.method = summary.method, 
        pmcorrect.method = pmcorrect.method, ids = summary.subset, 
        summary.param = summary.param, pmcorrect.param = pmcorrect.param)
    return(eset)
}

