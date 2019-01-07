qq2 <- function (pvector, ...) 
{
  if (!is.numeric(pvector)) 
    stop("Input must be numeric.")
  pvector <- pvector[!is.na(pvector) & !is.nan(pvector) & !is.null(pvector) & 
                       is.finite(pvector) & pvector < 1 & pvector > 0]
  o = -log10(sort(pvector, decreasing = FALSE))
  e = -log10(ppoints(length(pvector)))
  def_args <- list(pch = 20, xlim = c(0, max(e)), ylim = c(0, 
                                                           max(o)), xlab = "Expected -log(p-value)", 
                   ylab = "Observed -log(p-value)")
  dotargs <- list(...)
  tryCatch(do.call("plot", c(list(x = e, y = o), def_args[!names(def_args) %in% 
                                                            names(dotargs)], dotargs)), warn = stop)
  abline(0, 1, col = "red")
}
