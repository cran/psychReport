#' @title mathString
#'
#' @description Returns formatted string following addition/subtraction.
#'
#' @param str1 string
#' @param str2 string
#' @param operation "add" vs. "subtract"
#' @param numDigits number 0 (default)
#' @param unit "ms" vs. "mv" vs. "\%"
#'
#' @return NULL
#'
#' @examples
#' # Example 1:
#' string <- mathString("550 ms", "480 ms")
#'
#' # Example 2:
#' string <- mathString("2.34", "1.65", "add", numDigits = 2, unit = "mv")
#'
#' \dontrun{
#' # Example use in *.Rnw Sweave file
#' # \Sexpr{string} }
#'
#' @export
mathString <- function(str1, str2, operation = "subtract",
                       numDigits = 0, unit = "ms") {

  num1 <- as.numeric(regmatches(str1, gregexpr("[[:digit:]]+\\.*[[:digit:]]*", str1)))
  num2 <- as.numeric(regmatches(str2, gregexpr("[[:digit:]]+\\.*[[:digit:]]*", str2)))
  if (operation == "add") {
    result = num1 + num2
  } else if (operation == "subtract") {
    result = num1 - num2
  }

  return(numValueString(result, numDigits, unit))

}
