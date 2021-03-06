% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anova.R
\name{aovDispTable}
\alias{aovDispTable}
\title{aovDispTable}
\usage{
aovDispTable(aovObj, caption = sys.call())
}
\arguments{
\item{aovObj}{Output from aov or ezANOVA}

\item{caption}{Required for heading}
}
\description{
Display formatted ANOVA table in command window.
}
\examples{
# Example 1:
# create dataframe
dat <- createDF(nVP = 6, nTrl = 1,
                design = list("Comp" = c("comp", "incomp")))

dat <- addDataDF(dat, RT = list("Comp comp"   = c(500, 150, 100),
                                "Comp incomp" = c(520, 150, 100)))

aovObj <- aov(RT ~ Comp + Error(VP/(Comp)), dat)
aovDispTable(aovObj)

# or with ezANOVA
library(ez)
aovRT <- ezANOVA(dat, dv=.(RT), wid = .(VP), within = .(Comp), return_aov = TRUE, detailed = TRUE)
aovDispTable(aovRT)

}
