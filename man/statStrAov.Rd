% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anova.R
\name{statStrAov}
\alias{statStrAov}
\title{statStrAov}
\usage{
statStrAov(aovObj, effect)
}
\arguments{
\item{aovObj}{Output from aov or ezANOVA (NB. ezANOVA must be called with detailed = TRUE)}

\item{effect}{The effect required from the anova table}
}
\description{
Returns Latex formatted string from ANOVA required for R/knitr integration.
For example, \deqn{F(1, 20) = 8.45, p < 0.01, pes = 0.45}
Returns values to 2 sig decimal places and < 0.01, < 0.001 for p values.
}
\examples{
# Example 1:
# create dataframe and add data with 2(Comp: comp vs. incomp) and 2(Side: left vs. right)
dat <- createDF(nVP = 20, nTrl = 1,
                design = list("Comp" = c("comp", "incomp"),
                              "Side" = c("left", "right")))

dat <- addDataDF(dat, RT = list("Comp:Side comp:left"    = c(500, 150, 100),
                                "Comp:Side comp:right"   = c(500, 150, 100),
                                "Comp:Side incomp:left"  = c(520, 150, 100),
                                "Comp:Side incomp:right" = c(520, 150, 100)))

aovRT <- aov(RT ~ Comp*Side + Error(VP/(Comp*Side)), dat)
aovRT <- aovTable(aovRT)

aovString <- statStrAov(aovRT, "Comp")
aovString <- statStrAov(aovRT, "Comp:Side")


# or using ezANOVA
library(ez)
aovRT <- ezANOVA(dat, dv=.(RT), wid = .(VP), within = .(Comp, Side),
                 return_aov = TRUE, detailed = TRUE)
aovRT <- aovTable(aovRT)

aovString <- statStrAov(aovRT, "Comp")
aovString <- statStrAov(aovRT, "Comp:Side")

}
