% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aovRoundDigits.R
\name{aovRoundDigits}
\alias{aovRoundDigits}
\title{aovRoundDigits}
\usage{
aovRoundDigits(ezObj, nsmall = 2)
}
\arguments{
\item{ezObj}{Output from ezANOVA}

\item{nsmall}{Number of small digits to round to within ANOVA table}
}
\value{
dataframe
}
\description{
Round digits to n decimal places in ezANOVA table
}
\examples{
library(psychReport)
requiredPackages(c("dplyr", "ez"))
# Example 1:
# create dataframe with 2(Comp: comp vs. incomp) and 2(Side: left vs. right) factors/levels
dat <- createDF(nVP = 20,
                nTrl = 50,
                design = list("Comp" = c("comp", "incomp"),
                              "Side" = c("left", "right")))

dat <- addDataDF(dat,
                 RT = list(list(c("Comp:comp", "Side:left"), vals = c(500, 150, 150)),
                           list(c("Comp:comp", "Side:right"), vals = c(500, 150, 150)),
                           list(c("Comp:incomp", "Side:left"), vals = c(500, 150, 150)),
                           list(c("Comp:incomp", "Side:right"), vals = c(500, 150, 150))))

# aggregate dat across trials
datAggVP <- dat \%>\%
    group_by(VP, Comp, Side) \%>\%
    summarize(N  = n(),
              rt = mean(RT))

# repeated measures ANOVA using ezANOVA
aovRT <- ezANOVA(datAggVP, dv=.(rt), wid = .(VP), within = .(Comp, Side),
                 return_aov = TRUE, detailed = TRUE)
aovDispTable(aovRT)
aovRT <- aovRoundDigits(aovRT, 2)
aovDispTable(aovRT)

}