months <- 0:60
p_preg_per_month <- c("25" = 0.25, "30" = 0.15, "35" = 0.1, "40" = 0.05, "45" = 0.01)
p_success <- unlist(lapply(
  p_preg_per_month,
  function(p) pnbinom(months, 1, p)
))

#Now we just create a data frame suitable for passing to ggplot2 …
mfr_group <- paste(
  "MFR =",
  format(p_preg_per_month, digits = 2),
  "at age",
  names(p_preg_per_month)
)

mfr_group <- factor(mfr_group, levels = mfr_group)
preg_data <- data.frame(
  months = rep.int(months, length(mfr_group)) ,
  mfr_group = rep(mfr_group, each = length(months)),
  p_success = p_success
)

#and draw the plot.
library(ggplot2)
(p <- ggplot(preg_data, aes(months, p_success, colour = mfr_group)) +
    geom_point() +
    scale_x_continuous(breaks = seq.int(0, 60, 12)) +
    scale_y_continuous(breaks = seq.int(0, 1, 0.1), limits = c(0, 1)) +
    scale_colour_discrete("Monthly fecundity rate") +
    xlab("Months") +
    ylab("Probability of conception") #+
    #opts(panel.grid.major = theme_line(colour = "grey60"))
)