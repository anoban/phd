library("bayesplot")
library("ggplot2")
packageVersion("bayesplot") # 1.16.0.9000 - pre-release version

saved_plots <- readRDS("./Aki/LOGO-Cmdstan-1.40/output/ppc_loo_pit_ecdf.Rds")


# https://mc-stan.org/bayesplot/articles/articles-online-only/loo-pit-correlated-tests.html
# the reason why our plots don't have the red tails that Aki suggested is that the method argument to ppc_loo_pit_ecdf() must be set to "correlated" in order to get the
# red tails - BUT method is in fact set to "correlated" in our code????
# this happens with both versions 1.15.0 and the latest release version 1.16.0

# try the pre-release version
# nope - may be the plots should have been created with the pre-release version???

for (m in names(saved_plots)) ggplot2::ggsave(filename = sprintf("./Aki/LOGO-Cmdstan-1.40/output/ppc_loo_pit_ecdf_%s.png", m),
                                              plot = saved_plots[m], height = 6, width = 16, dpi = 500, units = "in")

# works - THE PLOTS SHOULD BE CREATED AND SERIALZIED USING THE PRE-RELEASE VERSION OF bayesplot (1.16.0.9000) FOR THE RED TAILS

# all the 8 species in the LOGO escalation set
read.csv("./Aki/LOGO-Cmdstan-1.40/output/logo/logo_escalated_folds.csv")[, "species"] |> gsub(pattern = '_', replacement = ' ') |> unique()
