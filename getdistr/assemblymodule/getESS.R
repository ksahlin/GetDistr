x = scan(file = "/Users/ksahlin/_tmp/mcmc.r", what = double(0))
y = as.mcmc(x)
ESS = effectiveSize(y)
sink(file = "/Users/ksahlin/_tmp/mcmc_ess_nr.txt")
ESS
sink()
