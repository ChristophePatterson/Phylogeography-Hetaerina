### R script to define gamma-distributed priors for the theta and tau parameters in G-PhoCS analyses.
### By Ivan Prates, June 2018 modified by Christophe Patterson for Hetaerina damselflies

# Given a mutation rate of:
u = 2.8*10^-9 # Broad use mutation rate used in other insect papers
# 1. For tau:

# The tree height in number of generations is:
gen_het = 33.08*10^6# Estimated divergence time between H. amer/clav and H. titia from Standings et al
gen_Amer_calv = 3.76*10^6 # Estimated divegence time between H. amer and H. calv from Standing et al 
gen_titia = 3.5*10^6# Estimated divergence time of H. titia
gen_titia_alt = 1.1*10^6# Estimate divergence time between North and South Atlantic populations of H. titia
gen_Amer = 2.2*10^6 # Estimated divergence times between sub populations of H. amer

# The expected tree height (in coalescent units) is:
r = u * gen_het 
gen_Amer_calv_r = gen_Amer_calv * u
gen_titia_r = gen_titia * u
gen_titia_alt_r = gen_titia_alt *u 
gen_Amer_r = gen_Amer * u

r_all <- c(r,gen_Amer_calv_r,gen_titia_r,gen_titia_alt_r,gen_Amer_r)

r_all*1/((u)/10^-6)

# 2. For theta:
N = 1000000 
N = 10000000
theta = 4*N*u
theta # From 0.018 to 0.044

# Then, implement gamma parameters as:
shape = c(1,2,5,10)
rate = c(20,50,200)

# Plot:
par(mfrow = c(length(rate), 1)) # number of plot rows = number of rate values
for (j in 1:length(rate)) {
  for (i in 1:length(shape)) { # shapes go together, so in inside loop
    x.max = qgamma(0.999, shape=shape[i], rate=rate[j])
    x = seq(from=0, to=x.max, by=x.max/1000)
    dens = dgamma(x, shape=shape[i], rate=rate[j])
    palette = viridis(n = length(shape))
    plot(x, dens, type='l', ylim = c(0,50), xlim = c(0.003,0.5), col=palette[i], lwd = 2, xlab = "theta")
    title(paste("beta =", rate[j])) # titles for each plot
    abline(v = r_all)
    legend(x = 0.35, y = 30 - 3*i, paste("alpha =", shape[i]), text.col=palette[i], bg = "transparent", bty = "n") # settings and position of legend
    par(new = T) # plot in same graph
  }
  par(new = F) # to stop plotting on top
}
par(new = F) # to stop plotting on top
