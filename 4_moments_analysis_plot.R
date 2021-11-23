# 4_moments_analysis_plot.R

# clear
rm(list = ls())

# specify k
k = 2

# load MN(k)-GRACH(1,1) parameter estimation results
load(paste0("MNG_",k,"_data.RData"))

# load sources
source("1_source_functions.R")

# compute conditional one-step moments
no_cores = detectCores()-1
cl = makeCluster(no_cores)
clusterExport(cl = cl, varlist = ls(), envir = environment())
clusterEvalQ(cl = cl, source("1_source_functions.R"))
var_t1 = pbsapply(2:length(eps), 
                  function(i) Moment.R(ith = 2, n = 1, param = opt.par, 
                                       eps = eps[1:i]),
                  cl = cl)
skew_t1 = pbsapply(2:length(eps), 
                   function(i) skew(n = 1, opt.par, eps = eps[1:i]),
                   cl = cl)
kurt_t1 = pbsapply(2:length(eps),
                   function(i) kurt(n = 1, opt.par, eps = eps[1:i]),
                   cl = cl)
stopCluster(cl)

# plot volatility
df = data.frame(date = index(eps[-1]), vol = sqrt(var_t1))
df = gather(data = df, key = "variable", value = "value", -date)

dpi = 1000
tiff(paste0("one-step_volatility.tif"), width = 7*dpi, 
     height = 3.5*dpi, res = dpi)

title = "Conditional One-Step Volatility"

ggplot(df, aes(x = date, y = value)) +
  geom_line(aes(color = variable, alpha = variable)) +
  scale_size_manual(values = 0.01) +
  scale_alpha_manual(values = 0.7) +
  scale_x_date(breaks = function(x)
    seq.Date(from = as.Date("1990-01-01"),
             to   = as.Date("2018-01-01"), by = "4 years"),
    date_labels = "%Y") +
  theme_economist() +
  scale_colour_economist() +
  labs(title = title, x = element_blank(), y = element_blank())+
  theme(text = element_text(size = 10, family = "serif", face = "italic"),
        plot.title   = element_text(size = 14),
        legend.position = "none")

model.plot = recordPlot()
dev.off()

# plot skewness
df = data.frame(date = index(eps[-1]), skew = skew_t1)
df = gather(data = df, key = "variable", value = "value", -date)

dpi = 1000
tiff(paste0("one-step_skewness.tif"), width = 7*dpi, 
     height = 3.5*dpi, res = dpi)

title = "Conditional One-Step Skewness"

ggplot(df, aes(x = date, y = value)) +
  geom_line(aes(color = variable, alpha = variable)) +
  scale_size_manual(values = 0.01) +
  scale_alpha_manual(values = 0.7) +
  scale_x_date(breaks = function(x)
    seq.Date(from = as.Date("1990-01-01"),
             to   = as.Date("2018-01-01"), by = "4 years"),
    date_labels = "%Y") +
  theme_economist() +
  scale_colour_economist() +
  labs(title = title, x = element_blank(), y = element_blank())+
  theme(text = element_text(size = 10, family = "serif", face = "italic"),
        plot.title   = element_text(size = 14),
        legend.position = "none")

model.plot = recordPlot()
dev.off()

# plot kurtosis
df = data.frame(date = index(eps[-1]), kurt = kurt_t1)
df = gather(data = df, key = "variable", value = "value", -date)

dpi = 1000
tiff(paste0("one-step_kurtosis.tif"), width = 7*dpi, 
     height = 3.5*dpi, res = dpi)

title = "Conditional One-Step Kurtosis"

ggplot(df, aes(x = date, y = value)) +
  geom_line(aes(color = variable, alpha = variable)) +
  scale_size_manual(values = 0.01) +
  scale_alpha_manual(values = 0.7) +
  ylim(2.5, 7) +
  scale_x_date(breaks = function(x)
    seq.Date(from = as.Date("1990-01-01"),
             to   = as.Date("2018-01-01"), by = "4 years"),
    date_labels = "%Y") +
  theme_economist() +
  scale_colour_economist() +
  labs(title = title, x = element_blank(), y = element_blank())+
  theme(text = element_text(size = 10, family = "serif", face = "italic"),
        plot.title   = element_text(size = 14),
        legend.position = "none")

model.plot = recordPlot()
dev.off()

rm(list = setdiff(ls(), c("k", "var_t1", "skew_t1", "kurt_t1")))

# save data
save.image(paste0("MNG_",k,"_one-step_moments.RData"))