# 8_rolling_ES_plot.R

# clear
rm(list = ls())

# specify p
k = 2
n = 10
p = 0.025

# load MN(k)-GRACH(1,1) n-step rolling parameter estimation results
load(paste0("MNG_",k,"_roll_par_",n,"_step","_data.RData"))

# load MN(k)-GRACH(1,1) n-step rolling simulated aggregated returns
load(paste0("MNG_",k,"_R_",n,"_step_simulation","_data.RData"))

# load sources
source("1_source_functions.R")

# parallel computing
no_cores = detectCores()-1
cl = makeCluster(no_cores)
clusterExport(cl = cl, varlist = ls(), envir = environment())
clusterEvalQ(cl = cl, source("1_source_functions.R"))

ES_JSU = pbsapply(1:nrow(rolling.param), function(y) 
  VaR_ES_JSU(p, n, as.vector(rolling.param[y,]), r.mat[,y], output = 2), 
  cl = cl)

ES_Sim = pbsapply(1:nrow(rolling.param), function(y) 
  VaR_ES_Sim(p, R.mat_tn[,y], output = 2),
  cl = cl)

stopCluster(cl)

range = "2008/2009"
ES_JSU_range = xts(ES_JSU, order.by = index(rolling.param))[range]
ES_Sim_range = xts(ES_Sim, order.by = index(rolling.param))[range]

# plot
df = data.frame(date = index(rolling.param[range]), 
                JSU = as.vector(ES_JSU_range), Sim = as.vector(ES_Sim_range))
df = gather(data = df, key = "variable", value = "value", -date)

dpi = 1000
tiff(paste0("ES_",n,".tif"), width = 7*dpi, height = 3.5*dpi, res = dpi)

title = paste0("Estimated 2.5% ES of rolling ",n,"-step aggregated returns")

ggplot(df, aes(x = date, y = value)) +
  geom_line(aes(color = variable, alpha = variable)) +
  scale_size_manual(values = c(0.01,0.01,0.01)) +
  scale_alpha_manual(values = c(0.7,0.7,0.7)) +
  ylim(min(df$value), 0) +
  scale_x_date(breaks = function(x)
    seq.Date(from = as.Date("2008-01-01"),
             to   = as.Date("2010-01-01"), by = "1 years"),
    date_labels = "%Y",
    limits = as.Date(c("2008-01-01","2010-01-01"))) +
  theme_economist() +
  scale_colour_economist() +
  labs(title = title, x = element_blank(), y = element_blank())+
  theme(text = element_text(size = 10, family = "serif", face = "italic"),
        plot.title   = element_text(size = 14),
        legend.title = element_blank())

model.plot = recordPlot()
dev.off()