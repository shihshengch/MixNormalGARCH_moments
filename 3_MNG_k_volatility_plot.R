# 3_MNG_k_volatility_plot.R

# clear
rm(list = ls())

# specify k
k = 2

# load MN(k)-GRACH(1,1) parameter estimation results
load(paste0("MNG_",k,"_data.RData"))

# load sources
source("1_source_functions.R")

# show estimation result
opt.result

# compute volatility
fitted = MNG.eps2var(opt.par, eps)
h.opt   = fitted$h
X.opt = t(fitted$X)
fitted.vol = data.frame(index(eps),sqrt(cbind(h.opt, X.opt)))
names(fitted.vol) = c("date", "vol", sprintf("j = %s", 1:k))

# plot
df = gather(data = subset(fitted.vol, select = -vol),
            key = "variable", value = "value", -date)
lambda = as.vector(MNG.obj(opt.par)$lambda)
alpha = (1-(1/k)*((1:k)-1))[order(lambda, decreasing = TRUE)]

dpi = 1000
tiff(paste0("MNG_",k,"_vol.tif"), width = 7*dpi, 
     height = 3.5*dpi, res = dpi)

title = expr(italic(paste("Evoluation of  ", sigma[jt], 
                          "  for MN(", !!k, ")-GARCH(1,1)")))

ggplot(df, aes(x = date, y = value)) +
  geom_line(aes(color = variable, alpha = variable, size = variable)) +
  scale_alpha_manual(values = alpha) +
  scale_size_manual( values = alpha) +
  scale_x_date(breaks = function(x)
    seq.Date(from = as.Date("1990-01-01"),
             to   = as.Date("2018-01-01"), by = "4 years"),
    date_labels = "%Y") +
  theme_economist() +
  scale_colour_economist() +
  labs(title = title, x = element_blank(), y = element_blank()) +
  theme(text = element_text(size = 10, family = "serif", face = "italic"),
        plot.title   = element_text(size = 14),
        legend.text  = element_text(size = 12),
        legend.title = element_blank(),
        legend.margin = margin(t = -0.3, unit = 'cm'),
        legend.position = "bottom")

model.plot = recordPlot()
dev.off()