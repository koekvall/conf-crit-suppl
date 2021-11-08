library(tidyverse)
res_mat <- rbind(readRDS("ber_sims_gamma21131021.Rds"),
                 readRDS("ber_sims_gamma16171021.Rds"))

tidy_dat <- as_tibble(res_mat) %>%
  pivot_longer(col = 1:6,
               names_to = c(".value", "method"),
               names_sep = "_",
               values_drop_na = TRUE)
tidy_dat %>%
  group_by(method, lam, n) %>%
  summarize(prob = 1 - mean(stat > qchisq(0.95, 2)), reps = n(), se = sqrt(prob * (1 - prob) / reps)) %>%
  ggplot(aes(x = lam, y = prob, group = method, col = method)) +
  geom_line() +
  geom_errorbar(aes(ymin = prob - se, ymax = prob + se)) +
  facet_grid(. ~ n) +
  theme_bw() +
  labs(x = "sigma", y = "coverage")


tidy_dat %>%
  filter(n == 80) %>%
  group_by(method, lam) %>%
  summarize(prob = 1 - mean(stat > qchisq(0.95, 2)), reps = n(), se = sqrt(prob * (1 - prob) / reps)) %>%
  ggplot(aes(x = lam, y = prob, group = method, col = method)) +
  geom_line() +
  geom_errorbar(aes(ymin = prob - se, ymax = prob + se)) +
  theme_bw() +
  labs(x = "sigma", y = "coverage")

PDF = TRUE
out_dir <- "~/" # end in "/"
if(PDF) pdf(paste(out_dir, "fig_cover_ber.pdf", sep = ""), width = 10, height = 4)
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                "#0072B2", "#D55E00", "#CC79A7")
par(cex.axis = 1.1, cex.lab = 1.1, cex = 1.1)
par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.5,1.1))
par(mfrow = c(1, 2))
summary_dat <- tidy_dat %>%
  group_by(method, lam, n) %>%
  summarize(prob = 1 - mean(stat > qchisq(0.95, 2)), reps = n(), se = sqrt(prob * (1 - prob) / reps))
attach(summary_dat)
for(jj in c(20, 80)){
  plot(x = lam[n == jj & method == "our"],
       y = prob[n == jj & method == "our"],
       type = "l",
       ylim = c(0.92, 0.99),
       lwd = 2,
       ylab = "estimated cover probability",
       xlab = expression(lambda),
       main = paste0("n = ", jj),
       col = cbbPalette[1])
  arrows(lam[n == jj & method == "our"],
         prob[n == jj & method == "our"],
         lam[n == jj & method == "our"],
         prob[n == jj & method == "our"] + 2 *  se[n == jj & method == "our"],
         length = 0.05,
         angle = 90,
         col = cbbPalette[1])
  arrows(lam[n == jj & method == "our"],
         prob[n == jj & method == "our"],
         lam[n == jj & method == "our"],
         prob[n == jj & method == "our"] - 2 *  se[n == jj & method == "our"],
         length = 0.05,
         angle = 90,
         col = cbbPalette[1])

  abline(h = 0.95, lwd = 1)

  lines(x = lam[n == jj & method == "lrt"],
       y = prob[n == jj & method == "lrt"],
       lwd = 2,
       col = cbbPalette[2],
       lty = 2)
  arrows(lam[n == jj & method == "lrt"],
         prob[n == jj & method == "lrt"],
         lam[n == jj & method == "lrt"],
         prob[n == jj & method == "lrt"] + 2 *  se[n == jj & method == "lrt"],
         length = 0.05,
         angle = 90,
         col = cbbPalette[2])
  arrows(lam[n == jj & method == "lrt"],
         prob[n == jj & method == "lrt"],
         lam[n == jj & method == "lrt"],
         prob[n == jj & method == "lrt"] - 2 *  se[n == jj & method == "lrt"],
         length = 0.05,
         angle = 90,
         col = cbbPalette[2])

  lines(x = lam[n == jj & method == "wald"],
        y = prob[n == jj & method == "wald"],
        lwd = 2,
        col = cbbPalette[3],
        lty = 3)
  arrows(lam[n == jj & method == "wald"],
         prob[n == jj & method == "wald"],
         lam[n == jj & method == "wald"],
         prob[n == jj & method == "wald"] + 2 *  se[n == jj & method == "wald"],
         length = 0.05,
         angle = 90,
         col = cbbPalette[3])
  arrows(lam[n == jj & method == "wald"],
         prob[n == jj & method == "wald"],
         lam[n == jj & method == "wald"],
         prob[n == jj & method == "wald"] - 2 *  se[n == jj & method == "wald"],
         length = 0.05,
         angle = 90,
         col = cbbPalette[3])
}
if(PDF) dev.off()
detach(summary_dat)
