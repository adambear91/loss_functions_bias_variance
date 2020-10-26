source("utilities.R", local = TRUE)

#### Loop through parameters for different losses ####
bs <- seq(0, 5, length.out = 501) # biases
ss <- seq(1, 6, length.out = 501) # sigmas
ns_power <- c(.5, 1, 2, 3)        # exponents for power loss
params <- bind_rows(
    crossing(b = bs, s = ss, n = ns_power, ltype = "power"),
    crossing(b = bs, s = ss, n = NA, ltype = c("log", "exp", "lm"))
)

future::plan("multiprocess")
results <- furrr::future_map_dbl(
    seq_len(nrow(params)),
    ~ do.call(calc_loss, c(params[., ])),
    .progress = TRUE)

results <- bind_cols(params, loss = results)

# Normalize losses for each function
results <- results %>%
    group_by(n, ltype) %>%
    mutate(loss_std = as.numeric(scale(loss))) %>% # standardize loss
    mutate(b_plus_s = b + s) %>%  # budget constraint
    group_by(b_plus_s, .add = TRUE) %>%
    mutate(tangent = ifelse(loss_std == min(loss_std), T, F)) %>%  # losses along tangent of isoclines
    ungroup()

# just power losses
results_pw <- filter(results, ltype == "power")
tangent_pw <- results_pw %>%
    filter(b_plus_s %% .25 == 0, between(b, .25, 4.75), between(s, 1.25, 5.75), tangent)

# non-power losses
results_other <- results %>%
    filter(ltype %in% c("log", "exp", "lm")) %>%
    mutate(ltype = factor(ltype, levels = c("lm", "log", "exp"), ordered = TRUE))
tangent_other <- results_other %>%
    filter(b_plus_s %% .25 == 0, between(b, .25, 4.75), between(s, 1.25, 5.75), tangent)


#### Calculate minimum loss for each beta value ####
pw_mins <- results_pw %>%
    filter(b %in% 0:5) %>%
    group_by(b, n) %>%
    mutate(is_min = loss_std == min(loss_std)) %>%
    filter(is_min) %>%
    ungroup()

nonpower_mins <- results_other %>%
    filter(b %in% 0:5) %>%
    group_by(b, ltype) %>%
    mutate(is_min = loss_std == min(loss_std)) %>%
    filter(is_min) %>%
    ungroup()


#### Plot of different loss functions ####
label_loss <- function(loss) paste("Power", loss)
xs <- seq(-10, 10, length.out = 1000)
power_losses <- map_dfc(
    list(l.5 = .5, l1 = 1, l2 = 2, l3 = 3),
    ~ abs(xs) ^ .
)
log_loss <- log(abs(xs))
exp_loss <- exp(.2 * abs(xs))
lm_loss <- -exp(-.5 * abs(xs)^2)
losses <-
    cbind(x = xs, power_losses, log_loss = log_loss, exp_loss = exp_loss, lm_loss = lm_loss) %>%
    gather(key = "loss_fn", value = "loss", -x) %>%
    group_by(loss_fn) %>%
    mutate(loss_std = as.numeric(scale(loss))) %>% # standardize
    mutate(loss_std = (loss_std - min(loss_std))) %>%
    mutate(loss_std = loss_std / max(loss_std))
losses$loss_fn <- ordered(
    losses$loss_fn, c("l.5", "l1", "l2", "l3",  "log_loss", "lm_loss", "exp_loss")
)

ggplot(losses, aes(x = x, y = loss_std)) +
    geom_line() +
    facet_wrap(~loss_fn, ncol = 2,
               labeller = as_labeller(
                   c(l.5 = "Power: .5",
                     l1 = "Power: 1",
                     l2 = "Power: 2",
                     l3 = "Power: 3",
                     log_loss = "Logarithmic",
                     exp_loss = "Exponential",
                     lm_loss = "Local Mass")
               )) +
    scale_x_continuous(name = NULL, breaks = NULL) +
    scale_y_continuous(name = NULL, breaks = NULL) +
    theme_classic() +
    theme(line = element_blank(), strip.text = element_text(size = 12))

#### Plot results ####
# Contour plot for power loss
ggplot(results_pw, aes(x = b, y = s)) +
    geom_raster(aes(fill = loss_std)) +
    geom_contour(aes(z = loss_std), color = "white", bins = 10) +
    geom_line(data = tangent_pw, size = 1.25, color = "green",
              aes(x = b, y = s)) +
    facet_wrap(vars(n), ncol = 1, strip.position = "top", labeller = labeller(.rows = label_loss)) +
    scale_fill_gradient(low = "blue", high = "red") +
    scale_x_continuous(name = NULL) +
    scale_y_continuous(name = NULL, breaks = 1:6) +
    coord_fixed() +
    theme(legend.position = "none",
          axis.text = element_text(size = 14),
          strip.text = element_blank())

# Line plot for power loss
ggplot(results_pw %>% filter(b %in% 0:5),
       aes(x = s, y = loss_std, color = factor(b))) +
    geom_line() +
    geom_point(
        aes(x = s, y = loss_std), data = pw_mins, inherit.aes = FALSE,
        color = "red", size = 1.5,
    ) +
    facet_wrap(vars(n), ncol = 1, labeller = labeller(.rows = label_loss)) +
    scale_color_discrete(name = expression(beta), h = c(120, 240)) +
    scale_x_continuous(name = NULL, breaks = 1:6) +
    scale_y_continuous(name = NULL, breaks = seq(-4, 4, by = 1)) +
    theme(legend.position = "right",
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 18),
          strip.text = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 2)))

# Contour plot for other losses
ggplot(results_other, aes(x = b, y = s)) +
    geom_raster(aes(fill = loss_std)) +
    geom_contour(aes(z = loss_std), color = "white", bins = 10) +
    geom_line(data = tangent_other, size = 1.25, color = "green",
              aes(x = b, y = s)) +
    facet_wrap(vars(ltype), ncol = 1) +
    scale_fill_gradient(low = "blue", high = "red") +
    scale_x_continuous(name = NULL) +
    scale_y_continuous(name = NULL, breaks = 1:6) +
    coord_fixed() +
    theme(legend.position = "none",
          axis.text = element_text(size = 14),
          strip.text = element_blank())

# Line plot for other losses
ggplot(results_other %>% filter(b %in% 0:5),
       aes(x = s, y = loss_std, color = factor(b))) +
    geom_line() +
    geom_point(
        aes(x = s, y = loss_std), data = nonpower_mins, inherit.aes = FALSE,
        color = "red", size = 1.5,
    ) +
    facet_wrap(vars(ltype), ncol = 1) +
    scale_color_discrete(name = expression(beta), h = c(120, 240)) +
    scale_x_continuous(name = NULL, breaks = 1:6) +
    scale_y_continuous(name = NULL, breaks = seq(-4, 4, by = 1)) +
    theme(legend.position = "right",
          axis.text = element_text(size = 14),
          legend.text = element_text(size = 18),
          strip.text = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 2)))
