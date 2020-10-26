library(tidyverse)
library(glmnet)
theme_set(theme_classic())

# power loss
calc_power_loss <- function(b, s, n) {
        xs <- seq(-50, 50, length.out = 1000)
        f_x <- dnorm(xs, b, s) / sum(dnorm(xs, b, s))
        as.numeric(f_x %*% (abs(xs) ^ n))
}

# log loss
calc_log_loss <- function(b, s, c = 1) {
        xs <- seq(-50, 50, length.out = 1000)
        f_x <- dnorm(xs, b, s) / sum(dnorm(xs, b, s))
        as.numeric(f_x %*% (c * log(abs(xs))))
}

# exponential loss
calc_exp_loss <- function(b, s, c = .2) {
        xs <- seq(-50, 50, length.out = 1000)
        f_x <- dnorm(xs, b, s) / sum(dnorm(xs, b, s))
        as.numeric(f_x %*% (exp(c * abs(xs))))
}


# local mass loss function (Brainard & Freeman, 1997)
calc_lm_loss <- function(b, s, c = .5) {
        xs <- seq(-50, 50, length.out = 1000)
        f_x <- dnorm(xs, b, s) / sum(dnorm(xs, b, s))
        as.numeric(f_x %*% -exp(-c * abs(xs)^2))
}


# Generic loss
calc_loss <- function(b, s, ltype, ...) {
        params <- list(...)
        if (ltype == "power") {
            return(calc_power_loss(b, s, params$n))
        } else if (ltype == "log") {
            return(calc_log_loss(b, s))
        } else if (ltype == "exp") {
            return(calc_exp_loss(b, s))
        } else if (ltype == "lm") {
            return(calc_lm_loss(b, s))
        } else if (ltype == "asym") {
            return(calc_asym_loss(b, s, params$n))
        } else {
            return(NA)
        }
}

# Modified version of colMeans function to deal with vectors
colMeansMod <- function(x) {
    if (is.matrix(x)) {
        return(colMeans(x))
    } else {
        return(x)
    }
}


# Return similarity matrix of k nearest agents based on training items
# If randomize = TRUE, similarity is determined randomly (useful for debugging)
getSimilarityMat <- function(data, agent_inds, train_inds, k, randomize = FALSE) {
    if (k >= length(agent_inds)) stop("`k` must be less than number of agents.")

    train_data <- data[agent_inds, train_inds]
    train_cors <- cor(t(train_data)) # correlations between agents
    diag(train_cors) <- NA # this ensures algorithm can't learn from itself

    if (randomize) {
        # assign k random similar individuals (excluding oneself)
        similarity_mat <- map(
            seq_len(length(agent_inds)),
            ~sample(agent_inds[-.], k)
        ) %>%
            unlist() %>%
            matrix(nrow = length(agent_inds), byrow = TRUE)
    } else {
        # create matrix with each row storing agent's k closest neighbors
        similarity_mat <- map(
            seq_len(length(agent_inds)),
            ~agent_inds[order(train_cors[., ], decreasing = TRUE)[seq_len(k)]]
        ) %>%
            unlist() %>%
            matrix(nrow = length(agent_inds), byrow = TRUE)
    }


    return(as.matrix(similarity_mat))
}

# Test recommendations based on `similarity_mat` with specified loss function
#   `power_loss` is the exponent for the power loss function
#   `epsilon` is the tolerance for a tolerance loss function
testRecs <- function(data, similarity_mat, agent_inds, test_inds, power_loss = NULL, epsilon = NULL) {
    test_data <- data[agent_inds, test_inds]

    if(!is.null(power_loss)) {
        if(!is.null(epsilon)) warning("Since `power_loss` is specified, `epsilon` is ignored.")
        # compute average agent's costs based on loss f'n
        # Cost F'n: mean(abs(pred - target)^power_loss)^(1/power_loss)
        #   where abs(pred - target)^power_loss is loss f'n
        cost <- map_dbl(
            seq_along(agent_inds),
            ~mean(abs(
                colMeansMod(test_data[match(similarity_mat[., ], agent_inds), ]) -
                    test_data[., ]
            )^power_loss)^(1/power_loss)
        ) %>%
        mean() # average over all agent's individual costs
    } else if(!is.null(epsilon)) {
        # Get mean of predictions outside of epsilon tolerance for absolute error
        cost <- map_dbl(
            seq_along(agent_inds),
            ~mean(abs(
                colMeansMod(test_data[match(similarity_mat[., ], agent_inds), ]) -
                    test_data[., ]
            ) > epsilon)
        ) %>%
            mean()
    } else {
        stop("`power_loss` or `epsilon` must be specified.")
    }

    return(cost)
}


# Fit knn model on training data and then test with specified loss or threshold
# Mode can be "power_loss" or "tolerance"
trainTest <- function(data, agent_inds, train_inds, test_inds, k, mode = "power_loss",
                      power_loss = NULL, t = NULL, random_similarity = FALSE) {
    similarity_mat <- getSimilarityMat(data, agent_inds, train_inds, k, randomize = random_similarity)

    if (mode == "power_loss") {
        cost <- testRecs(data, similarity_mat, agent_inds, test_inds, power_loss)
        results <- tibble(
            n_train = length(train_inds),
            n_test = length(test_inds),
            k = k,
            power_loss = power_loss,
            cost = cost
        )
    } else if (mode == "tolerance") {
        cost <- testRecs(data, similarity_mat, agent_inds, test_inds, epsilon = t)
        results <- tibble(
            n_train = length(train_inds),
            n_test = length(test_inds),
            k = k,
            epsilon = t,
            cost = cost
        )
    } else {
        stop("The mode you entered is not valid.")
    }

    # return a tibble so we can iterate over parameter combos
    return(results)
}

# Simulate a universe of agents & train/test items on set of parameters
simUniverse <- function(data, n_agents, n_train, n_test, mode = "power_loss",
                        power_losses = NULL, ts = NULL, random_similarity = FALSE) {
    agent_inds <- sample(seq_len(nrow(data)), n_agents)
    item_inds <- sample(seq_len(ncol(data)), n_train + n_test)
    train_inds <- item_inds[seq_len(n_train)]
    test_inds <- item_inds[-seq_len(n_train)]

    if (mode == "power_loss") {
    # parameters to loop through
        params <- crossing(k = seq_len(n_agents - 1), power_loss = power_losses)

        # get cost on all parameter combinations
        results <- map2_dfr(
            params$k,
            params$power_loss,
            ~trainTest(
                data,
                agent_inds, train_inds, test_inds, k = .x, power_loss = .y,
                random_similarity = random_similarity
            ),
        )
    } else if (mode == "tolerance") {
        params <- crossing(k = seq_len(n_agents - 1), t = ts)
        results <- map2_dfr(
            params$k,
            params$t,
            ~trainTest(
                data,
                agent_inds, train_inds, test_inds, k = .x, mode = "tolerance", t = .y,
                random_similarity = random_similarity
            ),
        )
    } else {
        stop("The mode you entered is not valid.")
    }

    return(results)
}


# Fit lasso regression and extract optimal lambda penalties
optimizeLasso <- function(data, n_agents, n_train, n_test, lambdas, power_losses = NULL, ts = NULL) {
    if (is.null(power_losses) & is.null(ts)) {
        stop("`power_losses` or `ts` must be specified.")
    }

    agent_inds <- sample(seq_len(nrow(data)), n_agents)
    item_inds <- sample(seq_len(ncol(data)), n_train + n_test)
    train_inds <- item_inds[seq_len(n_train)]
    test_inds <- item_inds[-seq_len(n_train)]

    # predict holdout agent's ratings on basis of all other agents' ratings
    # then extract optimal lasso penalty for different loss functions
    fit <- glmnet(
        t(data[agent_inds[-1], train_inds]), data[agent_inds[1], train_inds],
        alpha = 1, lambda = lambdas
    )
    predictions <- predict(fit, newx = t(data[agent_inds[-1], test_inds]), s = lambdas)
    target <- data[agent_inds[1], test_inds]

    # *minimum* lambda penalty that performs best for each loss f'n
    if (is.null(ts)) {
        opt_lambdas <- map_dbl(
            power_losses,
            ~lambdas[
                which(colMeans(abs(predictions - target)^.) == min(colMeans(abs(predictions - target)^.)))
            ] %>%
                mean()
        )

        return(
            tibble(
                n_agents = n_agents,
                n_train = n_train,
                n_test = n_test,
                power_loss = power_losses,
                opt_lambda = opt_lambdas
            )
        )
    } else {
        opt_lambdas <- map_dbl(
          ts,
          ~ lambdas[
            which(colMeans(abs(predictions - target) > .) == min(colMeans(abs(predictions - target) > .)))
          ] %>%
            mean() # average all minima (if there are multiple)
        )

        return(
            tibble(
                n_agents = n_agents,
                n_train = n_train,
                n_test = n_test,
                epsilon = ts,
                opt_lambda = opt_lambdas
            )
        )
    }
}
