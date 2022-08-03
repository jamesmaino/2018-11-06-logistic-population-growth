library(tidyverse) 
library(deSolve)

parameters <- c(r = 0.3, K = 100)

state <- c(N = 2)

logistic <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
        dN <- (1 - N / K) * (r * N)
        list(dN)
    })
}

times <- seq(0, 50, by = 0.2)

sol <- as.data.frame(
    ode(y = state, times = times, func = logistic, parms = parameters)
)

ggplot(sol) +
    geom_line(aes(time, N), col = "red") +
    theme_minimal() +
    theme(
        plot.background = element_rect(fill = rgb(.2, .21, .27)),
        text = element_text(colour = "grey"),
        axis.text = element_text(colour = "grey"),
        panel.grid = element_line(colour = "grey")
    )
ggsave("plots/ode1.png", width=8, height=5)

r <- 0.3

K <- 100

N0 <- 2

time <- seq(0, 50, by = 0.2)

sol <- data.frame(
    time,
    N = K * N0 / (N0 + (K - N0) * exp(-r * time))
)

ggplot(sol) +
    geom_line(aes(time, N), col = "red") +
    theme_minimal() +
    theme(
        plot.background = element_rect(fill = rgb(.2, .21, .27)),
        text = element_text(colour = "grey"),
        axis.text = element_text(colour = "grey"),
        panel.grid = element_line(colour = "grey")
    )
ggsave("plots/ode2.png", width=8, height=5)

