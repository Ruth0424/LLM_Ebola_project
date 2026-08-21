## PROMPT: Generate code using odin to model an outbreak of Ebola, using the compartmental model given in this paper: “Modeling spatial invasion of Ebola in West Africa” (Jeremy P. D’Silva, Marisa C. Eisenberg (2017)). 
## INTERFACE: VS Code
## MODEL: GPT 5.4 mini

# Ebola outbreak model using odin
# Based on da Silva et al. (2017), Toward optimal control of the SEIHRD model
# with application to Ebola transmission dynamics.
#
# Compartments:
# S: susceptible
# E: exposed
# I: infectious in the community
# H: hospitalised
# R: recovered
# D: dead before burial
# B: buried

if (!requireNamespace("odin", quietly = TRUE)) {
	install.packages("odin", repos = "https://cloud.r-project.org")
}

library(odin)

ebola_dsilva_seihrdb <- odin::odin({
	initial(S) <- N0 - E0 - I0 - H0 - R0 - D0 - B0
	initial(E) <- E0
	initial(I) <- I0
	initial(H) <- H0
	initial(R) <- R0
	initial(D) <- D0
	initial(B) <- B0

	# Transmission parameters
	beta_I <- user(0.55)
	beta_H <- user(0.30)
	beta_D <- user(0.45)

	# Progression and removal rates
	sigma <- user(1 / 9.6)
	rho <- user(0.20)
	gamma_I <- user(0.15)
	gamma_H <- user(0.30)
	mu_I <- user(0.08)
	mu_H <- user(0.15)
	psi <- user(0.20)

	# Initial conditions and population size
	N0 <- user(10000)
	E0 <- user(0)
	I0 <- user(1)
	H0 <- user(0)
	R0 <- user(0)
	D0 <- user(0)
	B0 <- user(0)

	N <- S + E + I + H + R + D + B
	lambda <- (beta_I * I + beta_H * H + beta_D * D) / N

	deriv(S) <- -lambda * S
	deriv(E) <- lambda * S - sigma * E
	deriv(I) <- sigma * E - (rho + gamma_I + mu_I) * I
	deriv(H) <- rho * I - (gamma_H + mu_H) * H
	deriv(R) <- gamma_I * I + gamma_H * H
	deriv(D) <- mu_I * I + mu_H * H - psi * D
	deriv(B) <- psi * D

	output(force_of_infection) <- lambda
	output(incidence) <- lambda * S
	output(recoveries) <- gamma_I * I + gamma_H * H
	output(deaths) <- mu_I * I + mu_H * H
	output(burial_flow) <- psi * D
	output(total_infectious) <- I + H + D
	output(total_population) <- N
})

model <- ebola_dsilva_seihrdb$new(
	N0 = 10000,
	E0 = 0,
	I0 = 1,
	H0 = 0,
	R0 = 0,
	D0 = 0,
	B0 = 0,
	beta_I = 0.55,
	beta_H = 0.30,
	beta_D = 0.45,
	sigma = 1 / 9.6,
	rho = 0.20,
	gamma_I = 0.15,
	gamma_H = 0.30,
	mu_I = 0.08,
	mu_H = 0.15,
	psi = 0.20
)

times <- seq(0, 180, by = 1)
solution <- as.data.frame(model$run(times))

peak_idx <- which.max(solution$total_infectious)
cat("Peak infectious burden:", round(solution$total_infectious[peak_idx], 2),
		"on day", round(solution$t[peak_idx], 1), "\n")
cat("Final state after 180 days:\n")
print(tail(solution, 1))

if (!dir.exists("output")) {
	dir.create("output")
}

png("output/dsilva_seihrdb.png", width = 1600, height = 900, res = 180)
matplot(
	solution$t,
	solution[, c("S", "E", "I", "H", "R", "D", "B")],
	type = "l",
	lty = 1,
	lwd = 2,
	col = c("blue", "orange", "red", "purple", "green", "black", "brown"),
	xlab = "Time (days)",
	ylab = "Population",
	main = "Ebola outbreak model with da Silva et al. (2017) compartmental structure"
)
legend(
	"right",
	legend = c("Susceptible", "Exposed", "Infectious", "Hospitalised", "Recovered", "Dead", "Buried"),
	col = c("blue", "orange", "red", "purple", "green", "black", "brown"),
	lty = 1,
	lwd = 2
)
dev.off()

cat("Plot saved to output/dsilva_seihrdb.png\n")
