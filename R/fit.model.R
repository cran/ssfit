fit.model <-
function(p, q, n, r, starting_values, h_vector, data_true, sim_data, features, n_iter, print_results = TRUE) {

#print_results: if TRUE, prints estimates as it iterates

##########

#sim_data <- function(n, phi, p) {
#	return(dataset)
#   		}

#features <- function(data) {
#	return(stat)
#	}


# p parameters
# s = 1, ..., p

# q features

# displacement

h <- list()

for (s in 1:p) {
	h[[s]] <- h_vector[s]
	}

# starting values

theta_0 <- list()

for (s in 1:p) {
	theta_0[[s]] <- starting_values[s]
	}

#################

# number of design points needed

k <- 2^ceiling(log2(p + 1))

# D: contrast matrix (part of a hadamard matrix)

library(survey)

# hadamard(n) gives a version of a Hadamard matrix of order greater than n, but with 0's instead of -1's
# 2*hadamard(2)-1 gives a 4x4 Hadamard matrix

if (p == 1) {
	D <- matrix(c(1, -1), nrow = 2, ncol = 1)
	}
if (p > 1) {
	D <- 2 * hadamard(k - 1) - 1
	D <- D[,2:(p + 1)]
	}

#################

sim <- list()
z <- list()

results <- matrix(NA, nrow = n_iter, ncol = p)

# k: number of design points
# d = 1, ..., k

# k x 1 vector of ones

k_ones <- matrix(1, nrow = k, ncol = 1)

design.points <- list()

# initialise table of z_data (observed statistics (features))

z_data <- matrix(NA, nrow = 1, ncol = q)

for (j in 1:q) {
	z_data[1, j] <- features(data_true)[j]
	}

# ITERATE

for (iter in 1:n_iter) {

# initialise matrices with simulations at each of the design points
# initialise tables of z's

for (d in 1:k) {
	sim[[d]] <- list()
	z[[d]] <- matrix(NA, nrow = r, ncol = q)
	}


# simulate at each of the design points

for (d in 1:k) {
	design.points[[d]] <- list()
	for (s in 1:p) {
		design.points[[d]][[s]] <- theta_0[[s]] + (D[d, s] * h[[s]])
		}
	for (i in 1:r) {
		sim[[d]][[i]] <- sim_data(n, design.points[[d]])
		}
	}

# features: function that calculates the statistics / features

# calculate z's

for (d in 1:k) {
	for (i in 1:r) {
		z[[d]][i,] <- features(sim[[d]][[i]])
		}
	}

# add a N(0, 0.0001^2) if feature is zero

#for (d in 1:k) {
#	for (i in 1:r) {
#		for (j in 1:q) {
#			if (z[[d]][i, j] == 0) {
#				z[[d]][i, j] <- rnorm(1, 0, 0.0001)
#				}
#			}
#		}
#	}

# vector of ones

r_ones <- matrix(1, nrow = r, ncol = 1)

# zbar: 1 x q matrices with averages over simulations

zbar <- list()

for (d in 1:k) {
	zbar[[d]] <- (t(r_ones) %*% z[[d]]) / r
	}

# ztilde: r x q matrices with deviations from column average

ztilde <- list()
sigma <- list()

for (d in 1:k) {
	ztilde[[d]] <- matrix(NA, nrow = r, ncol = q)
	for (i in 1:r) {
		for (j in 1:q) {
			ztilde[[d]][i, j] <- z[[d]][i, j] - zbar[[d]][1, j]
			}
		}
	sigma[[d]] <- (t(ztilde[[d]]) %*% ztilde[[d]]) / (r-1)
	}

W <- t(zbar[[1]])

for (d in 2:k) {
	W <- cbind(W, t(zbar[[d]]))
	}

# V: matrix of differences

V <- W %*% D

# sigma_vv: covariance matrix of V (q x q)

sum1 <- sigma[[1]]
for (d in 2:k) {
	sum1 <- sum1 + sigma[[d]]
	}
sigma_vv <- sum1 / k

# get inverse of sigma_vv and calculate l's and lamda

c <- 1	# fix varince to 1

sigma_inv <- solve(sigma_vv)

lambda <- list()
l <- list()
y <- list()

for (s in 1:p) {
	lambda[[s]] <- sqrt((t(V[,s]) %*% sigma_inv %*% V[,s]) / c)
	l[[s]] <- (sigma_inv %*% V[,s]) / as.numeric(lambda[[s]])
	y[[s]] <- t(l[[s]]) %*% V[,s]
	}

# put l's into a matrix

L <- matrix(NA, nrow = q, ncol = p)

for (s in 1:p) {
	L[,s] <- l[[s]]
	}

ybar <- list()

for (s in 1:p) {
	ybar[[s]] <- list()
	}

for (d in 1:k) {
	for (s in 1:p) {
		ybar[[s]][[d]] <- t(l[[s]]) %*% t(zbar[[d]])
		}
	}

# create a Ybar matrix

Ybar <- matrix(NA, nrow = p, ncol = k)

for (s in 1:p) {
	for (d in 1:k) {
		Ybar[s, d] <- ybar[[s]][[d]]
		}
	}

# y_data: t(l) * z_data

y_data <- list()

for (s in 1:p) {
	y_data[[s]] <- t(l[[s]]) %*% t(z_data)
	}
	
# create column vector Y_data

Y_data <- matrix(NA, nrow = p, ncol = 1)

for (s in 1:p) {
	Y_data[s, 1] <- y_data[[s]]
	}

A <- (2 / k) * (t(V) %*% L)

for (u in 1:dim(A)[1]) {
	for (v in 1:dim(A)[2]) {
		A[u, v][u != v] <- 0
		}
	}

H <- matrix(0, nrow = p, ncol = p)

for (s in 1:p) {
	H[s, s] <- h[[s]]
	}

# estimates of parameters

theta_hat <- list()

theta_hat <- matrix(NA, nrow = p, ncol = 1)
th_0 <- matrix(unlist(theta_0))
Y_d <- matrix(unlist(y_data))

theta_hat <- th_0 + (2 * H %*% solve(A) %*% (Y_d - ((Ybar %*% k_ones) / k)))

for (s in 1:p) {
	theta_0[[s]] <- theta_hat[s]
	}

results[iter,] <- t(theta_hat)

if (print_results == TRUE) {
	print(results[iter,])
	}

}	# end for


# calculate variance

var_theta <- 4 * H %*% solve(A) %*% t(L) %*% sigma_vv %*% L %*% H %*% solve(A)

return(list("estimates"=results, "var_estimates"=var_theta, "L"=L, "sigma"=sigma_vv))

}
