# ------------------------------------------------------------------------------- #
# ------------------------------ HETEROTRPHIC ENERGY ---------------------------- #
# ------------------------------------------------------------------------------- #
# FROM HOOGENBOOM & CONNOLLY

Re <- function(u, W, v = 0.0104) {
	u * W / v
}

h_m <- function(u, W, D = 2 * 10^-5, a = 9.48, b = 0.849) {	
	a * (Re(u, W) ^ b) * D / W
}

P <- function(t) {
	h_m(u, W) * A * (P_vs(t) - P_vf)
}

P_day <- function(t) {					# Where has A (area) gone??
	h_m(u, W) * (P_vs(t) - P_vf)
}

P_vs <- function(t) {	
	
	P_max = O2(Re(u, W), x_P, alpha_P, beta_P) * P_vf/100
	R_d =   O2(Re(u, W), x_R, alpha_R, beta_R)  *P_vf/100
	
	I_k = 100
	temp = (P_max - R_d) * tanh(I(t) / I_k) + R_d
	return(temp)
}

P_vf <- 280.52

I <- function(t) {						# Is L = 12?
	L = 12	

	temp = I_max * sin(pi * t / L)^2
	
	return(temp)
}

O2 <- function(Re, x, alpha, beta) {
	(100 + x) + beta * exp(-alpha * Re/10000)
}

# ------------------------------------------------------------------------------- #
# --------------------------------- CONSTANTS ----------------------------------- #
# ------------------------------------------------------------------------------- #

# Phytosunthesis constants: Acropora nasuta
x_P <- 6.79
alpha_P <- 0.406
beta_P <- 28.26
x_R <-  -1.71
alpha_R <- 0.364
beta_R <- -17.74

# Phytosunthesis constants: Montipora
#x_P <- 1.7
#alpha_P <- 0.39
#beta_P <- 23
#x_R <- -1.4
#alpha_R <- 0.81
#beta_R <- -52

# ------------------------------------------------------------------------------- #
# ---------------------------------- INTEGRATION -------------------------------- #
# ------------------------------------------------------------------------------- #

# Do you integrate over 12 or 24 hours?  If 24, 

W <- 10						# cm
u <- 14
I_max <- 1800

int_P_day <- integrate(P_day, 0, 12)$value
int_R_day <- abs(P_day(0)) * 12

# Figure
ss <- seq(0, 12, 0.1)
plot(ss, P_day(ss), type = "l")
abline(h = 0, lty = 2)

# ------------------------------------------------------------------------------- #

W <- 10						# cm
u_vec <- seq(2, 40, 0.5)		# cm/s
E_vec <- seq(50, 1800, 50)
P_mat <- matrix(NA, length(E_vec), length(u_vec))
R_mat <- matrix(NA, length(E_vec), length(u_vec))

for (i in 1:length(E_vec)) {
	for (j in 1:length(u_vec)) {
		I_max <- E_vec[i]
		u <- u_vec[j]
		
		P_mat[i,j] <- integrate(P_day, 0, 12)$value
		R_mat[i,j] <- abs(P_day(0)) * 12  # Placeholder for Respiration
		
	}
}

contour(E_vec, u_vec, (P_mat / R_mat), axes = F, xlab = "Maximum daily irradiance", ylab = "Average flow velocity", main = paste(W, "cm"))
contour(E_vec, u_vec, (P_mat / R_mat), levels = c(1), col = "red", add = T, lwd = 2)
abline(v = 1000, lty = 2)

axis(1, at = seq(0, 1500, 500), las = 2)
axis(2, at = seq(0, 110, 20), las = 2)

# ------------------------------------------------------------------------------- #
# -------------------- MEAN and MAX YEARLY WATER VELOCITY ----------------------- #
# ------------------------------------------------------------------------------- #

#vel_dat = read.delim("wave_db_hxua.txt")
d_reef <- c(0, 4, 8, 12, 16, 20, 40, 60, 80, 100, 120)
wat_mat <- matrix(NA, length(d_reef), 3)
wat_mat_max <- matrix(NA, length(d_reef), 3)

#wat_mean = tapply(vel_dat$u[vel_dat$dist == d_reef[d]] * 0.6366197, list(vel_dat$year[vel_dat$dist == d_reef[d]]), mean)

for (d in 1:length(d_reef)) {
	
	if (d_reef[d] != 120) {
		#wat_mean = tapply(vel_dat$u[vel_dat$dist == d_reef[d]] * (2/pi), list(vel_dat$year[vel_dat$dist == d_reef[d]]), mean)
		
		wat_mean <- tapply(vel_dat$duration[vel_dat$dist == d_reef[d]] * vel_dat$u[vel_dat$dist == d_reef[d]] * (2/pi), list(vel_dat$year[vel_dat$dist == d_reef[d]]), sum) / tapply(vel_dat$duration[vel_dat$dist == d_reef[d]], list(vel_dat$year[vel_dat$dist == d_reef[d]]), sum)		
		
		wat_max <- tapply(vel_dat$u[vel_dat$dist == d_reef[d]], list(vel_dat$year[vel_dat$dist == d_reef[d]]), max)
		
	} else {
		#wat_mean1 = tapply(vel_dat$u[vel_dat$dist == 100] * (2/pi), list(vel_dat$year[vel_dat$dist == 100]), mean)
		#wat_mean2 = tapply(vel_dat$u[vel_dat$dist == 200] * (2/pi), list(vel_dat$year[vel_dat$dist == 200]), mean)

		wat_mean1 <- tapply(vel_dat$duration[vel_dat$dist == 100] * vel_dat$u[vel_dat$dist == 100] * (2/pi), list(vel_dat$year[vel_dat$dist == 100]), sum) / tapply(vel_dat$duration[vel_dat$dist == 100], list(vel_dat$year[vel_dat$dist == 100]), sum)	
		wat_mean2 <- tapply(vel_dat$duration[vel_dat$dist == 200] * vel_dat$u[vel_dat$dist == 200] * (2/pi), list(vel_dat$year[vel_dat$dist == 200]), sum) / tapply(vel_dat$duration[vel_dat$dist == 200], list(vel_dat$year[vel_dat$dist == 200]), sum)			

		wat_mean <- wat_mean1 - (wat_mean1 - wat_mean2) * 0.2
		wat_max1 <- tapply(vel_dat$u[vel_dat$dist == 100], list(vel_dat$year[vel_dat$dist == 100]), max)
		wat_max2 <- tapply(vel_dat$u[vel_dat$dist == 200], list(vel_dat$year[vel_dat$dist == 200]), max)
		wat_max <- wat_max1 - (wat_max1 - wat_max2) * 0.2
	}
	
	wat_mat[d,] <- c(sort(wat_mean)[19], sort(wat_mean)[2], sort(wat_mean)[36])
	wat_mat_max[d,] <- c(sort(wat_max)[19], sort(wat_max)[2], sort(wat_max)[36])
	
}

# ------------------------------------------------------------------------------- #
# ----------------------------- FIGURE ENVIRO GRADIENT -------------------------- #
# ------------------------------------------------------------------------------- #

ss0 <- sub2$tens_MN[sub2$distance_m == 0]
ss40 <- sub2$tens_MN[sub2$distance_m == 40]
ss80 <- sub2$tens_MN[sub2$distance_m == 80]
ss120 <- sub2$tens_MN[sub2$distance_m == 120]


I_max = 1000
par(mar = c(5, 4, 4, 5) + 0.1)
plot(d_reef, wat_mat[,1], type = "n", ylim = c(0, 14), axes = F, xlab = "Distance from crest, m", ylab = NA)
polygon(c(d_reef, rev(d_reef)), c(wat_mat[,2], rev(wat_mat[,3])), col = "grey", border = "grey")
lines(d_reef, wat_mat[,1], lty = 2)
text(90, 0.8, "Median", cex = 0.8)
polygon(c(d_reef, rev(d_reef)), c(wat_mat_max[,2], rev(wat_mat_max[,3])), col = "grey", border = "grey")
lines(d_reef, wat_mat_max[,1], lty = 2)
text(90, 2.1, "Max", cex = 0.8)

axis(1, las = 2)
axis(2, at = seq(0, 8, 2), las = 2)

polygon(c(d_reef, rev(d_reef)), c(rep(13.5, 11), rev(rep(12.5, 11))), col = "grey", border = "grey")
lines(d_reef, rep(13, 11), lty = 2)
text(90, 13.6, "Irradiance", cex = 0.8)

axis(2, at = c(11, 12, 13, 14), labels = c(0, 500, 1000, 1500), las = 2)

polygon(
	c(c(0, 40, 80, 120), rev(c(0, 40, 80, 120))), 
	4 + c(
		log10(10^6 * c(sort(ss0)[round(length(ss0) * 0.025)], sort(ss40)[round(length(ss40) * 0.025)], sort(ss80)[round(length(ss80) * 0.025)], sort(ss120)[round(length(ss120) * 0.025)])), 
		rev(log10(10^6 * c(sort(ss0)[round(length(ss0) * 0.975)], sort(ss40)[round(length(ss40) * 0.975)], sort(ss80)[round(length(ss80) * 0.975)], sort(ss120)[round(length(ss120) * 0.975)])))
	), 
	col = "grey", border = "grey")

lines(c(0, 40, 80, 120), 4 + log10(10^6 * c(sort(ss0)[round(length(ss0) * 0.5)], sort(ss40)[round(length(ss40) * 0.5)], sort(ss80)[round(length(ss80) * 0.5)], sort(ss120)[round(length(ss120) * 0.5)])), lwd = 1, lty = 2)
text(90, 9.7, "Strength", cex = 0.8)

axis(4, at = c(8, 9, 10, 11), labels = expression(10^4, 10^5, 10^6, 10^7), las = 2)
abline(v = c(0, 60, 120), col = c("red", "green", "blue"))

# ------------------------------------------------------------------------------- #
# ------------------------------ HETEROTRPHIC ENERGY ---------------------------- #
# ------------------------------------------------------------------------------- #

# Flow effect in S&J 1991 converted to units of cm-2 d-1 based on data in Palardy et al. 2005 giving feeding rates of 0.6 to 3.5 zooplankton per cm-2 h-1 for 3 species at 5 - 10cm s-1 flow. Using carbon per prey  of 0.15ug from Ribes et al. 1998

x_vals = c(2.9,5.1,6.0,6.0,7.1,10.0,10.0,14.9,15.9,37.8,44.9) # flow speeds from S&J 1991
y_vals = c(0.32,0.09,0.19,0.62,0.40,0.26,2.00,3.80,3.73,5.96,4.30) # prey capture in ug C cm-2 d-1
xy_1 = list(x=x_vals,y=y_vals) # convert x and y values into a list for use with 'nls'

f_nls2 = nls(y~(fmax*(x-a))/(sub_x + (x-a)),xy_1,start=list(fmax=10,sub_x=20,a=5))
# call to nls to fit three parameter rectangular hyperbola function to data

summary(f_nls2)

# PLOTS
#pp<-seq(0,50,length=100) # vector of values for 'smooth' model predictions 
#plot(x_vals,y_vals)
#lines(pp,predict(f_nls2,list(x=pp)))

H_max = 24.9
u_k = 21.327
a = 4.259

H_d = function(u) {
	H_max * (u - a) / (u_k + (u - a))
}

# ------------------------------------------------------------------------------- #
# ------------------------------------- LRO ------------------------------------- #
# ------------------------------------------------------------------------------- #

energy = function(x, d, kappa = 0.5) {
	W <<- 2 * sqrt((10^x)/pi) * 100
	x <<- x
	
	# Energy acquisition at position d over the gradient for colony area x
	u <<- wat_mat[d,1] * 100
	E_P = integrate(P_day, 0, 12)$value * 44.61 * 21.83 * 10^x * 365 / 1000 		# these end numbers are conversion factors
	E_R = abs(P_day(0)) * 12            * 44.61 * 21 * 10^x * 365 / 1000
	E = (E_P - E_R) #+ E_H  # kJ cm-2 d-1 -- note I moved '* 10^x * 365 / 1000' above

	# Energy acquisition at position 1 (reef crest) for colony area x
	u <<- wat_mat[1,1] * 100
	E_P = integrate(P_day, 0, 12)$value * 44.61 * 21.83 # these end numbers are conversion factors
	E_R = abs(P_day(0)) * 12            * 44.61 * 21
	E_crest = (E_P - E_R) * 10^x * 365 / 1000  # kJ cm-2 d-1

	# absolute energy allocation considering reproductive maturity
	if (W < 7) {
		E_growth_abs = E_P - E_R # kJ cm-2 d-1
		E_repro_abs  = 0
	} else {
		E_growth_abs = kappa * E
		E_repro_abs  = E - (kappa * E)

		#E_growth_abs = kappa * E_P - E_R		# SEAN ADD FOR ARC PROPOSAL
		#E_repro_abs  = (1 - kappa) * E_P #E - (kappa * E)		
	}
	
	ini_area = 10^x
	ini_r = sqrt(ini_area / pi)	
	inc_area_abs = (E_growth_abs / 1) / 10000 #was 1 before sabatical...? 2.56  kJ, to grow a cm^2
	inc_r_abs = sqrt((ini_area + inc_area_abs)/pi) - ini_r
	new_area_abs = log10(ini_area + inc_area_abs)
	#fecundity_abs = 10^(1.03 + 1.28 * log10((10^x) * 10000)) * E_repro_abs / ((1.23 * 10^x) * 10000)
	fecundity_abs = 672 * E_repro_abs / 1.24 #before sabbatical 672 and 1.24 ## HERE

	return(c(E, E_crest, inc_r_abs, inc_area_abs, new_area_abs, fecundity_abs))
}

hall = function(x, a = 1.03, b = 1.28) {
	10^(a + b * log10((10^x) * 10000))/(9*10^-2)
}

stimson = function(x) {
	0.1
}

### ### ### ### ### 
### TEST FIGURES ### 
### ### ### ### ### 

par(mfcol = c(3, 3))

store_c = c()
for (x in seq(-3, 0.7, 0.01)) {
	store_c = rbind(store_c, c(log10_area = x, energy(x, 1, 0.1)))
}
store_m = c()
for (x in seq(-3, 0.7, 0.01)) {
	store_m = rbind(store_m, c(log10_area = x, energy(x, 8, 0.1)))
}
store_b = c()
for (x in seq(-3, 0.7, 0.01)) {
	store_b = rbind(store_b, c(log10_area = x, energy(x, 11, 0.1)))
}

plot(store_c[,1], log10(store_c[,4]), type = "l", ylim = c(-4, 0), xlab = "colony size, log10, m^2", ylab = "added radius, log10, m", main = "kappa = 0.1", col = "red")
abline(v = log10(pi*(7/200)^2), col = "grey")
lines(store_m[,1], log10(store_m[,4]), lty = 1, col = "orange")
lines(store_b[,1], log10(store_m[,4]), lty = 1, col = "black")

plot(store_c[,1], log10(store_c[,5]), type = "l", ylim = c(-4, 0.5), xlab = "colony size, log10, m^2", ylab = "added area, log10, m", col = "red")
abline(v = log10(pi*(7/200)^2), col = "grey")
lines(store_m[,1], log10(store_m[,5]), lty = 1, col = "orange")
lines(store_b[,1], log10(store_b[,5]), lty = 1, col = "black")

plot(store_c[,1], log10(store_c[,7]), type = "l", ylim = c(2, 7), xlab = "colony size, log10, m^2", ylab = "volume of gametes, log10, mm^3", col = "red")
abline(v = log10(pi*(7/200)^2), col = "grey")
lines(store_m[,1], log10(store_m[,7]), lty = 1, col = "orange")
lines(store_b[,1], log10(store_b[,7]), lty = 1, col = "black")

###

store_c = c()
for (x in seq(-3, 0.7, 0.01)) {
	store_c = rbind(store_c, c(log10_area = x, energy(x, 1, 0.5)))
}
store_m = c()
for (x in seq(-3, 0.7, 0.01)) {
	store_m = rbind(store_m, c(log10_area = x, energy(x, 8, 0.5)))
}
store_b = c()
for (x in seq(-3, 0.7, 0.01)) {
	store_b = rbind(store_b, c(log10_area = x, energy(x, 11, 0.5)))
}

plot(store_c[,1], log10(store_c[,4]), type = "l", ylim = c(-4, 0), xlab = "colony size, log10, m^2", ylab = "added radius, log10, m", main = "kappa = 0.5", col = "red")
abline(v = log10(pi*(7/200)^2), col = "grey")
lines(store_m[,1], log10(store_m[,4]), lty = 1, col = "orange")
lines(store_b[,1], log10(store_m[,4]), lty = 1, col = "black")

plot(store_c[,1], log10(store_c[,5]), type = "l", ylim = c(-4, 0.5), xlab = "colony size, log10, m^2", ylab = "added area, log10, m", col = "red")
abline(v = log10(pi*(7/200)^2), col = "grey")
lines(store_m[,1], log10(store_m[,5]), lty = 1, col = "orange")
lines(store_b[,1], log10(store_b[,5]), lty = 1, col = "black")

plot(store_c[,1], log10(store_c[,7]), type = "l", ylim = c(2, 7), xlab = "colony size, log10, m^2", ylab = "volume of gametes, log10, mm^3", col = "red")
abline(v = log10(pi*(7/200)^2), col = "grey")
lines(store_m[,1], log10(store_m[,7]), lty = 1, col = "orange")
lines(store_b[,1], log10(store_b[,7]), lty = 1, col = "black")

###

store_c = c()
for (x in seq(-3, 0.7, 0.01)) {
	store_c = rbind(store_c, c(log10_area = x, energy(x, 1, 0.9)))
}
store_m = c()
for (x in seq(-3, 0.7, 0.01)) {
	store_m = rbind(store_m, c(log10_area = x, energy(x, 8, 0.9)))
}
store_b = c()
for (x in seq(-3, 0.7, 0.01)) {
	store_b = rbind(store_b, c(log10_area = x, energy(x, 11, 0.9)))
}

plot(store_c[,1], log10(store_c[,4]), type = "l", ylim = c(-4, 0), xlab = "colony size, log10, m^2", ylab = "added radius, log10, m", main = "kappa = 0.9", col = "red")
abline(v = log10(pi*(7/200)^2), col = "grey")
lines(store_m[,1], log10(store_m[,4]), lty = 1, col = "orange")
lines(store_b[,1], log10(store_m[,4]), lty = 1, col = "black")

plot(store_c[,1], log10(store_c[,5]), type = "l", ylim = c(-4, 0.5), xlab = "colony size, log10, m^2", ylab = "added area, log10, m", col = "red")
abline(v = log10(pi*(7/200)^2), col = "grey")
lines(store_m[,1], log10(store_m[,5]), lty = 1, col = "orange")
lines(store_b[,1], log10(store_b[,5]), lty = 1, col = "black")

plot(store_c[,1], log10(store_c[,7]), type = "l", ylim = c(2, 7), xlab = "colony size, log10, m^2", ylab = "volume of gametes, log10, mm^3", col = "red")
abline(v = log10(pi*(7/200)^2), col = "grey")
lines(store_m[,1], log10(store_m[,7]), lty = 1, col = "orange")
lines(store_b[,1], log10(store_b[,7]), lty = 1, col = "black")

# ------------------------------------------------------------------------------- #
# -------------------------------- KAPPA ---------------------------------------- #
# ------------------------------------------------------------------------------- #

bg_s = 0.9
storm_f = 1
I_max = 1800
d_reef = c(0, 4, 8, 12, 16, 20, 40, 60, 80, 100, 120)
sub_str = c(rep(ss0, 7), ss0 - (ss0 - ss80)/2, ss80, ss80 - (ss80 - ss120)/2, ss120)
#sub_str = rep(0.3, length(d_reef))

keep_lro_abs = c()

#for (kappa in c(0.5)) {
for (kappa in c(0.05, 0.25, 0.5, 0.75, 0.95)) {
		
	f_mat = matrix(NA, length(d_reef), 50)
	g_mat = matrix(NA, length(d_reef), 50)
	x_mat = matrix(NA, length(d_reef), 51)
	n_mat = matrix(NA, length(d_reef), 51)
	
	for (d in 1:length(d_reef)) {
		if (d_reef[d] != 120) {
			wat_max = tapply(vel_dat$u[vel_dat$dist == d_reef[d]], list(vel_dat$year[vel_dat$dist == d_reef[d]]), max)
		} else {
			wat_max1 = tapply(vel_dat$u[vel_dat$dist == 100], list(vel_dat$year[vel_dat$dist == 100]), max)
			wat_max2 = tapply(vel_dat$u[vel_dat$dist == 200], list(vel_dat$year[vel_dat$dist == 200]), max)
			wat_max = wat_max1 - (wat_max1 - wat_max2) * 0.2
		}
				
		x = dmt_eq(sub_str[d] * 1000000, wat_max * storm_f)
		gfit = optim(c(10, 10), gfcn)		
	
		x_store_abs = c(-3)
		n_store_abs = c(1)
		g_store_abs = c()
		f_store_abs = c()

		for (i in 1:50) {
			x_abs = x_store_abs[i]
			n_abs = n_store_abs[i]
	
			allo_abs = energy(x_abs, d, kappa)
			
			n_store_abs = c(n_store_abs, n_abs * sx(x_abs))
			g_store_abs = c(g_store_abs, allo_abs[4])
			x_store_abs = c(x_store_abs, allo_abs[5])
			f_store_abs = c(f_store_abs, allo_abs[6])	

		}
		f_mat[d,] = f_store_abs
		g_mat[d,] = g_store_abs
		x_mat[d,] = x_store_abs
		n_mat[d,] = n_store_abs

	}
	keep_lro_abs = rbind(keep_lro_abs, c((apply(f_mat * n_mat[,1:50], 1, sum))))
	
}

# ------------------------------------------------------------------------------- #

par(mfrow = c(1, 1), mar = c(5, 4, 3, 3))
plot(1, 1, ylim = c(-50, 300), xlim = c(0, 120), type = "n", xlab = "Distance from crest, m", ylab = "Lifetime reproductive output (1000s eggs)", axes = F) #, yaxs = "i")

for (i in 1:5) {
	lines(d_reef, (keep_lro_abs[i,])/1000, lty = i, col = "black")
}

axis(1, las = 2)
axis(2, at = seq(0, 300, 100), las = 2)
abline(v = c(0, 60, 120), col = c("red", "orange", "black"))

dat = read.delim("data/lizard_csf.txt")
dat = dat[dat$species == "hya",]
temp = c()
temp_95 = c()
for (i in seq(0, 70, 10)) {	
	temp = c(temp, sum(dat$area_cm_2[dat$dist >= i & dat$dist < (i + 10)])/(8*20*10000))
	#temp = c(temp, mean((dat$area_cm_2[dat$dist >= i & dat$dist < (i + 10)])))
	#temp_95 = c(temp_95, sd(dat$area_cm_2[dat$dist >= i & dat$dist < (i + 10)]))
}
temp = c(temp, rep(0, 4))
#temp_95 = c(temp_95, rep(0, 4))
rect(seq(1, 111, 10), rep(-50, 12), seq(9, 119, 10), -50+500*temp, border = NA, col = "grey")
#segments(seq(5, 115, 10), -5+5*(temp+temp_95)/2500, seq(5, 115, 10), -5+5*(temp-temp_95)/2500)
axis(4, at = c(-50, 75, 200), labels = c(0, 0.25, 0.5), las = 2, col = "grey")

# ------------------------------------------------------------------------------- #
# -------------------------------- SCENARIOS ------------------------------------ #
# ------------------------------------------------------------------------------- #

bg_s = 0.9
storm_f = 1
I_max = 1800
d_reef = c(0, 4, 8, 12, 16, 20, 40, 60, 80, 100, 120)
sub_str = c(rep(ss0, 7), ss0 - (ss0 - ss80)/2, ss80, ss80 - (ss80 - ss120)/2, ss120)
#sub_str = rep(0.3, length(d_reef))
kappa = 0.5

keep_lro_abs = c()

for (scn in c(1:4)) {

	storm_f = 1
	I_max = 1800
	sub_str = c(rep(ss0, 7), ss0 - (ss0 - ss80)/2, ss80, ss80 - (ss80 - ss120)/2, ss120)

	if (scn == 2) {
		sub_str = rep(ss0, length(d_reef))
	}
		
	f_mat = matrix(NA, length(d_reef), 50)
	g_mat = matrix(NA, length(d_reef), 50)
	x_mat = matrix(NA, length(d_reef), 51)
	n_mat = matrix(NA, length(d_reef), 51)
	
	for (d in 1:length(d_reef)) {
		if (d_reef[d] != 120) {
			wat_max = tapply(vel_dat$u[vel_dat$dist == d_reef[d]], list(vel_dat$year[vel_dat$dist == d_reef[d]]), max)
		} else {
			wat_max1 = tapply(vel_dat$u[vel_dat$dist == 100], list(vel_dat$year[vel_dat$dist == 100]), max)
			wat_max2 = tapply(vel_dat$u[vel_dat$dist == 200], list(vel_dat$year[vel_dat$dist == 200]), max)
			wat_max = wat_max1 - (wat_max1 - wat_max2) * 0.2
		}

		if (scn == 3) {
			#wat_max1 = tapply(vel_dat$u[vel_dat$dist == 100], list(vel_dat$year[vel_dat$dist == 100]), max)
			#wat_max2 = tapply(vel_dat$u[vel_dat$dist == 200], list(vel_dat$year[vel_dat$dist == 200]), max)
			#wat_max = wat_max1 - (wat_max1 - wat_max2) * 0.2
			wat_max = tapply(vel_dat$u[vel_dat$dist == d_reef[1]], list(vel_dat$year[vel_dat$dist == d_reef[1]]), max)
		}
				
		x = dmt_eq(sub_str[d] * 1000000, wat_max * storm_f)
		gfit = optim(c(10, 10), gfcn)		
	
		x_store_abs = c(-3)
		n_store_abs = c(1)
		g_store_abs = c()
		f_store_abs = c()

		for (i in 1:50) {
			x_abs = x_store_abs[i]
			n_abs = n_store_abs[i]
	
			allo_abs = energy(x_abs, d, kappa)
			if (scn == 2 | scn == 4) {
				allo_abs = energy(x_abs, 1, kappa)				
			}
			
			n_store_abs = c(n_store_abs, n_abs * sx(x_abs))
			g_store_abs = c(g_store_abs, allo_abs[4])
			x_store_abs = c(x_store_abs, allo_abs[5])
			f_store_abs = c(f_store_abs, allo_abs[6])	

		}
		f_mat[d,] = f_store_abs
		g_mat[d,] = g_store_abs
		x_mat[d,] = x_store_abs
		n_mat[d,] = n_store_abs

	}
	keep_lro_abs = rbind(keep_lro_abs, c((apply(f_mat * n_mat[,1:50], 1, sum))))
	
}

# ------------------------------------------------------------------------------- #

par(mfrow = c(1, 1), mar = c(5, 4, 3, 3))
plot(1, 1, ylim = c(1, 3), xlim = c(0, 120), type = "n", xlab = "Distance from crest, m", ylab = "Lifetime reproductive output (1000s eggs)", axes = F) #, yaxs = "i")

lines(d_reef, log10(keep_lro_abs[1,]/1000), lty = 1)
lines(d_reef, log10(keep_lro_abs[2,]/1000), lty = 2)
lines(d_reef, log10(keep_lro_abs[3,]/1000), lty = 3)
lines(d_reef, log10(keep_lro_abs[4,]/1000), lty = 4)

axis(1, las = 2)
axis(2, at = c(log10(seq(20, 90, 10)), log10(seq(200, 900, 100))), labels = NA, las = 2, tck = -0.03)
axis(2, at = c(1, 2, 3), labels = c(10, 100, 1000), las = 2, tick = T)
abline(v = c(0, 60, 120), col = c("red", "orange", "black"))

#dat = read.delim("lizard_csf.txt")
#dat = dat[dat$species == "hya",]
#temp = c()
#temp_95 = c()
#for (i in seq(0, 70, 10)) {	
#	temp = c(temp, sum(dat$area_cm_2[dat$dist >= i & dat$dist < (i + 10)])/(8*20*10000))
#	#temp = c(temp, mean((dat$area_cm_2[dat$dist >= i & dat$dist < (i + 10)])))
#	#temp_95 = c(temp_95, sd(dat$area_cm_2[dat$dist >= i & dat$dist < (i + 10)]))
#}
#temp = c(temp, rep(0, 4))
##temp_95 = c(temp_95, rep(0, 4))
#rect(seq(1, 111, 10), rep(-5, 12), seq(9, 119, 10), -5+50*temp, border = NA, col = "grey")
##segments(seq(5, 115, 10), -5+5*(temp+temp_95)/2500, seq(5, 115, 10), -5+5*(temp-temp_95)/2500)
#axis(4, at = c(-5, 1.25, 7.5), labels = c(0, 0.25, 0.5), las = 2, col = "grey")


# ------------------------------------------------------------------------------- #
# ------------------------------ MEAN COLONY AGE -------------------------------- #
# ------------------------------------------------------------------------------- #

#n_mat_b = n_mat >= 0.01
##temp = unlist(lapply(apply(n_mat_b == 1, 1, which), max))
#
#temp = c()
#for (i in 1:length(d_reef)) {
#	a = n_mat[i, max(which(n_mat[i,] > 0.5))]
#	b = n_mat[i, max(which(n_mat[i,] > 0.5)) + 1]
#	
#	temp = c(temp, max(which(n_mat[i,] > 0.5)) + 1 - (0.5 - b) / (a - b))
#}
#
#plot(d_reef, temp, xlim = c(0, 100), ylim = c(0, 8), type = "l", xlab = "Distance from crest, m", ylab = "Mean colony age, y")

# ------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------- #

par(mfrow = c(3, 1))
t = 25

#par(mar = c(5, 4, 2, 2))
plot(seq(0, t, 1), n_mat[1,1:(t+1)], type = "n", ylim = c(-0.2, 1), xlab = "Years", ylab = "Proportion survivors", col = "red", las = 2, axes = F)
abline(h = 0)
lines(seq(0, t, 1), n_mat[1,1:(t+1)], lty = 1, col = "red")
lines(seq(0, t, 1), n_mat[8,1:(t+1)], lty = 1, col = "orange")
lines(seq(0, t, 1), n_mat[11,1:(t+1)], lty = 1, col = "black")

axis(1, las = 2)
axis(2, at = seq(0, 1, 0.5), labels = seq(0, 1, 0.5), las = 2)

temp = c()
for (i in 1:length(d_reef)) {
	a = n_mat[i, max(which(n_mat[i,] > 0.5))]
	b = n_mat[i, max(which(n_mat[i,] > 0.5)) + 1]
	mn = max(which(n_mat[i,] > 0.5)) + 1 - (0.5 - b) / (a - b)

	c = x_mat[i, max(which(n_mat[i,] > 0.5))]
	d = x_mat[i, max(which(n_mat[i,] > 0.5)) + 1]
	mn_sz = x_mat[i, max(which(n_mat[i,] > 0.5)) + 1] + (c - d) * (0.5 - b) / (a - b)

	a = n_mat[i, max(which(n_mat[i,] > 0.975))]
	b = n_mat[i, max(which(n_mat[i,] > 0.975)) + 1]
	lo = max(which(n_mat[i,] > 0.975)) + 1 - (0.975 - b) / (a - b)

	c = x_mat[i, max(which(n_mat[i,] > 0.975))]
	d = x_mat[i, max(which(n_mat[i,] > 0.975)) + 1]
	lo_sz = x_mat[i, max(which(n_mat[i,] > 0.975)) + 1] + (c - d) * (0.975 - b) / (a - b)

	a = n_mat[i, max(which(n_mat[i,] > 0.025))]
	b = n_mat[i, max(which(n_mat[i,] > 0.025)) + 1]
	hi = max(which(n_mat[i,] > 0.025)) + 1 - (0.025 - b) / (a - b)

	c = x_mat[i, max(which(n_mat[i,] > 0.025))]
	d = x_mat[i, max(which(n_mat[i,] > 0.025)) + 1]
	hi_sz = x_mat[i, max(which(n_mat[i,] > 0.025)) + 1] + (c - d) * (0.025 - b) / (a - b)

	temp = rbind(temp, cbind(dist = i, mn, hi, lo, mn_sz, hi_sz, lo_sz))

}

points(temp[c(1, 8, 11),2], c(-0.15, -0.1, -0.05), col = c("red", "orange", "black"), pch = 20)
segments(temp[c(1, 8, 11),3], c(-0.15, -0.1, -0.05), temp[c(1, 8, 11),4], c(-0.15, -0.1, -0.05), col = c("red", "orange", "black"))

###

#plot(seq(1, 50, 1), g_mat[1,1:50], ylim = c(0, 0.15), type = "l", xlab = "time", ylab = "Proportion growth", col = "red", las = 2)
#lines(seq(1, 50, 1), g_mat[8,1:50], lty = 1, col = "green")
#lines(seq(1, 50, 1), g_mat[11,1:50], lty = 1, col = "blue")

#plot(seq(1, 50, 1), x_mat[1,1:50], type = "l", xlab = "time", ylab = "Colony size, m^2", col = "red", las = 2)
#lines(seq(1, 50, 1), x_mat[8,1:50], lty = 1, col = "green")
#lines(seq(1, 50, 1), x_mat[11,1:50], lty = 1, col = "blue")

#par(mar = c(5, 4, 2, 2))
#plot(seq(0, t, 1), (n_mat[1,1:(t+1)] * 10^x_mat[1,1:(t+1)]), type = "l", xlab = "Years", ylab = "Mean size of survivors, m^2", col = "red", las = 2, axes = F, ylim = c(0, 0.01))
#lines(seq(0, t, 1), (n_mat[8,1:(t+1)] * 10^x_mat[8,1:(t+1)]), lty = 1, col = "orange")
#lines(seq(0, t, 1), (n_mat[11,1:(t+1)] * 10^x_mat[11,1:(t+1)]), lty = 1, col = "black")
#axis(1, las = 2)
#axis(2, at = seq(0, 0.01, 0.005), las = 2)

plot(seq(0, t, 1), x_mat[1,1:(t+1)], type = "l", xlab = "Years", ylab = "Mean size of survivors, m^2", col = "red", las = 2, axes = F, ylim = c(-3, 1))
lines(seq(0, t, 1), x_mat[8,1:(t+1)], lty = 1, col = "orange")
lines(seq(0, t, 1), x_mat[11,1:(t+1)], lty = 1, col = "black")
axis(1, las = 2)
axis(2, at = seq(-3, 1, 1), labels = 10^seq(-3, 1, 1), las = 2)

#par(mar = c(5, 4, 2, 2))
plot(seq(0, t, 1), (f_mat[1,1:(t+1)] * n_mat[1,1:(t+1)])/1000, type = "l", xlab = "Years", ylab = "Mean fecundity of cohort, mm^3", col = "red", las = 2, axes = F, ylim = c(0, 55))
lines(seq(0, t, 1), (f_mat[8,1:(t+1)] * n_mat[8,1:(t+1)])/1000, lty = 1, col = "orange")
lines(seq(0, t, 1), (f_mat[11,1:(t+1)] * n_mat[11,1:(t+1)])/1000, lty = 1, col = "black")
axis(1, las = 2)
axis(2, at = seq(0, 50, 25), las = 2)	


#dat = read.delim("lizard_csf.txt")
#dat = dat[dat$species == "hya",]
#plot(dat$dist, log10(dat$area_cm_2/10000), pch = 20, cex = 0.5, col = "grey", xlim = c(0, 120), axes = F, xlab = "Distance from crest, m", ylab = "Log10 size, m^2")
##lines(d_reef, temp[,5], lty = 1)
##lines(d_reef, temp[,6], lty = 2)
##lines(d_reef, temp[,7], lty = 2)
##axis(1, las = 2)
##axis(2, at = seq(-3, 0, 1), labels = c(0.001, 0.01, 0.1, 1), las = 2)
##
#temp = c()
#for (i in seq(0, 70, 10)) {	
#	temp = c(temp, median(dat$area_cm_2[dat$dist >= i & dat$dist < (i + 10)]))
#}
#lines(seq(5, 75, 10), log10(temp/10000), col = "red")


#plot(seq(1, 50, 1), (e_mat[1,1:50]/10^x_mat[1,1:50]), type = "l", ylim = c(0, 1.3), xlab = "time", ylab = "Energy, %O2 m^-2", col = "red", las = 2)
#lines(seq(1, 50, 1), (e_mat[8,1:50]/10^x_mat[8,1:50]), lty = 1, col = "green")
#lines(seq(1, 50, 1), (e_mat[11,1:50]/10^x_mat[11,1:50]), lty = 1, col = "blue")

par(mar = c(5, 4, 2, 2))
plot(1, 1, ylim = c(3, 5), xlim = c(0, 120), type = "n", xlab = "Distance from crest, m", ylab = "Lifetime repro output, mm^3", axes = F) #, yaxs = "i")

lines(d_reef, log10(keep[1,]), lty = 2)
#lines(d_reef, log10(keep[2,]), lty = 3)
lines(d_reef, log10(keep[3,]), lty = 3)
lines(d_reef, log10(keep[4,]), lty = 1)

axis(1, las = 2)
axis(2, at = c(3, 4, 5), labels = expression(10^3, 10^4, 10^5), las = 2)
abline(v = c(0, 60, 120), col = c("red", "orange", "black"))

# ------------------------------------------------------------------------ #
# ------------------------------------------------------------------------ #
# ------------------------------------------------------------------------ #

bg_s = 0.9
storm_f = 1
d_reef = c(0, 4, 8, 12, 16, 20, 40, 60, 80, 100, 120)
sub_str = c(rep(ss0, 7), ss0 - (ss0 - ss80)/2, ss80, ss80 - (ss80 - ss120)/2, ss120)
#sub_str = rep(0.3, length(d_reef))

I_max = 1800

size_store = c()
fecu_store = c()
grow_store = c()
r_store = c()
ener_store = c()
d_grow_store = c()

for (d in 1:length(d_reef)) {
	
	n = 1
	size_temp = c()
	fecu_temp = c()
	grow_temp = c()
	r_temp = c()
	ener_temp = c()
	
	for (x in seq(-3, 1, 0.01)) {

		#W = 2 * sqrt((10^x)/pi) * 100
		#u = wat_mat[d,1] * 100
		#E = integrate(P_day, 0, 24)$value * 10^x
		#
		#u = wat_mat[1,1] * 100
		#E_max = integrate(P_day, 0, 24)$value * 10^x
			
		temp = energy(x, d, kappa = 0.5)
		
		#if (d == 1 & x == -2.5) {
		#	EE = E_max
		#}
		
		size_temp = c(size_temp, x)
		fecu_temp = c(fecu_temp, temp[6])
		grow_temp = c(grow_temp, temp[5])
		r_temp = c(r_temp, temp[3])
		ener_temp = c(ener_temp, temp[1])
	}
	
	size_store = rbind(size_store, size_temp)
	fecu_store = rbind(fecu_store, fecu_temp)
	grow_store = rbind(grow_store, grow_temp)
	r_store = rbind(r_store, r_temp)
	d_grow_store = rbind(d_grow_store, log10(10^grow_temp - 10^seq(-2.5, 1, 0.1)))
	ener_store = rbind(ener_store, ener_temp)

}

par(mfrow = c(2, 2), mar = c(5, 4, 2, 2))

#### SURV

plot(0, 0, axes = F, type = "n", xlim = c(0, 120), ylim = c(0, 1), xlab = "Distance from crest, m", ylab = "Survivorship")
lines(d_reef, s_store[,which(seq(-2.5, 1, 0.1) == -2)], lty = 4)
lines(d_reef, s_store[,which(seq(-2.5, 1, 0.1) == -1)], lty = 3)
lines(d_reef, s_store[,which(seq(-2.5, 1, 0.1) == 0)], lty = 2)
lines(d_reef, s_store[,which(seq(-2.5, 1, 0.1) == 1)], lty = 1)

axis(1, las = 2)
axis(2, at = c(0, 0.5, 1), las = 2)
abline(v = c(0, 60, 120), col = c("red", "orange", "black"))

#plot(seq(-2.5, 1, 0.1), s_store[which(d_reef == 0),], type = "l", axes = F, xlim = c(-2.5, 1), ylim = c(0,1), xlab = "Colony size, m^2", ylab = "Survival probability", col = "red")
#lines(seq(-2.5, 1, 0.1), s_store[which(d_reef == 60),], col = "orange")
#lines(seq(-2.5, 1, 0.1), s_store[which(d_reef == 120),], col = "black")
#axis(1, at = seq(-2, 1, 1), labels = c(0.01, 0.1, 1, 10), las = 2)
#axis(2, at = seq(0, 1, 0.5), las = 2)

### ENERGY

#contour(x = d_reef, y = seq(-2.5, 1, 0.1), z = (ener_store/10^size_store), axes = F, xlab = "Distance from crest, m", ylab = "Colony size, m^2", main = "Energy surplus per area, %O2.cm^-2.d^-1", cex.main = 0.7)

plot(0, 0, axes = F, type = "n", xlim = c(0, 120), ylim = c(1.5, 5), xlab = "Distance from crest, m", ylab = "Log10 Energy surplus, kJyr-1")
lines(d_reef, log10(ener_store[,which(seq(-2.5, 1, 0.1) == -2)]), lty = 4)
lines(d_reef, log10(ener_store[,which(seq(-2.5, 1, 0.1) == -1)]), lty = 3)
lines(d_reef, log10(ener_store[,which(seq(-2.5, 1, 0.1) == 0)]), lty = 2)
lines(d_reef, log10(ener_store[,which(seq(-2.5, 1, 0.1) == 1)]), lty = 1)

axis(1, las = 2)
axis(2, at = seq(2, 5, 1), labels = expression(10^2, 10^3, 10^4, 10^5), las = 2)
abline(v = c(0, 60, 120), col = c("red", "orange", "black"))
tt = seq(-2.5, 1, 0.001)[as.vector(unlist(lapply(apply(ener_store < E_growth, 1, which), min)))]

#plot(size_store[1,], (ener_store[1,]/10^size_store[1,]), xlim = c(-2.5, 1), ylim = c(0, 3), xlab = "Colony size, m^2", ylab = "Energy surplus per area, %O2.cm^-2.d^-1", axes = F, type = "l", col = "red")
#axis(1, at = seq(-2, 1, 1), labels = c(0.01, 0.1, 1, 10), las = 2)
#axis(2, at = seq(0, 3, 1), las = 2)
#lines(size_store[8,], (ener_store[8,]/10^size_store[8,]), lty = 1, col = "orange")
#lines(size_store[11,], (ener_store[11,]/10^size_store[11,]), lty = 1, col = "black")

#### fecu

plot(0, 0, axes = F, type = "n", xlim = c(0, 120), ylim = c(3.5, 7.7), xlab = "Distance from crest, m", ylab = "Log10 fecundity, no. eggs")
lines(d_reef, log10(fecu_store[,which(seq(-2.5, 1, 0.1) == -2)]), lty = 4)
lines(d_reef, log10(fecu_store[,which(seq(-2.5, 1, 0.1) == -1)]), lty = 3)
lines(d_reef, log10(fecu_store[,which(seq(-2.5, 1, 0.1) == 0)]), lty = 2)
lines(d_reef, log10(fecu_store[,which(seq(-2.5, 1, 0.1) == 1)]), lty = 1)

axis(1, las = 2)
axis(2, at = seq(4, 7, 1), labels = expression(10^4, 10^5, 10^6, 10^7), las = 2)
abline(v = c(0, 60, 120), col = c("red", "orange", "black"))

#
#contour(x = d_reef, y = seq(-2.5, 1, 0.1), z = log10(fecu_store),axes = F, xlab = "Distance from crest, m", ylab = "Colony size, m^2", main = "Log10 fecundity (no. eggs)", cex.main = 0.7)
##contour(x = d_reef, y = seq(-2.5, 1, 0.1), z = log10(fecu_store/10^size_store), levels = seq(4, 5.4, 0.2),axes = F, xlab = "Distance from crest, m", ylab = "Colony size, m^2", main = "Log10 fecundity per area, mm^3m^-2", cex.main = 0.7)
#axis(1, las = 2)
#axis(2, at = seq(-2, 1, 1), labels = c(0.01, 0.1, 1, 10), las = 2)
#abline(v = c(0, 60, 120), col = c("red", "orange", "black"))

#plot(size_store[1,], log10(fecu_store[1,]/10^size_store[1,]), xlim = c(-2.5, 1), ylim = c(4, 6), xlab = "Colony size, m^2", ylab = "Log10 fecundity per area, mm^3m^-2", axes = F, type = "l", col = "red")
#axis(1, at = seq(-2, 1, 1), labels = c(0.01, 0.1, 1, 10), las = 2)
#axis(2, at = seq(4, 6, 1), labels = expression(10^4, 10^5, 10^6), las = 2)
#lines(size_store[8,], log10(fecu_store[8,]/10^size_store[8,]), lty = 1, col = "orange")
#lines(size_store[11,], log10(fecu_store[11,]/10^size_store[11,]), lty = 1, col = "black")

#### grow

plot(0, 0, axes = F, type = "n", xlim = c(0, 120), ylim = c(-2.5, 0), xlab = "Distance from crest, m", ylab = "Log10 added radius, m")
lines(d_reef, log10(r_store[,which(seq(-2.5, 1, 0.1) == -2)]), lty = 4)
lines(d_reef, log10(r_store[,which(seq(-2.5, 1, 0.1) == -1)]), lty =3)
lines(d_reef, log10(r_store[,which(seq(-2.5, 1, 0.1) == 0)]), lty = 2)
lines(d_reef, log10(r_store[,which(seq(-2.5, 1, 0.1) == 1)]), lty = 1)

axis(1, las = 2)
axis(2, at = seq(-2, 0, 1), labels = expression(0.01, 0.1, 1), las = 2)
abline(v = c(0, 60, 120), col = c("red", "orange", "black"))


#contour(x = d_reef, y = seq(-2.5, 1, 0.1), z = (d_grow_store), levels = log10(c(0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1)), labels = c(0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1), axes = F, xlab = "Distance from crest, m", ylab = "Colony size, m^2", main = "Growth (added area), m^2", cex.main = 0.7)
#axis(1, las = 2)
#axis(2, at = seq(-2, 1, 1), labels = c(0.01, 0.1, 1, 10), las = 2)
#abline(v = c(0, 60, 120), col = c("red", "orange", "black"))

#plot(size_store[1,], grow_store[1,], ylim = c(-2.5, 1), xlim = c(-2.5, 1), xlab = "Colony size t, m^2", ylab = "Colony size t+1, m^2", axes = F, type = "l", col = "red")
#axis(1, at = seq(-2, 1, 1), labels = c(0.01, 0.1, 1, 10), las = 2)
#axis(2, at = seq(-2, 1, 1), labels = c(0.01, 0.1, 1, 10), las = 2)
#abline(0, 1, lty = 2)
#lines(size_store[8,], grow_store[8,], lty = 1, col = "orange")
#lines(size_store[11,], grow_store[11,], lty = 1, col = "black")

####
####

par(mfcol = c(2, 2), mar = c(5, 4, 2, 2))

plot(0, 0, axes = F, type = "n", xlim = c(-3, 1), ylim = c(-3, 0), xlab = "Colony size, m^2", ylab = "Added radius, m")
lines(seq(-3, 1, 0.01), log10(r_store[which(d_reef == 0),]), col = "red")
lines(seq(-3, 1, 0.01), log10(r_store[which(d_reef == 60),]), col = "orange")
lines(seq(-3, 1, 0.01), log10(r_store[which(d_reef == 120),]), col = "black")
abline(v = log10(pi*(7/200)^2), col = "grey")
abline(h = log10(0.1/2), lty=2)

axis(1, at = c(log10(seq(0.002, 0.009, 0.001)), log10(seq(0.02, 0.09, 0.01)), log10(seq(0.2, 0.9, 0.1)), log10(seq(2, 9, 1))), labels = NA, las = 2, tck = -0.03)
axis(1, at = c(-3, -2, -1, 0, 1), labels = c(0.001, 0.01, 0.1, 1, 10), las = 2, tick = T)
axis(2, at = seq(-3, 0, 1), labels = expression(0.001, 0.01, 0.1, 1), las = 2)

####

plot(0, 0, axes = F, type = "n", xlim = c(-3, 1), ylim = c(3, 7), xlab = "Colony size, m^2", ylab = "Log10 fecundity, x10^3 eggs")
lines(seq(-3, 1, 0.01), log10(fecu_store[which(d_reef == 0),]), col = "red")
lines(seq(-3, 1, 0.01), log10(fecu_store[which(d_reef == 60),]), col = "orange")
lines(seq(-3, 1, 0.01), log10(fecu_store[which(d_reef == 120),]), col = "black")
abline(v = log10(pi*(7/200)^2), col = "grey")
lines(seq(log10(pi*(7/200)^2), 1, 0.01), log10(hall(seq(log10(pi*(7/200)^2), 1, 0.01))), lty = 2)

#axis(1, at = seq(-2, 1, 1), labels = expression(0.01, 0.1, 1, 10), las = 2)
axis(1, at = c(log10(seq(0.002, 0.009, 0.001)), log10(seq(0.02, 0.09, 0.01)), log10(seq(0.2, 0.9, 0.1)), log10(seq(2, 9, 1))), labels = NA, las = 2, tck = -0.03)
axis(1, at = c(-3, -2, -1, 0, 1), labels = c(0.001, 0.01, 0.1, 1, 10), las = 2, tick = T)

axis(2, at = seq(3, 7, 1), labels = expression(10^3, 10^4, 10^5, 10^6, 10^7), las = 2)


######


bg_s = 0.9
storm_f = 1
d_reef = c(0, 4, 8, 12, 16, 20, 40, 60, 80, 100, 120)
sub_str = c(rep(ss0, 7), ss0 - (ss0 - ss80)/2, ss80, ss80 - (ss80 - ss120)/2, ss120)
#sub_str = rep(0.3, length(d_reef))

I_max = 1800

size_store = c()
ffecu_store = c()
grow_store = c()
rr_store = c()
ener_store = c()
d_grow_store = c()

for (kappa in c(0.05, 0.25, 0.5, 0.75, 0.95)) {
	for (d in c(1)) {
	
		n = 1
		size_temp = c()
		fecu_temp = c()
		grow_temp = c()
		r_temp = c()
		ener_temp = c()
	
		for (x in seq(-3, 1, 0.01)) {

			#W = 2 * sqrt((10^x)/pi) * 100
			#u = wat_mat[d,1] * 100
			#E = integrate(P_day, 0, 24)$value * 10^x
			#
			#u = wat_mat[1,1] * 100
			#E_max = integrate(P_day, 0, 24)$value * 10^x
			
			temp = energy(x, d, kappa = kappa)
		
			#if (d == 1 & x == -2.5) {
			#	EE = E_max
			#}
		
			size_temp = c(size_temp, x)
			fecu_temp = c(fecu_temp, temp[6])
			grow_temp = c(grow_temp, temp[5])
			r_temp = c(r_temp, temp[3])
			ener_temp = c(ener_temp, temp[1])
		}
	
		size_store = rbind(size_store, size_temp)
		ffecu_store = rbind(ffecu_store, fecu_temp)
		grow_store = rbind(grow_store, grow_temp)
		rr_store = rbind(rr_store, r_temp)
		#d_grow_store = rbind(d_grow_store, log10(10^grow_temp - 10^seq(-2.5, 1, 0.1)))
		#ener_store = rbind(ener_store, ener_temp)

	}
}

####
####
par(mfcol = c(2, 1), mar = c(5, 4, 2, 2))

plot(0, 0, axes = F, type = "n", xlim = c(-3, 1), ylim = c(-0.5, 0.5), xlab = "Colony size, m^2", ylab = "Added radius, m")
lines(seq(-3, 1, 0.01), (rr_store[1,]), col = "black")
lines(seq(-3, 1, 0.01), (rr_store[2,]), col = "black")
lines(seq(-3, 1, 0.01), (rr_store[3,]), col = "black")
lines(seq(-3, 1, 0.01), (rr_store[4,]), col = "black")
lines(seq(-3, 1, 0.01), (rr_store[5,]), col = "black")
abline(v = log10(pi*(7/200)^2), col = "grey")
abline(h = (0.1/2), lty=2)

axis(1, at = c(log10(seq(0.002, 0.009, 0.001)), log10(seq(0.02, 0.09, 0.01)), log10(seq(0.2, 0.9, 0.1)), log10(seq(2, 9, 1))), labels = NA, las = 2, tck = -0.03)
axis(1, at = c(-3, -2, -1, 0, 1), labels = c(0.001, 0.01, 0.1, 1, 10), las = 2, tick = T)
axis(2)
#axis(2, at = seq(-3, 0, 1), labels = expression(0.001, 0.01, 0.1, 1), las = 2)

####

plot(0, 0, axes = F, type = "n", xlim = c(-3, 1), ylim = c(3, 7), xlab = "Colony size, m^2", ylab = "Log10 fecundity, x10^3 eggs")
#lines(seq(-2.5, 1, 0.01), log10(fecu_store[which(d_reef == 0),]), col = "red")
#lines(seq(-2.5, 1, 0.01), log10(fecu_store[which(d_reef == 60),]), col = "orange")
#lines(seq(-2.5, 1, 0.01), log10(fecu_store[which(d_reef == 120),]), col = "black")
lines(seq(-3, 1, 0.01), log10(ffecu_store[1,]), col = "black")
lines(seq(-3, 1, 0.01), log10(ffecu_store[2,]), col = "black")
lines(seq(-3, 1, 0.01), log10(ffecu_store[3,]), col = "black")
lines(seq(-3, 1, 0.01), log10(ffecu_store[4,]), col = "black")
lines(seq(-3, 1, 0.01), log10(ffecu_store[5,]), col = "black")
abline(v = log10(pi*(7/200)^2), col = "grey")
lines(seq(log10(pi*(7/200)^2), 1, 0.01), log10(hall(seq(log10(pi*(7/200)^2), 1, 0.01))), lty = 2)
#lines(seq(log10(pi*(7/200)^2), 1, 0.01), log10(hall(seq(log10(pi*(7/200)^2), 1, 0.01), 1.03 + 0.15, 1.28 + 0.05)), lty = 2)
#lines(seq(log10(pi*(7/200)^2), 1, 0.01), log10(hall(seq(log10(pi*(7/200)^2), 1, 0.01), 1.03 + 0.15, 1.28 - 0.05)), lty = 2)

#axis(1, at = seq(-2, 1, 1), labels = expression(0.01, 0.1, 1, 10), las = 2)
axis(1, at = c(log10(seq(0.002, 0.009, 0.001)), log10(seq(0.02, 0.09, 0.01)), log10(seq(0.2, 0.9, 0.1)), log10(seq(2, 9, 1))), labels = NA, las = 2, tck = -0.03)
axis(1, at = c(-3, -2, -1, 0, 1), labels = c(0.001, 0.01, 0.1, 1, 10), las = 2, tick = T)

axis(2, at = seq(3, 7, 1), labels = expression(10^3, 10^4, 10^5, 10^6, 10^7), las = 2)




########## SEAN ASK FOR ARC  T versus T+1

plot(0, 0, axes = F, type = "n", xlim = c(-3, 1), ylim = c(-3, 1), xlab = "Colony size t, m^2", ylab = "Colony size t+1, m^2")
lines(seq(-3, 1, 0.01), (grow_store[1,]), col = "black")
lines(seq(-3, 1, 0.01), (grow_store[2,]), col = "black")
lines(seq(-3, 1, 0.01), (grow_store[3,]), col = "black")
lines(seq(-3, 1, 0.01), (grow_store[4,]), col = "black")
lines(seq(-3, 1, 0.01), (grow_store[5,]), col = "black")
abline(v = log10(pi*(7/200)^2), col = "grey")
abline(0, 1, lty=2)

axis(1, at = c(log10(seq(0.002, 0.009, 0.001)), log10(seq(0.02, 0.09, 0.01)), log10(seq(0.2, 0.9, 0.1)), log10(seq(2, 9, 1))), labels = NA, las = 2, tck = -0.03)
axis(1, at = c(-3, -2, -1, 0, 1), labels = c(0.001, 0.01, 0.1, 1, 10), las = 2, tick = T)
axis(2, at = c(log10(seq(0.002, 0.009, 0.001)), log10(seq(0.02, 0.09, 0.01)), log10(seq(0.2, 0.9, 0.1)), log10(seq(2, 9, 1))), labels = NA, las = 2, tck = -0.03)
axis(2, at = c(-3, -2, -1, 0, 1), labels = c(0.001, 0.01, 0.1, 1, 10), las = 2, tick = T)
#axis(2, at = seq(-3, 0, 1), labels = expression(0.001, 0.01, 0.1, 1), las = 2)


k = c(0.05, 0.25, 0.5, 0.75, 0.95)
store = c()
for (i in 1:5) {
	
	store = rbind(store, cbind(rep(k[i], dim(grow_store)[2]), seq(-3, 1, 0.01), grow_store[i,]))
	
}

write.table(store, "output/for_sean.txt", col.names = F, row.names = F, quote = F, sep = "\t")






