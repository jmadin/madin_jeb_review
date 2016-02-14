# ---------------------- MADIN & CONNOLLY 2006 -------------------- #
library(VGAM)

rfcn <- function(p) {												# Normal Variable Standard Dev
	a1<-p[1]
	a2<-p[2]
	a3<-p[3]
	return(sum(-log(dnorm(y, a1 + a2*x, a3))))
}

gfcn <- function(p) {
	k1<-p[1]
	k2<-p[2]
	sum(-log(dgumbel(x, k1, k2)))
}

g_ir <- function(Ir) {												# from Bosscher & Meesters, ***
	tanh(Ir / 351.6)
}

dmt_eq <- function(sig, uw, pw = 1025) {								# From Madin & Connolly, 2006
	sig / (pw * uw^2)
}

survive <- function(csf, size) {
	mu  = vfit$par[1] + vfit$par[2] * (size)
	sig = vfit$par[3]
	
	return((1 - pgumbel(10^csf, gfit$par[1], gfit$par[2])) * dnorm(csf, mu, sig))
}

sx <- function(x) {
	if (length(x) > 1) {
		keep = c()
		for (i in 1:length(x)) {
			keep = c(keep, integrate(survive, -5, 5, size=x[i], stop.on.error = F)$value)
		}
	} else {
		keep = integrate(survive, -5, 5, size = x, stop.on.error = F)$value
	}
	return(keep * bg_s)
}

ss_opt <- function(f, M) {
	f <<- f
	Mo <<- M
	Mo[1,] <<- M[1,] + f * 10^ss$mids

	ev <<- eigen(Mo)$vectors[,1] / sum(eigen(Mo)$vectors[,1])
	return(as.real(sum((ev - ss$density)^2)))
}

fx <- function(x) {
	temp = 1.03 + 1.28 * log10(10^x * 10000)
	temp[x < log10(100/10000)] = NA
	return(temp)
}

# ----------------------------------------------------------------- #
# ---------------------------- DATA ------------------------------- #
# ----------------------------------------------------------------- #

wav <- vel_dat <- read.delim("output/data/wave_db_hxua.csv")
csf_dat <- read.delim("data/lizard_csf.txt")
gth_dat <- read.delim("data/hya_growth.txt")/10000 # m^2
gth_dat <- gth_dat[gth_dat[,2] > 0,][-33,]

pop_size <- csf_dat[csf_dat[,1] == "hya", 3]/10000
pop_size[pop_size < 10/10000] <- 10/10000

# substrate strength

sub <- read.delim("data/penetr_final.txt", as.is = T)

sub$distance_m <- sub$site
sub$distance_m[sub$distance_m == "crest"] <- 0
sub$distance_m[sub$distance_m == "back_crest"] <- 40
sub$distance_m[sub$distance_m == "flat"] <- 80
sub$distance_m[sub$distance_m == "back"] <- 120
sub$distance_m <- as.numeric(sub$distance_m)

sub <- sub[order(sub$distance_m),]
sub2 <- sub[sub$depth_mm <= 100,]

ss0 <- sub2$tens_MN[sub2$distance_m == 0]
ss40 <- sub2$tens_MN[sub2$distance_m == 40]
ss80 <- sub2$tens_MN[sub2$distance_m == 80]
ss120 <- sub2$tens_MN[sub2$distance_m == 120]

#par(mfrow = c(3, 1))


# ------------------------- WATER VELOCITY ---------------------------------

uu0 <- wav$u[wav$dist == 0]
uu4 <- wav$u[wav$dist == 4]
uu8 <- wav$u[wav$dist == 8]
uu12 <- wav$u[wav$dist == 12]
uu16 <- wav$u[wav$dist == 16]
uu20 <- wav$u[wav$dist == 20]
uu40 <- wav$u[wav$dist == 40]
uu60 <- wav$u[wav$dist == 60]
uu80 <- wav$u[wav$dist == 80]
uu120 <- wav$u[wav$dist == 100]

par(mar = c(2, 3, 1, 3))

plot(0, 0, ylab = NA, xlab = NA, ylim = c(1, 10), xlim = c(0, 120), axes = F, type = "n")

axis(1, c(0, 60, 120))
axis(2, 7 + seq(0, 3, 1), labels = seq(0, 3, 1), las = 2)

#axis(4, - 1 + log10(c(seq(10^4, 10^5, 10^4), seq(10^5, 10^6, 10^5), seq(10^6, 10^7, 10^6))), labels = F, tck = -0.02)
axis(4, seq(4, 7, 1), labels = expression(10^4, 10^5, 10^6, 10^7), las = 2)

#axis(2, 6 + log10(c(seq(1, 10, 1), seq(10, 100, 10), seq(200, 1000, 100))), labels = F, tck = -0.02)
axis(2, 1 + c(0, 1, 2, 3), labels = expression(10^0, 10^1, 10^2, 10^3), las = 2)


polygon(c(c(0, 4, 8, 12, 16, 20, 40, 60, 80, 120), rev(c(0, 4, 8, 12, 16, 20, 40, 60, 80, 120))), 7 + c(c(sort(uu0)[round(length(uu0) * 0.025)], sort(uu4)[round(length(uu4) * 0.025)], sort(uu8)[round(length(uu8) * 0.025)], sort(uu12)[round(length(uu12) * 0.025)], sort(uu16)[round(length(uu16) * 0.025)], sort(uu20)[round(length(uu20) * 0.025)], sort(uu40)[round(length(uu40) * 0.025)], sort(uu60)[round(length(uu60) * 0.025)], sort(uu80)[round(length(uu80) * 0.025)], sort(uu120)[round(length(uu120) * 0.025)]), rev(c(sort(uu0)[round(length(uu0) * 0.975)], sort(uu4)[round(length(uu4) * 0.975)], sort(uu8)[round(length(uu8) * 0.975)], sort(uu12)[round(length(uu12) * 0.975)], sort(uu16)[round(length(uu16) * 0.975)], sort(uu20)[round(length(uu20) * 0.975)], sort(uu40)[round(length(uu40) * 0.975)], sort(uu60)[round(length(uu60) * 0.975)], sort(uu80)[round(length(uu80) * 0.975)], sort(uu120)[round(length(uu120) * 0.975)]))), border = "grey", col = "grey")

lines(c(0, 4, 8, 12, 16, 20, 40, 60, 80, 120), 7 + (c(sort(uu0)[round(length(uu0) * 0.5)], sort(uu4)[round(length(uu4) * 0.5)], sort(uu8)[round(length(uu8) * 0.5)], sort(uu12)[round(length(uu12) * 0.5)], sort(uu16)[round(length(uu16) * 0.5)], sort(uu20)[round(length(uu20) * 0.5)], sort(uu40)[round(length(uu40) * 0.5)], sort(uu60)[round(length(uu60) * 0.5)], sort(uu80)[round(length(uu80) * 0.5)], sort(uu120)[round(length(uu120) * 0.5)])), lty = 2)

####

polygon(c(c(0, 40, 80, 120), rev(c(0, 40, 80, 120))), c(log10(10^6 * c(sort(ss0)[round(length(ss0) * 0.025)], sort(ss40)[round(length(ss40) * 0.025)], sort(ss80)[round(length(ss80) * 0.025)], sort(ss120)[round(length(ss120) * 0.025)])), rev(log10(10^6 * c(sort(ss0)[round(length(ss0) * 0.975)], sort(ss40)[round(length(ss40) * 0.975)], sort(ss80)[round(length(ss80) * 0.975)], sort(ss120)[round(length(ss120) * 0.975)])))), col = "grey", border = "grey")

lines(c(0, 40, 80, 120), log10(10^6 * c(sort(ss0)[round(length(ss0) * 0.5)], sort(ss40)[round(length(ss40) * 0.5)], sort(ss80)[round(length(ss80) * 0.5)], sort(ss120)[round(length(ss120) * 0.5)])), lwd = 1, lty = 2)


# SUBSTRATE MEDIAN
rep <- 10^6
off <- -3

dmt0 <- (sample(ss0, rep, replace = T) * 10^6) / (sample(uu0, rep, replace = T)^2 * 1025)
dmt40 <- (sample(ss40, rep, replace = T) * 10^6) / (sample(uu40, rep, replace = T)^2 * 1025)
dmt80 <- (sample(ss80, rep, replace = T) * 10^6) / (sample(uu80, rep, replace = T)^2 * 1025)
dmt120 <- (sample(ss120, rep, replace = T) * 10^6) / (sample(uu120, rep, replace = T)^2 * 1025)


polygon(c(c(0, 40, 80, 120), rev(c(0, 40, 80, 120))), 0 + c(log10(c(sort(dmt0)[round(length(dmt0) * 0.5)], sort(dmt40)[round(length(dmt40) * 0.5)], sort(dmt80)[round(length(dmt80) * 0.5)], sort(dmt120)[round(length(dmt120) * 0.5)])), rev(log10(c(sort(dmt0)[round(length(dmt0) * 0.025)], sort(dmt40)[round(length(dmt40) * 0.025)], sort(dmt80)[round(length(dmt80) * 0.025)], sort(dmt120)[round(length(dmt120) * 0.025)])))), col = "grey", border = "grey")

lines(c(0, 40, 80, 120), 0 + log10(c(sort(dmt0)[round(length(dmt0) * 0.5)], sort(dmt40)[round(length(dmt40) * 0.5)], sort(dmt80)[round(length(dmt80) * 0.5)], sort(dmt120)[round(length(dmt120) * 0.5)])), lty = 2)

#abline(v = 0)
#abline(v = 60, lty = 2)
#abline(v = 120, lty = 3)





# -------------------------

ss0 <- median(sub2$tens_MN[sub2$distance_m == 0])
ss40 <- median(sub2$tens_MN[sub2$distance_m == 40])
ss80 <- median(sub2$tens_MN[sub2$distance_m == 80])
ss120 <- median(sub2$tens_MN[sub2$distance_m == 120])

sub_str <- c(rep(ss0, 7), ss0 - (ss0 - ss80)/2, ss80, ss80 - (ss80 - ss120)/2, ss120)



# ----------------------------------------------------------------- #
# -------------------------- SURVIVAL ----------------------------- #
# ----------------------------------------------------------------- #

x <- log10(csf_dat[csf_dat[,1] == "hya", 3]/10000)
y <- log10(csf_dat[csf_dat[,1] == "hya", 4])
vfit <- optim(c(2, 1, 1), rfcn)           							# Fit hyacinthus CSF/Size data

d_reef <- c(0, 4, 8, 12, 16, 20, 40, 60, 80, 100, 120)

#quartz(width = 7, height = 3)


plot(d_reef, sub_str, type = "l", axes = F)


bg_s <- 0.9
storm_f <- 1
d_reef <- c(0, 4, 8, 12, 16, 20, 40, 60, 80, 100, 120)
#sub_str = rep(0.3, length(d_reef))
sub_str <- c(rep(ss0, 7), ss0 - (ss0 - ss80)/2, ss80, ss80 - (ss80 - ss120)/2, ss120)

s_store <- c()	
for (d in 1:length(d_reef)) {
	if (d_reef[d] != 120) {
		wat_max <- tapply(vel_dat$u[vel_dat$dist == d_reef[d]], list(vel_dat$year[vel_dat$dist == d_reef[d]]), max)
	} else {
		wat_max1 <- tapply(vel_dat$u[vel_dat$dist == 100], list(vel_dat$year[vel_dat$dist == 100]), max)
		wat_max2 <- tapply(vel_dat$u[vel_dat$dist == 200], list(vel_dat$year[vel_dat$dist == 200]), max)
		wat_max <- wat_max1 - (wat_max1 - wat_max2) * 0.2
	}
	x <- dmt_eq(sub_str[d] * 1000000, wat_max * storm_f)
	gfit <- optim(c(2.6, 1), gfcn)		
	
	s_store <- rbind(s_store, sx(seq(-2.5, 1, 0.1)))	
}


par(mfrow = c(1, 2))
image(x = d_reef, y = seq(-2.5, 1, 0.1), z = s_store, axes = F, xlab = "Distance from crest, m", ylab = "Colony size, m^2", col = rev(grey(1:10/10)))
axis(1, las = 2)
axis(2, at = seq(-2, 1, 1), labels = c(0.01, 0.1, 1, 10), las = 2)
abline(v = c(0, 60, 120), col = c("red", "orange", "green"))

plot(seq(-2.5, 1, 0.1), s_store[which(d_reef == 0),], type = "l", axes = F, xlim = c(-2.5, 1), ylim = c(0,1), xlab = "Colony size, m^2", ylab = "Survival probability", col = "red")
lines(seq(-2.5, 1, 0.1), s_store[which(d_reef == 60),], col = "orange")
lines(seq(-2.5, 1, 0.1), s_store[which(d_reef == 120),], col = "green")

axis(1, at = seq(-2, 1, 1), labels = c(0.01, 0.1, 1, 10), las = 2)
axis(2, at = seq(0, 1, 0.5), las = 2)
#abline(v = c(0, 60, 120), col = c("red", "orange", "green"))
#text(115, 1.2, "B")

















