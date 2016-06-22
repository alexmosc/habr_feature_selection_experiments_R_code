library(scatterplot3d)
library(rgl)
library(geometry)
library(gbm)




library(onion)
data(bunny)

#XYZ point plot
open3d()
points3d(bunny, col=8, size=0.1)

bunny_dat <- as.data.frame(bunny)

inputs <- append(as.numeric(bunny_dat$x),
		   as.numeric(bunny_dat$y))

for (i in 1:10){
	
	naming <- paste('input_noise_', i, sep = '')
	bunny_dat[, eval(naming)] <- inputs[sample(length(inputs), nrow(bunny_dat), replace = T)]
	
	
}


### greedy search with filter function

library(FSelector)

sampleA <- bunny_dat[sample(nrow(bunny_dat), nrow(bunny_dat) / 2, replace = F), c("x",
			    "y",
			    "input_noise_1",  
			    "input_noise_2",
			    "input_noise_3",
			    "input_noise_4",
			    "input_noise_5",
			    "input_noise_6", 
			    "input_noise_7",
			    "input_noise_8",
			    "input_noise_9",
			    "input_noise_10",
			    "z")]

sampleB <- bunny_dat[!row.names(bunny_dat) %in% rownames(sampleA), c("x",
				  "y",
				  "input_noise_1",  
				  "input_noise_2",
				  "input_noise_3",
				  "input_noise_4",
				  "input_noise_5",
				  "input_noise_6", 
				  "input_noise_7",
				  "input_noise_8",
				  "input_noise_9",
				  "input_noise_10",
				  "z")]

linear_fit <- function(subset){
	dat <- sampleA[, c(subset, "z")]
	
	lm_m <- lm(formula = z ~ x * y,
		    data = dat,
		    model = T)
	
	lm_predict <- predict(lm_m,
				 newdata = sampleB)
	
	r_sq_validate <- 1 - sum((sampleB$z - lm_predict) ^ 2) / sum((sampleB$z - mean(sampleB$z)) ^ 2)
	
	print(subset)
	print(r_sq_validate)
	return(r_sq_validate)
}

linear_fit_f <- function(subset){
	dat <- sampleA[, c(subset, "z")]
	
	lm_m <- lm(formula = z ~.,
		    data = dat,
		    model = T)
	
	print(subset)
	print(summary(lm_m)$fstatistic[[1]])
	return(summary(lm_m)$fstatistic[[1]])
}

subset <- backward.search(attributes = names(sampleA)[1:(ncol(sampleA) - 1)], eval.fun = linear_fit_f)
#subset <- forward.search(attributes = names(sampleA)[1:(ncol(sampleA) - 1)], eval.fun = linear_fit_f)


# single predictor importance based on correlation

correlation_arr <- data.frame()

for (i in 1:12){
	
	correlation_arr[i, 1] <- colnames(sampleA)[i]
	correlation_arr[i, 2] <- cor(sampleA[, i], sampleA[, 'z'])
	
}

plot(correlation_arr$V2, type = 's')


### importance of predictors using GBM

gbm_dat <- bunny_dat[, c("x",
		     "y",
		     "input_noise_1",  
		     "input_noise_2",
		     "input_noise_3",
		     "input_noise_4",
		     "input_noise_5",
		     "input_noise_6", 
		     "input_noise_7",
		     "input_noise_8",
		     "input_noise_9",
		     "input_noise_10",
		     "z")]

gbm_fit <- gbm(formula = z ~.,
		 distribution = "gaussian",
		 data = gbm_dat,
		 n.trees = 500,
		 interaction.depth = 12,
		 n.minobsinnode = 100,
		 shrinkage = 0.0001,
		 bag.fraction = 0.9,
		 train.fraction = 0.7,
		 n.cores = 6)

gbm.perf(object = gbm_fit, 
	  plot.it = TRUE, 
	  oobag.curve = F, 
	  overlay = TRUE)

summary(gbm_fit)





### simulated annealing search with linear model interactions stats

library(scales)
library(GenSA)

sampleA <- bunny_dat[sample(nrow(bunny_dat), nrow(bunny_dat) / 2, replace = F), c("x",
											     "y",
											     "input_noise_1",  
											     "input_noise_2",
											     "input_noise_3",
											     "input_noise_4",
											     "input_noise_5",
											     "input_noise_6", 
											     "input_noise_7",
											     "input_noise_8",
											     "input_noise_9",
											     "input_noise_10",
											     "z")]

sampleB <- bunny_dat[!row.names(bunny_dat) %in% rownames(sampleA), c("x",
									      "y",
									      "input_noise_1",  
									      "input_noise_2",
									      "input_noise_3",
									      "input_noise_4",
									      "input_noise_5",
									      "input_noise_6", 
									      "input_noise_7",
									      "input_noise_8",
									      "input_noise_9",
									      "input_noise_10",
									      "z")]

#calculate parameters
predictor_number <- dim(sampleA)[2] - 1
sample_size <- dim(sampleA)[1]
par_v <- runif(predictor_number, min = 0, max = 1)
par_low <- rep(0, times = predictor_number)
par_upp <- rep(1, times = predictor_number)


############### the fitness function
sa_fit_f <- function(par){
	
	indexes <- c(1:predictor_number)
	
	for (i in 1:predictor_number){
		if (par[i] >= threshold) {
			indexes[i] <- i
		} else {
			indexes[i] <- 0
		}
	}
	
	local_predictor_number <- 0
	for (i in 1:predictor_number){
		if (indexes[i] > 0) {
			local_predictor_number <- local_predictor_number + 1
		}
	}
	
	if (local_predictor_number > 0) {
		
		sampleAf <- as.data.frame(sampleA[, c(indexes[], dim(sampleA)[2])])
		
		lm_m <- lm(formula = z ~ .^2,
			    data = sampleAf,
			    model = T)
		
		lm_predict <- predict(lm_m,
					 newdata = sampleB)
		
		r_sq_validate <- 1 - sum((sampleB$z - lm_predict) ^ 2) / sum((sampleB$z - mean(sampleB$z)) ^ 2)

	} else  {
		r_sq_validate <- 0
	}

	return(-r_sq_validate)
}

# fittness based on p-values
sa_fit_f2 <- function(par){
	
	indexes <- c(1:predictor_number)
	
	for (i in 1:predictor_number){
		if (par[i] >= threshold) {
			indexes[i] <- i
		} else {
			indexes[i] <- 0
		}
	}
	
	local_predictor_number <- 0
	for (i in 1:predictor_number){
		if (indexes[i] > 0) {
			local_predictor_number <- local_predictor_number + 1
		}
	}
	
	if (local_predictor_number > 0) {
		
		sampleAf <- as.data.frame(sampleA[, c(indexes[], dim(sampleA)[2])])
		
		lm_m <- lm(formula = z ~ .^2,
			    data = sampleAf,
			    model = T)
		
		mean_val <- mean(summary(lm_m)$coefficients[, 4])
		
	} else  {
		mean_val <- 1
	}
	
	return(mean_val)
}


#stimating threshold for variable inclusion

threshold <- 0.5 # do it simply

#run feature selection

start <- Sys.time()

sao <- GenSA(par = par_v, fn = sa_fit_f2, lower = par_low, upper = par_upp
	      , control = list(
	      	#maxit = 10
	      	max.time = 60
	      	, smooth = F
	      	, simple.function = F))

trace_ff <- data.frame(sao$trace)$function.value
plot(trace_ff, type = "l")
percent(- sao$value)
final_vector <- c((sao$par >= threshold), T)
names(sampleA)[final_vector]
final_sample <- as.data.frame(sampleA[, final_vector])

Sys.time() - start


# check model
lm_m <- lm(formula = z ~ .^2,
	    data = sampleA[, final_vector],
	    model = T)

summary(lm_m)






### greedy search with wrapper GBM validation fitness

library(FSelector)
library(gbm)

sampleA <- bunny_dat[sample(nrow(bunny_dat), nrow(bunny_dat) / 2, replace = F), c("x",
											     "y",
											     "input_noise_1",  
											     "input_noise_2",
											     "input_noise_3",
											     "input_noise_4",
											     "input_noise_5",
											     "input_noise_6", 
											     "input_noise_7",
											     "input_noise_8",
											     "input_noise_9",
											     "input_noise_10",
											     "z")]

sampleB <- bunny_dat[!row.names(bunny_dat) %in% rownames(sampleA), c("x",
									      "y",
									      "input_noise_1",  
									      "input_noise_2",
									      "input_noise_3",
									      "input_noise_4",
									      "input_noise_5",
									      "input_noise_6", 
									      "input_noise_7",
									      "input_noise_8",
									      "input_noise_9",
									      "input_noise_10",
									      "z")]

gbm_fit <- function(subset){
	dat <- sampleA[, c(subset, "z")]
	
	gbm_fit <- gbm(formula = z ~.,
			 distribution = "gaussian",
			 data = dat,
			 n.trees = 50,
			 interaction.depth = 10,
			 n.minobsinnode = 100,
			 shrinkage = 0.1,
			 bag.fraction = 0.9,
			 train.fraction = 1,
			 n.cores = 7)
	
	gbm_predict <- predict(gbm_fit,
				 newdata = sampleB,
				 n.trees = 50)
	
	r_sq_validate <- 1 - sum((sampleB$z - gbm_predict) ^ 2) / sum((sampleB$z - mean(sampleB$z)) ^ 2)
	
	print(subset)
	print(r_sq_validate)
	return(r_sq_validate)
}


subset <- backward.search(attributes = names(sampleA)[1:(ncol(sampleA) - 1)], eval.fun = gbm_fit)
subset <- forward.search(attributes = names(sampleA)[1:(ncol(sampleA) - 1)], eval.fun = gbm_fit)







############## simlated annelaing search with mutual information fitness function

library(infotheo) # measured in nats, converted to bits
library(scales)
library(GenSA)

sampleA <- bunny_dat[sample(nrow(bunny_dat), nrow(bunny_dat) / 2, replace = F), c("x",
											     "y",
											     "input_noise_1",  
											     "input_noise_2",
											     "input_noise_3",
											     "input_noise_4",
											     "input_noise_5",
											     "input_noise_6", 
											     "input_noise_7",
											     "input_noise_8",
											     "input_noise_9",
											     "input_noise_10",
											     "z")]

sampleB <- bunny_dat[!row.names(bunny_dat) %in% rownames(sampleA), c("x",
									      "y",
									      "input_noise_1",  
									      "input_noise_2",
									      "input_noise_3",
									      "input_noise_4",
									      "input_noise_5",
									      "input_noise_6", 
									      "input_noise_7",
									      "input_noise_8",
									      "input_noise_9",
									      "input_noise_10",
									      "z")]

# discretize all variables

dat <- sampleA

disc_levels <- 5

for (i in 1:13){
	
	naming <- paste(names(dat[i]), 'discrete', sep = "_")
	dat[, eval(naming)] <- discretize(dat[, eval(names(dat[i]))], disc = "equalfreq", nbins = disc_levels)[,1]
	
}

sampleA <- dat[, 14:26]

#calculate parameters
predictor_number <- dim(sampleA)[2] - 1
sample_size <- dim(sampleA)[1]
par_v <- runif(predictor_number, min = 0, max = 1)
par_low <- rep(0, times = predictor_number)
par_upp <- rep(1, times = predictor_number)


#load functions to memory
shuffle_f_inp <- function(x = data.frame(), iterations_inp, quantile_val_inp){
	
	mutins <- c(1:iterations_inp)
	
	for (count in 1:iterations_inp){
		
		xx <- data.frame(1:dim(x)[1])
		
		for (count1 in 1:(dim(x)[2] - 1)){
			
			y <- as.data.frame(x[, count1])
			y$count <- sample(1 : dim(x)[1], dim(x)[1], replace = F)
			y <- y[order(y$count), ]
			
			xx <- cbind(xx, y[, 1])
			
		}
		
		mutins[count] <- multiinformation(xx[, 2:dim(xx)[2]])
		
	}
	
	quantile(mutins, probs = quantile_val_inp)
	
}


shuffle_f <- function(x = data.frame(), iterations, quantile_val){
	
	height <- dim(x)[1]
	
	mutins <- c(1:iterations)
	
	for (count in 1:iterations){
		
		x$count <- sample(1 : height, height, replace = F)
		
		y <- as.data.frame(c(x[dim(x)[2] - 1], x[dim(x)[2]]))
		y <- y[order(y$count), ]
		
		x[dim(x)[2]] <- NULL
		x[dim(x)[2]] <- NULL
		x$dep <- y[, 1]
		
		rm(y)
		
		receiver_entropy <- entropy(x[, dim(x)[2]])
		received_inf <- mutinformation(x[, 1 : dim(x)[2] - 1], x[, dim(x)[2]])
		corr_ff <- received_inf / receiver_entropy
		
		mutins[count] <- corr_ff
		
	}
	
	quantile(mutins, probs = quantile_val)
	
}

############### the fitness function
fitness_f <- function(par){
	
	indexes <- c(1:predictor_number)
	
	for (i in 1:predictor_number){
		if (par[i] >= threshold) {
			indexes[i] <- i
		} else {
			indexes[i] <- 0
		}
	}
	
	local_predictor_number <- 0
	for (i in 1:predictor_number){
		if (indexes[i] > 0) {
			local_predictor_number <- local_predictor_number + 1
		}
	}
	
	if (local_predictor_number > 1) {
		
		sampleAf <- as.data.frame(sampleA[, c(indexes[], dim(sampleA)[2])])
		
		pred_entrs <- c(1:local_predictor_number)
		
		for (count in 1:local_predictor_number){
			
			pred_entrs[count] <- entropy(sampleAf[count])
			
		}
		
		max_pred_ent <- sum(pred_entrs) - max(pred_entrs)
		
		pred_multiinf <- multiinformation(sampleAf[, 1:dim(sampleAf)[2] - 1])
		
		pred_multiinf <- pred_multiinf - shuffle_f_inp(sampleAf, iterations_inp, quantile_val_inp)
		
		if (pred_multiinf < 0){
			
			pred_multiinf <- 0
			
		}
		
		pred_mult_perc <- pred_multiinf / max_pred_ent
		
		inf_corr_val <- shuffle_f(sampleAf, iterations, quantile_val)
		
		receiver_entropy <- entropy(sampleAf[, dim(sampleAf)[2]])
		received_inf <- mutinformation(sampleAf[, 1:local_predictor_number], sampleAf[, dim(sampleAf)[2]])
		
		if (inf_corr_val - (received_inf / receiver_entropy) < 0){
			
			fact_ff <- (inf_corr_val - (received_inf / receiver_entropy)) * (1 - pred_mult_perc)
		} else {
			
			fact_ff <- inf_corr_val - (received_inf / receiver_entropy)
			
		}
		
	} else if (local_predictor_number == 1) {
		
		sampleAf<- as.data.frame(sampleA[, c(indexes[], dim(sampleA)[2])])
		
		inf_corr_val <- shuffle_f(sampleAf, iterations, quantile_val)
		
		receiver_entropy <- entropy(sampleAf[, dim(sampleAf)[2]])
		received_inf <- mutinformation(sampleAf[, 1:local_predictor_number], sampleAf[, dim(sampleAf)[2]])
		
		fact_ff <- inf_corr_val - (received_inf / receiver_entropy)
		
	} else  {
		fact_ff <- 0
	}
	return(fact_ff)
}


########## estimating threshold for variable inclusion

iterations = 1
quantile_val = 1

iterations_inp = 1
quantile_val_inp = 1

levels_arr <- numeric()
for (i in 1:predictor_number){
	levels_arr[i] <- length(unique(sampleA[, i]))
}

mean_levels <- mean(levels_arr)

optim_var_num <- log(x = sample_size / 100, base = round(mean_levels, 0))

if (optim_var_num / predictor_number < 1){
	
	threshold <- 1 - optim_var_num / predictor_number
	
} else {
	
	threshold <- 0.5
	
}


#run feature selection

start <- Sys.time()

sao <- GenSA(par = par_v, fn = fitness_f, lower = par_low, upper = par_upp
	      , control = list(
	      	#maxit = 10
	      	max.time = 1200
	      	, smooth = F
	      	, simple.function = F))

trace_ff <- data.frame(sao$trace)$function.value
plot(trace_ff, type = "l")
percent(- sao$value)
final_vector <- c((sao$par >= threshold), T)
names(sampleA)[final_vector]
final_sample <- as.data.frame(sampleA[, final_vector])

Sys.time() - start




######################## THE END



















