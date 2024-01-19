# Set the working directory to the specified path
setwd('~/Documents/IMB/Johanna/')
set.seed(584930)

# Load required libraries
library("RColorBrewer")  # Library for color palettes
library(pROC)  # Library for ROC curve analysis

# Read the dataset from a CSV file
feats <- read.csv(file = 'data/3did_Interactome_ProtCID_AF_IUPred_metrics_mc_set_for_ML_without_synthetic_constructs.csv')

# Define a range of numeric variables (columns) in the dataset
numvars <- 2:18

# Name variables for plots
varnames <- c("3did n structures","3did n interchain","3did n intrachain","3did score",                   
              "3did z score","3did n res-res contacts","3did fraction interchain structures","HuRI z score",                  
              "BioPlex z score","ProtCID n crystal forms","ProtCID min seq identity","ProtCID n proteins",       
              "AF2 model confidence","AF2 DockQ","Avg IUPred Interface full","Avg IUPred Interface cryst",    
              "DDI approved")

names(varnames) <- c("Num.structures","Num.interchain.structures","Num.intrachain.structures","X3did.score",                   
                     "X3did.z.score","Num.res.res.contacts","Fraction.structures.interchain","HuRI.z.score",                  
                     "BioPlex.z.score","Num.crystal.forms.in.cluster","Min.seq.identity.in.cluster","Num.proteins.in.cluster",       
                     "Model.confidence","DockQ","Avg.IUPred.Iface.full.seq","Avg.IUPred.Iface.cryst.seq",    
                     "DDI.approved")


# Find variables with missing values in the dataset
missvars <- which(apply(feats, 2, function(x) any(is.na(x))))

# Generate a correlation heatmap and save it as a PDF file
pdf('model/figs/corr0.pdf', width = 6, height = 6)
heatmap(
  cor(feats[, numvars], use = 'pairwise.complete.obs'),
  labRow = sapply(names(feats[,numvars]), function(x) varnames[x]), 
  scale = 'none',
  col = rev(brewer.pal(11, "RdBu")),
  margins = c(15, 15),
  labCol = F
)
dev.off()

# Define a subset of numeric variables to focus on
rednumvars <- which(names(feats) %in% c(
  "Num.intrachain.structures", "Num.interchain.structures", "Fraction.structures.interchain", "X3did.z.score",
  "HuRI.z.score", "BioPlex.z.score", "Min.seq.identity.in.cluster", "Model.confidence", "DockQ",
  "Avg.IUPred.Iface.full.seq", "DDI.approved"
))

# Filter variables with missing values from the reduced numeric variables
fullrednumvars <- rednumvars[!rednumvars %in% missvars]

# Filter variables with missing values from the full numeric variables
fullnumvars <- numvars[!numvars %in% missvars]

# Generate another correlation heatmap for the selected variables and save it as a PDF
pdf('model/figs/corr.pdf', width = 6, height = 6)
heatmap(
  cor(feats[, rednumvars], use = 'pairwise.complete.obs'),
  labRow = sapply(names(feats[,rednumvars]), function(x) varnames[x]), 
  scale = 'none',
  col = rev(brewer.pal(11, "RdBu")),
  margins = c(16, 16),
  labCol = F
)
dev.off()

# Perform logistic regression with selected variables
summary(glm(DDI.approved ~ ., binomial(link = "logit"), feats[, numvars]))
summary(glm(DDI.approved ~ ., binomial(link = "logit"), feats[, fullrednumvars]))

# Stepwise variable selection using logistic regression
step(glm(DDI.approved ~ ., binomial(link = "logit"), feats[, fullrednumvars]))

# Select a subset of variables for logistic regression
selvars <- c('DDI.approved', "X3did.z.score", "DockQ", "Avg.IUPred.Iface.full.seq")

# Perform logistic regression with the selected variables
summary(glm(DDI.approved ~ ., binomial(link = "logit"), feats[, selvars]))

# Create a pairs plot for the selected variables and save it as a PDF
pdf('model/figs/pairs.pdf')
pairs(feats[, selvars], cex.labels = .9)
dev.off()

# Define and fit logistic regression models
# minimal model
mod0 <- glm(
  formula = DDI.approved ~ X3did.z.score,
  family = binomial(link = "logit"),
  data = feats
)
summary(mod0)

# full model
mod1 <- glm(
  formula = DDI.approved ~ .,
  family = binomial(link = "logit"),
  data = feats[, fullrednumvars]
)
summary(mod1)

# reduced model
mod2 <- glm(
  formula = DDI.approved ~ .,
  family = binomial(link = "logit"),
  data = feats[, selvars]
)
summary(mod2)

# Perform ANOVA on the fitted models
anova(mod0, mod2, mod1, test = 'Chisq')

# Create a data frame containing fitted values and DDI type
scores <- data.frame(
  DDI.type = feats$DDI.type,
  minimal_model = mod0$fitted.values,
  full_model = mod1$fitted.values,
  reduced_model = mod2$fitted.values
)

# Sort DDI type values and write the data to a CSV file
scores$sorted_DDI_type <- sapply(lapply(strsplit(scores$DDI.type, '_'), sort), paste, collapse = '_')
write.csv(scores, file = 'model/results/fitted_values.csv', row.names = F)

# Leave-one-out (LOO) cross-validation
mod0$data$loo <- mod1$data$loo <- mod2$data$loo <- NA

for (i in 1:nrow(feats)) {
  submod0 <- glm(
    formula = DDI.approved ~ X3did.z.score,
    family = binomial(link = "logit"),
    data = feats[-i, ]
  )
  mod0$data$loo[i] <- predict.glm(submod0, newdata = feats[i,])
  
  submod1 <- glm(
    formula = DDI.approved ~ .,
    family = binomial(link = "logit"),
    data = feats[-i, fullrednumvars]
  )
  mod1$data$loo[i] <- predict.glm(submod1, newdata = feats[i, fullrednumvars])
  
  submod2 <- glm(
    formula = DDI.approved ~ .,
    family = binomial(link = "logit"),
    data = feats[-i, selvars]
  )
  mod2$data$loo[i] <- predict.glm(submod2, newdata = feats[i, selvars])
  
}

# Generate ROC curves and save them as a PDF
pdf('model/figs/roc.pdf')
par(mar = c(3, 3, 3, 1), mgp = c(2, .5, 0), tck = -.01, cex = 1.5)
pROC_obj0 <- with(mod0$data, roc(DDI.approved, loo, plot = T))
pROC_obj1 <- with(mod1$data, roc(DDI.approved, loo, plot = T, add = T, col = 'red'))
pROC_obj2 <- with(mod2$data, roc(DDI.approved, loo, plot = T, add = T, col = 'orange'))
legend(x = .5, y = .5, legend = c('full', 'reduced', 'Z score'), col = c('red', 'orange', 'black'), lty = 1)
dev.off()

#custom ROC plot
{
pdf('model/figs/roc2.pdf', width = 2, height = 2)
par(mar = c(3, 3, 1, 0), mgp = c(2, .5, 0), cex = .9)
with(coords(pROC_obj0), plot(sensitivity~specificity, xlim=c(1,0), type='l', bty='l', las=1))
with(coords(pROC_obj1), lines(sensitivity~specificity, type='l', col='red'))
with(coords(pROC_obj2), lines(sensitivity~specificity,type='l', col='orange'))
legend(x = .7, y = .5, legend = c('full', 'reduced', 'Z score'), col = c('red', 'orange', 'black'), lty = 1, cex=.7)
dev.off()
}


# Perform ROC tests using bootstrap method
roc.test(pROC_obj0, pROC_obj1, method = 'bootstrap')
roc.test(pROC_obj0, pROC_obj2, method = 'bootstrap')
roc.test(pROC_obj1, pROC_obj2, method = 'bootstrap')

# Get the best threshold's sensitivity and specificity
coords(pROC_obj0,"best")
coords(pROC_obj1,"best")
coords(pROC_obj2,"best")

# Perform the Hosmer-Lemeshow goodness-of-fit test for logistic regression models
glmtoolbox::hltest(mod0)
glmtoolbox::hltest(mod1)
glmtoolbox::hltest(mod2)

# Calculate specific quantiles and Z-scores
plogis(2.3, location = -mod0$coefficients[1] / mod0$coefficients[2], scale = 1 / mod0$coefficients[2])
qlogis(.8, location = -mod0$coefficients[1] / mod0$coefficients[2], scale = 1 / mod0$coefficients[2])
qlogis(.9, location = -mod0$coefficients[1] / mod0$coefficients[2], scale = 1 / mod0$coefficients[2])

{
# Analyze empirical data for minimal model and create a PDF plot
percs <- sapply(qlogis(mod0$fitted.values), function(x) sum(x > quantile(qlogis(mod0$fitted.values), 1:9/10)))
mean.x <- aggregate(mod0$data$X3did.z.score, list(percs), mean)
mean.y <- aggregate(mod0$data$DDI.approved, list(percs), mean)

pdf('model/figs/empirical_minimal.pdf', width = 2, height = 2)
par(mar = c(3, 3, 1, 0), mgp = c(2, .5, 0), cex = .9)
curve(
  plogis(x, location = -mod0$coefficients[1] / mod0$coefficients[2], scale = 1 / mod0$coefficients[2]),
  -2, 8, las =1,
  ylim = c(0, 1),
  ylab = 'DDI Approval Rate',
  xlab = '3did Z score',
  bty = 'l'
)
segments(x0 = mean.x$x, y0 = mean.y$x, y1 = plogis(mean.x$x, location = -mod0$coefficients[1] / mod0$coefficients[2], scale = 1 / mod0$coefficients[2]))
points(x = mean.x$x, y = mean.y$x, col = 'white', pch = 20)
points(x = mean.x$x, y = mean.y$x)
abline(v = 2.3, col = 'blue', lty = 2)
abline(v = qlogis(.8, location = -mod0$coefficients[1] / mod0$coefficients[2], scale = 1 / mod0$coefficients[2]), col = 'purple', lty = 2)
abline(v = qlogis(.9, location = -mod0$coefficients[1] / mod0$coefficients[2], scale = 1 / mod0$coefficients[2]), col = 'red', lty = 2)
dev.off()
}

{
# Define and analyze empirical data with the full model
percs <- sapply(qlogis(mod1$fitted.values), function(x) sum(x > quantile(qlogis(mod1$fitted.values), 1:9/10)))
mean.x <- aggregate(qlogis(mod1$fitted.values), list(percs), mean)
mean.y <- aggregate(mod1$data$DDI.approved, list(percs), mean)
table(percs)

# Create a PDF plot for the empirical data with the full model
pdf('model/figs/empirical_full.pdf', width = 2, height = 2)
par(mar = c(3, 3, 1, 0), mgp = c(2, .5, 0), cex = .9)
curve(
  plogis(x),
  -11, 5,
  ylim = c(0, 1),
  ylab = 'DDI Approval Rate',
  xlab = 'Full DDI score',
  bty = 'l', las=1
)
segments(x0 = mean.x$x, y0 = mean.y$x, y1 = plogis(mean.x$x))
points(x = mean.x$x, y = mean.y$x, col = 'white', pch = 20)
points(x = mean.x$x, y = mean.y$x)
# abline(v = qlogis(.8), col = 'purple', lty = 2)
# abline(v = qlogis(.9), col = 'red', lty = 2)
dev.off()
}



{
# Define and analyze empirical data with the reduced model
percs <- sapply(qlogis(mod2$fitted.values), function(x) sum(x > quantile(qlogis(mod2$fitted.values), 1:9/10)))
mean.x <- aggregate(qlogis(mod2$fitted.values), list(percs), mean)
mean.y <- aggregate(mod2$data$DDI.approved, list(percs), mean)
table(percs)

# Create a PDF plot for the empirical data with the reduced model
pdf('model/figs/empirical_reduced.pdf', width = 2, height = 2)
par(mar = c(3, 3, 1, 0), mgp = c(2, .5, 0), cex = .9)
curve(
  plogis(x),
  -5, 5,
  ylim = c(0, 1),
  ylab = 'DDI Approval Rate',
  xlab = 'Reduced DDI score',
  bty = 'l', las=1
)
segments(x0 = mean.x$x, y0 = mean.y$x, y1 = plogis(mean.x$x))
points(x = mean.x$x, y = mean.y$x, col = 'white', pch = 20)
points(x = mean.x$x, y = mean.y$x)
# abline(v = qlogis(.8), col = 'purple', lty = 2)
# abline(v = qlogis(.9), col = 'red', lty = 2)
dev.off()
}


# Calculate odds ratios and confidence intervals
oric <- function(model, alpha = 0.05) {
  coefs <- summary(model)$coefficients
  quan <- 1 - alpha / 2
  return(data.frame(
    avg = exp(coefs[, 1]),
    lower = exp(coefs[, 1] - qnorm(quan) * coefs[, 2]),
    upper = exp(coefs[, 1] + qnorm(quan) * coefs[, 2])
  ))
}

# Calculate odds ratios and confidence intervals for mod1
oric(mod1)
summary(mod1)

oric(mod0)
