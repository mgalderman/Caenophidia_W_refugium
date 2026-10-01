library(multcomp)
### Set up density values ----

flTE <- read.table('~/Library/CloudStorage/GoogleDrive-schieldlab@gmail.com/My Drive/projects/caenophidia_W_refugium/analysis/toxicity_hypothesis//fl-TE_density.txt',header=T)

auto.den <- flTE$auto.flTE/flTE$auto.bp
z.den <- flTE$z.flTE/flTE$z.bp
w.den <- flTE$w.flTE/flTE$w.bp

### Perform ANOVA and Tukey post-hoc tests ----

# Rearrange density data into long format

long_format <- list(auto = auto.den, chrZ = z.den, chrW = w.den)
long_format <- stack(long_format)
colnames(long_format) <- c("density","chromosome")
print(long_format)

density_aov <- aov(density ~ chromosome, data = long_format)
summary(density_aov)

# Examine the normality of data and residuals

hist(density_aov$residuals)
hist(long_format$density)
plot(density_aov$residuals)
boxplot(density ~ chromosome, data = long_format)

# Tukey's post-hoc

post_test <- glht(density_aov, linfct = mcp(chromosome = "Tukey"))
summary(post_test)

### Plot comparison of distributions ----

boxplot(auto.den, z.den, w.den, names = c('Auto','Z','W'),ylab='Full-length TE density')
segments(x0 = 1, y0 = .00005, x1 = 2.8, y1 = .00005)
text(x = 2, y = .000051, labels = "*", cex = 2)

segments(x0 = 2, y0 = .00004, x1 = 2.8, y1 = .00004)
text(x = 2.4, y = .000041, labels = "*", cex = 2)

# output at 4 x 3

### Plot simple visualization of toxicity indexes ----

tox <- c(0.0003,-0.0269,-0.0258,-0.0241,-0.0165,-0.0172,-0.0052,-0.0129,-0.0321)
tox <- sort(tox)
birds <- c(0.06,0.19,0.28,0.33,0.37,0.83)
all <- c(tox,birds)

plot(tox,pch=20)
plot(all,pch=20,xlab='Rank order',ylab='Toxicity index')

# output at 4 x 4
