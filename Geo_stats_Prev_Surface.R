#Creating mock survey sites and surface layers to predict a prevalence surface layer. 

#1. create dummy survye sites. 

library(geosphere)
library(MASS)
library(dplyr)

set.seed(101)

# 1. Simulate cluster centers
n_clusters <- 6
cluster_centers <- data.frame(
  cluster = 1:n_clusters,
  lon = runif(n_clusters, -1.5, -1.0),
  lat = runif(n_clusters, 50.6, 50.8)
)

# 2. Simulate survey points around clusters
n_points <- 300
survey_points <- data.frame(
  cluster = sample(1:n_clusters, n_points, replace = TRUE)
) %>%
  left_join(cluster_centers, by = "cluster") %>%
  mutate(
    lon = lon + rnorm(n(), sd = 0.04),  # jitter around center
    lat = lat + rnorm(n(), sd = 0.05)
  )

# 3. Distance from cluster center (radial)
survey_points <- survey_points %>%
  rowwise() %>%
  mutate(
    dist_km = distGeo(
      p1 = c(lon, lat),
      p2 = cluster_centers[cluster, c("lon", "lat")]
    ) / 1000
  ) %>%
  ungroup()

# 4. Normalize distance and longitude
survey_points <- survey_points %>%
  mutate(
    norm_dist = dist_km / max(dist_km),
    norm_lon = (lon - min(lon)) / (max(lon) - min(lon))
  )

# 5. Add spatially structured noise
coords <- cbind(survey_points$lon, survey_points$lat)
D <- as.matrix(dist(coords))
phi <- 0.01
nugget_variance <- 0.1  # for example, variance of nugget effect
spatial_variance <- 0.5
Sigma <- spatial_variance * exp(-D / phi) + diag(nugget_variance, n_points)
spatial_effect <- mvrnorm(1, rep(0, n_points), Sigma)

# 6. Combine both effects in linear predictor
linpred <- -2 + 2.5 * survey_points$norm_dist + 2 * survey_points$norm_lon + spatial_effect
prob_pos <- plogis(linpred)

# 7. Simulate outcome
surveyed <- sample(50:200, n_points, replace = TRUE)
positive <- rbinom(n_points, size = surveyed, prob = prob_pos)

# 8. Final dataset
survey_data <- survey_points %>%
  mutate(
    x = lon,
    y = lat,
    surveyed = surveyed,
    positive = positive,
    prob_pos = prob_pos,
    weights = 1
  ) %>%
  dplyr::select(lon, lat, x, y, cluster, dist_km, norm_dist, norm_lon,
                surveyed, positive, prob_pos, weights)

# Preview
head(survey_data)

survey_data$prev <- survey_data$positive / survey_data$surveyed

write.csv(survey_data, "survey_data_PrevSurface.csv")

# Convert to SpatialPointsDataFrame
library(sp)
coordinates(survey_data) <- ~lon + lat
library(raster)
# Set the same CRS as the raster (if known)
crs(survey_data) <- crs("+proj=longlat +datum=WGS84")

# Define bounding box as extent
bbox_extent <- extent(-1.5, -1.0, 50.6, 50.8)

# Crop survey_data
survey_data <- crop(survey_data, bbox_extent)

library(tmap)

tmap_mode("plot")  # or "view" for interactive mode
#tmap_mode("view")

tm_shape(survey_data) +
  tm_symbols(col = "prev",
             palette = "viridis",
             size = 0.1,
             border.col = "black",
             title.col = "Prevalence",
             style = "cont") +
  tm_layout(legend.outside = TRUE,
            legend.outside.position = "right",
            frame = TRUE) +   # show axis ticks and labels
  tm_compass(type = "8star", position = c("left", "top")) +
  tm_scale_bar(position = c("left", "bottom"))

###now to make the land usgae layer. 

library(raster)
library(sp)

# Your cluster centers (from your earlier survey generation)
# For example, you can recreate them here or extract from survey_data

library(raster)
library(sp)

bbox_iow <- c(xmin = -1.5, xmax = -1.0, ymin = 50.6, ymax = 50.8)
res <- 0.001

# Create base raster with all cells initially agricultural (1)
landuse_raster <- raster(xmn = bbox_iow["xmin"], xmx = bbox_iow["xmax"],
                         ymn = bbox_iow["ymin"], ymx = bbox_iow["ymax"],
                         res = res,
                         crs = CRS("+proj=longlat +datum=WGS84"))
values(landuse_raster) <- 1  # Agricultural everywhere initially

# Cluster centers
set.seed(101)
n_clusters <- 6
cluster_centers <- data.frame(
  cluster = 1:n_clusters,
  lon = runif(n_clusters, -1.5, -1.0),
  lat = runif(n_clusters, 50.6, 50.8))
  
cluster_points <- SpatialPoints(
  coords = cluster_centers[, c("lon", "lat")],
  proj4string = CRS(projection(landuse_raster))
)

# Urban radius in degrees (~0.03 degrees ~3 km)
urban_radius <- 3

# Create urban mask raster initialized with zeros
urban_mask <- raster(landuse_raster)
values(urban_mask) <- 0

# Mark urban areas in urban_mask (value = 1)
for(i in 1:length(cluster_points)) {
  cells <- cellsFromExtent(
    urban_mask,
    extent(
      cluster_points[i,]@coords[1] - urban_radius,
      cluster_points[i,]@coords[1] + urban_radius,
      cluster_points[i,]@coords[2] - urban_radius,
      cluster_points[i,]@coords[2] + urban_radius
    )
  )
  xy_cells <- xyFromCell(urban_mask, cells)
  pt <- coordinates(cluster_points[i,])[1,]
  dists <- spDistsN1(xy_cells, pt, longlat = TRUE)
  urban_cells <- cells[dists <= urban_radius]
  urban_mask[urban_cells] <- 1
}

# Create a vector with x coordinates for all raster cells to build forest gradient
xy_all <- xyFromCell(landuse_raster, 1:ncell(landuse_raster))
x_coords <- xy_all[,1]

# Normalize x_coords from 0 (west) to 1 (east) across raster extent
x_norm <- (x_coords - bbox_iow["xmin"]) / (bbox_iow["xmax"] - bbox_iow["xmin"])

# Define forest probability increasing from west (0) to east (1)
forest_prob <- x_norm

# Increase agricultural probability near urban areas (buffered)
# For example, agricultural probability = 0.9 near urban, fading to 0.5 far away

# Distance to nearest urban cell (in degrees)
dist_to_urban <- distanceFromPoints(landuse_raster, cluster_points)

# Normalize distance: 0 at urban, max distance ~max(dist_to_urban) (convert to 0-1)
max_dist <- max(values(dist_to_urban), na.rm = TRUE)
dist_norm <- values(dist_to_urban) / max_dist
dist_norm[is.na(dist_norm)] <- 1  # treat NA as farthest

# Agricultural probability increases near urban (inverse of normalized distance)
agri_prob <- 0.3 + (1 - dist_norm) * 0.4  # from 0.5 to 0.9 near urban

# Now combine probabilities:
# Assign class by max probability: Agricultural (1), Forest (3), Urban (2)

# Urban cells are already known
urban_cells <- which(values(urban_mask) == 1)

# Create a vector for classes, start with agricultural (1)
landuse_values <- rep(1, ncell(landuse_raster))

# Assign forest probabilistically (if forest_prob > agri_prob)
forest_cells <- which(forest_prob > agri_prob)
landuse_values[forest_cells] <- 3

# Urban cells override to 2
landuse_values[urban_cells] <- 2

# Assign back to raster
landuse_final <- landuse_raster
values(landuse_final) <- landuse_values

# Plot with colors
landuse_cols <- c("yellow", "gray", "darkgreen")
landuse_legend <- c("Agricultural (1)", "Urban (2)", "Forest (3)")

plot(landuse_final, col = landuse_cols, legend = FALSE, main = "Land Use with Agricultural Gradient near Urban and Forest Gradient East")
legend("topright", legend = landuse_legend, fill = landuse_cols)

#plot the points on the raster

library(ggplot2)
library(raster)
library(dplyr)

# Convert raster to data frame for ggplot
r_df <- as.data.frame(landuse_final, xy = TRUE)
colnames(r_df) <- c("x", "y", "landuse")

survey_df <- as.data.frame(survey_data)

# Define colors and labels again (if not already)
landuse_cols <- c("1" = "orange", "2" = "grey40", "3" = "darkgreen")
landuse_legend <- c("Rural", "Urban", "Forest")

tmap_mode("plot")

tm_basemap() +
  tm_shape(landuse_final) +
  tm_raster(style = "cat",
            palette = landuse_cols,
            labels = landuse_legend,
            title = "Land Use") +
  tm_shape(survey_data) +
  tm_symbols(col = "prev",
             palette = "viridis",
             size = 0.1,
             border.col = "black",
             title.col = "Prevalence",
             style = "cont") +
  tm_layout(
    legend.outside = TRUE,
    legend.outside.position = "right",  # use "bottom" if preferred
    legend.title.size = 1.2,
    legend.text.size = 0.8,
    legend.stack = "vertical"           # stack legends vertically
  ) +
tm_compass(type = "8star", position = c("left", "top")) +
  tm_scale_bar(position = c("left", "bottom"))

#####Extract Data points#####

library(raster)
library(sp)

# Assuming landuse_final raster with values 1=Agri, 2=Urban, 3=Forest
# And water_raster with 4=Water, NA elsewhere

# Create binary rasters for each feature
urban_dist_raster <- landuse_final
values(urban_dist_raster) <- ifelse(values(landuse_final) == 2, 1, NA)
dist_to_urban <- distance(urban_dist_raster)

forest_dist_raster <- landuse_final
values(forest_dist_raster) <- ifelse(values(landuse_final) == 3, 1, NA)
dist_to_forest <- distance(forest_dist_raster)

agri_dist_raster <- landuse_final
values(agri_dist_raster) <- ifelse(values(landuse_final) == 1, 1, NA)
dist_to_agri <- distance(agri_dist_raster)

# Extract distance values at survey points
# survey_data assumed to be SpatialPointsDataFrame
survey_coords <- coordinates(survey_data)

survey_data$dist_to_urban <- extract(dist_to_urban, survey_coords)
survey_data$dist_to_forest <- extract(dist_to_forest, survey_coords)
survey_data$dist_to_agri <- extract(dist_to_agri, survey_coords)

# Convert to data frame to view
survey_df <- as.data.frame(survey_data)

head(survey_df[, c("dist_to_urban", "dist_to_forest", "dist_to_agri")])

#explore the relationship between the dependent and independent variables. 

library(ggplot2)

#logit transformation
survey_df$logit_prev <- with(survey_df, log((positive + 0.5) / (surveyed - positive + 0.5)))

ggplot(survey_df, aes(x = dist_to_urban, y = logit_prev)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  labs(title = "Logit vs. Distance to Urban", x = "Distance to Urban (m)", y = "Logit")

ggplot(survey_df, aes(x = dist_to_forest, y = logit_prev)) +
  geom_point() +
  geom_smooth(method = "lm", color = "forestgreen") +
  labs(title = "Logit vs. Distance to Forest", x = "Distance to Forest (m)", y = "Logit")

ggplot(survey_df, aes(x = dist_to_agri, y = logit_prev)) +
  geom_point() +
  geom_smooth(method = "lm", color = "darkorange") +
  labs(title = "Logit vs. Distance to Agriculture", x = "Distance to Agriculture (m)", y = "Logit")

library(gstat)
library(sp)
library(sf)

coordinates(survey_df) <- ~x + y
proj4string(survey_df) <- CRS("+proj=longlat +datum=WGS84")  # set projection

# Plot
library(gstat)
library(ggplot2)
library(dplyr)

# Empirical variogram
vgm_emp <- gstat::variogram(logit_prev ~ 1, data = survey_df)

# Fit theoretical variogram model (exponential in this case)
vgm_fit <- fit.variogram(vgm_emp, model = vgm("Exp"))

# Extract model parameters
nugget <- vgm_fit$psill[vgm_fit$model == "Nug"]
partial_sill <- vgm_fit$psill[vgm_fit$model != "Nug"]
sill <- nugget + partial_sill
range_val <- vgm_fit$range[vgm_fit$model != "Nug"]

# Create fitted curve from the model
vgm_fit_line <- data.frame(
  dist = seq(min(vgm_emp$dist), max(vgm_emp$dist), length.out = 100)
)

vgm_fit_line$gamma <- variogramLine(vgm_fit, dist_vector = vgm_fit_line$dist)$gamma

ggplot(vgm_emp, aes(x = dist, y = gamma)) +
  geom_point(color = "black", size = 2) +
  geom_line(data = vgm_fit_line, aes(x = dist, y = gamma), color = "red", linewidth = 1) +
  labs(title = "Empirical and Fitted Semivariogram",
       x = "Distance",
       y = "Semivariance") +
  annotate("text", x = range_val * 1, y = sill * 0.6,
           label = paste0("Nugget: ", round(nugget, 3)), hjust = 0, color = "black") +
  annotate("text", x = range_val * 1, y = sill * 0.55,
           label = paste0("Partial Sill: ", round(partial_sill, 3)), hjust = 0, color = "black") +
  annotate("text", x = range_val * 1, y = sill * 0.5,
           label = paste0("Range: ", round(range_val, 3)), hjust = 0, color = "black")

#there is spatial autocorrelation within the dataset so a geostats model is required.

#run an model with the covariates using standard regression. 
#Then run a semi-variogram on the residuals. 

lm_fit_ur <- lm(logit_prev ~ dist_to_urban, data = survey_df)
summary(lm_fit_ur)
lm_fit_fo <- lm(logit_prev ~ dist_to_forest, data = survey_df)
summary(lm_fit_fo)
lm_fit_ag <- lm(logit_prev ~ dist_to_agri, data = survey_df)
summary(lm_fit_ag)

lm_fit <- glm(cbind(positive, surveyed - positive) ~ dist_to_urban,
                         data = survey_df,
                         family = binomial)
summary(lm_fit)

#now running a semi-variogram on the residuals to explore is there is remaing statial correlation. 

# Extract residuals and add to your dataframe
survey_df$residuals <- residuals(lm_fit)

library(gstat)

# Empirical semivariogram of residuals
vgm_emp_resid <- gstat::variogram(residuals ~ 1, data = survey_df)

# Fit a model to the semivariogram
#vgm_fit_resid <- fit.variogram(vgm_emp_resid, model = vgm("Exp"))

# Create prediction line from the fitted variogram
#vgm_fit_line_resid <- data.frame(
#  dist = seq(min(vgm_emp_resid$dist), max(vgm_emp_resid$dist), length.out = 100)
#)
#vgm_fit_line_resid$gamma <- variogramLine(vgm_fit_resid, dist_vector = vgm_fit_line_resid$dist)$gamma

# Plot
library(ggplot2)

ggplot(vgm_emp_resid, aes(x = dist, y = gamma)) +
  geom_point(color = "black", size = 2) +
  labs(title = "Semivariogram of Model Residuals",
       x = "Distance",
       y = "Semivariance")

#still looks slightly like there is some residual spatial auto correlation. 

library(geoR)
library(PrevMap)

# Fit a geostatistical logistic model

# Extract data and create components
data <- as.data.frame(lm_fit$data)
data$weights <- lm_fit$prior.weights  # typically "surveyed"
data$infected <- lm_fit$y        # "positive"

# Logit transformation with 0.5 adjustment for numerical stability
data$logitp <- log((data$infected + 0.5) / (data$weights - data$infected + 0.5))

# Original formula (usually ~1, or with covariates)
f <- formula(lm_fit)

# For MLE on logit-prevalence
f_linear <- update(f, logitp ~ .)

# For final binomial spatial model
f_binom <- update(f, infected ~ .)

data$x <- coordinates(survey_df)[,1]
data$y <- coordinates(survey_df)[,2]

fit_MLE <- linear.model.MLE(
  formula = f_linear,
  coords = ~ x + y,
  data = data,
  start.cov.pars = c(30, 0.2),  # you can adjust this
  kappa = 0.5,
  messages = FALSE
)

# Combine fixed effect and spatial params
par0 <- c(coef(lm_fit), coef(fit_MLE)[-c(1:length(coef(lm_fit)))]) 

# Initial values for spatial process (phi, tau2/sigma2)
init_pars <- c(par0["phi"], par0["tau^2"] / par0["sigma^2"])

#control setting for MCML
c.mcmc <- control.mcmc.MCML(
  n.sim = 10000, 
  burnin = 2000, 
  thin = 8,
  h = 1.65 / (nrow(data) ^ (1 / 6))  # bandwidth for proposal dist
)

#Fit spatial binomial GLM
fit_spatial <- binomial.logistic.MCML(
  formula = f_binom,
  units.m = ~ weights,
  coords = ~ x + y,
  par0 = par0,
  data = data,
  control.mcmc = c.mcmc,
  kappa = 0.5,
  start.cov.pars = init_pars,
  plot.correlogram = FALSE,
  messages = FALSE
)

summary(fit_spatial)
practicalRange(cov.model = "exponential", phi = coef(fit_spatial)[["phi"]])

#Now to make a surface layer to predict prevalence onto

library(sf)
library(raster)

# Define bbox in lon/lat (WGS84)
bbox_extent <- extent(-1.5, -1.0, 50.6, 50.8)

# Convert to sf polygon
bbox_poly <- st_as_sf(as(bbox_extent, "SpatialPolygons"))
st_crs(bbox_poly) <- 4326  # Set CRS to WGS84

# Now generate the grid with 5 km spacing in meters
pred_points <- st_make_grid(bbox_poly, cellsize = 0.0025, what = "centers")

# Convert to sf object
pred_points_sf <- st_sf(geometry = pred_points)

# Plot
plot(st_geometry(bbox_poly), col = NA, border = "red")
plot(st_geometry(pred_points_sf), add = TRUE, cex = 0.1, pch = 1, col = "blue")

tm_shape(bbox_poly) +
  tm_borders(col = "red", lwd = 2) +
  tm_shape(pred_points_sf) +
  tm_symbols(size = 0.1, col = "blue", shape = 1) +
  tm_layout(
    legend.show = FALSE,
    frame = TRUE
  ) +
  tm_compass(type = "8star", position = c("left", "top")) +
  tm_scale_bar(position = c("left", "bottom"))

library(raster)
library(sf)

# Extract values at prediction points
pred_points_sf$dist_to_urban <- extract(dist_to_urban, pred_points_sf)
pred_points_sf$dist_to_forest <- extract(dist_to_forest, pred_points_sf)
pred_points_sf$dist_to_agri <- extract(dist_to_agri, pred_points_sf)

pred_points_df <- as.data.frame(pred_points_sf)
pred_points_df <- pred_points_df[2:4]

# Generate spatial predictions
c.mcmc <- control.mcmc.MCML(n.sim = 10000, burnin = 2000, thin = 8,
                            h = (1.65) / (nrow(loaloa) ^ (1/6)))


pred_prevalence <- spatial.pred.binomial.MCML(object = fit_spatial, grid.pred = st_coordinates(pred_points_sf), 
                           predictors = pred_points_df, 
                           control.mcmc = c.mcmc, 
                           type = "marginal", 
                           scale.predictions = "prevalence", 
                           quantiles = c(0.025, 0.5, 0.975), 
                           standard.errors = TRUE,
                           messages = F)

# Generate a raster and plot the results
df <- data.frame(st_coordinates(pred_points_sf), 
                 mean = pred_prevalence$prevalence$predictions * 100)

#reverse logit to get prevalence.
df$prev <- plogis(df$mean)

prevalence <- rasterFromXYZ(df, crs = crs(pred_points_sf))

plot(prevalence$prev)

# Plot
tmap_mode("plot")

tm_basemap("OpenStreetMap") +
  tm_shape(prevalence$prev) +
  tm_raster(style = "cont",
            title = "Predicted Prevalence") +
  tm_shape(survey_data) +
  tm_symbols(col = "prev",
             palette = "viridis",
             size = 0.1,
             border.col = "black",
             title.col = "Survey Prevalence",
             style = "cont") +
  tm_layout(
    legend.outside = TRUE,
    legend.outside.position = "right",  # or "bottom"
    legend.title.size = 1.2,
    legend.text.size = 0.8
  )
