# clean environment
rm(list=ls())

#Loading Packages 
easypackages::packages ("sf", "spdep", "tmap", "spatialreg", 'spgwr',
                        "mapview", "RColorBrewer", "tidyverse",
                        "cowplot", "leafsync", "leaflet.extras2", "ape",
                        "xtable")

#Working directory
setwd('~/Documents/Spatial/Term project/Datasets')

#Read data
prices <- read.csv('HousingPrices-Amsterdam-August-2021.csv')
crime <- read.csv('CrimeSectionsClean.csv', sep = ';')
age <- read.csv('age.csv', sep=';')

#Clean age data
age$Population <- age$people / sum(age$people)
age$AvgAge <- age$average
age  <- subset(age, select = -c(average, people)) 

#Remove missing/outliers
prices <- na.omit(prices) # Four missing price
prices <- prices[!prices$X == "213" & !prices$X == '628', ] # Two in Haarlem, not Amsterdam
crime <- subset(crime, select = -c(X, X12.17.jaar, X18.24.jaar, X25.plus)) #Just take totals

#Generate prices per square meter
prices$PriceSQM <- prices$Price/prices$Area

#Geo data
sections <- read.csv('CitySectionsClean.csv', sep = ';')
sections.sf <- st_as_sf(sections, wkt = 'WKT_LNG_LAT', coords = c('LNG', 'LAT') ) 
prices.sf <- st_as_sf(prices, coords = c('Lon', 'Lat'))

#Spatial join houses into city sections
df <- st_join(sections.sf, prices.sf)
df <- subset(df, select = -c(Oppervlakte_m2, Address, Zip, Room) )

#Mapping data
pricemap <- mapview(prices, xcol = "Lon", ycol = "Lat", crs=4326, zcol = 'PriceSQM',
                    at = seq(2000, 12000, 1000), layer.name = 'Price per Square Meter') + 
  mapview(sections.sf, layer.name = 'City Sections')

pricemap

#Aggregate data
df_sections <- aggregate(Price ~ Stadsdeelnaam + LNG + LAT, data=df, mean)
df_sections <- merge(df_sections, crime, by = 'Stadsdeelnaam')
df_sections <- merge(df_sections, age, by = 'Stadsdeelnaam')

#Prices per section
price_df <- merge(df, df_sections, by = 'Stadsdeelnaam' )
price_df <- merge(prices, price_df, by = 'X' )
price_df <- subset(price_df, select = -c(LNG.x, LAT.x, Price.x, Area.y, LNG.y, LAT.y, Price.y))

#Clean dataframe
price_df <- subset(price_df, select = c(X, Price, Area.x, Room, Total, AvgAge, Population, Lon, Lat))

#Linear regression
price_eq <- Price ~ Total + AvgAge + Population + Area.x + Room
price_reg <- lm(price_eq, data = price_df)

#Examining linear regression
vif(price_reg)
summary(price_reg)

AIC(price_reg)
price_reg

#Making matrix of price distributions
price.dists <- as.matrix(dist(cbind(price_df$Lon, price_df$Lat)))
price.dists.inv <- 1/price.dists
diag(price.dists.inv) <- 0

#Making sure no duplicated coords
pricecoord <- coordinates(cbind(price_df$Lon,price_df$Lat))
pricecoord = pricecoord + rnorm(nrow(pricecoord), mean = 0, sd = 0.01)
nb<-tri2nb(pricecoord, row.names = NULL)

#Morans I
moran.test(price_df$Price, nb2listw(nb))

mc_global <- moran.mc(price_reg$residuals, nb2listw(nb), 2999, alternative="greater")
plot(mc_global)

#Spatial lag model
prices_nbq_w = nb2listw(nb) # weights
spa_lagmodel = lagsarlm(price_eq, data = price_df, listw= prices_nbq_w, zero.policy = TRUE)

summary(spa_lagmodel)

#Error model
spa_errmodel <- errorsarlm(price_eq, data=price_df, listw= prices_nbq_w, tol.solve=1.0e-30)
summary(spa_errmodel)

#GWR coord binding
gwrcoord <- cbind(price_df$Lon, price_df$Lat)

#Finding fixed bandwidth

'''
fbw <- gwr.sel(price_eq, 
               data = price_df,
               longlat = TRUE,
               method = 'aic',
               coords = gwrcoord,
               gweight = gwr.Gauss, 
               verbose = T)
'''

#Saved bandwidth for speed rerunning code
fbw <- 1.297611

#GWR
gwr_prices = gwr(price_eq,
                data = price_df,
                longlat = TRUE,
                bandwidth = fbw,
                gweight = gwr.Gauss,               
                coords = gwrcoord,
                hatmatrix=TRUE, 
                se.fit=TRUE)

#GWR results
gwr_prices
results <- as.data.frame(gwr_prices$SDF)
results_df <- merge(price_df, results, by.x = c('Lon', 'Lat'), by.y = c('coord.x', 'coord.y'))

#GWR categorise crime
results$crime = results$Total / results$Total_se
results$crime_cat <- cut(results$crime,
                                   breaks=c(min(results$crime), -1.96, 1.96, max(results$crime)),
                                   labels=c("sig","nonsig", "sig"))

#Population
results$pop = results$Population / results$Population_se
results$pop_cat <- cut(results$Population,
                         breaks=c(min(results$Population), -1.96, 1.96, max(results$Population)),
                         labels=c("sig","nonsig", "sig"))

#Age
results$age = results$AvgAge / results$AvgAge_se
results$age_cat <- cut(results$AvgAge,
                         breaks=c(min(results$AvgAge), -1.96, 1.96, max(results$AvgAge)),
                         labels=c("sig","nonsig", "sig"))

#Room
results$rooms = results$Room / results$Room_se
results$rooms_cat <- cut(results$Room,
                         breaks=c(min(results$Room), -1.96, 1.96, max(results$Room)),
                         labels=c("sig","nonsig", "sig"))

#Area
results$area = results$Area.x / results$Area.x_se
results$area_cat <- cut(results$Area.x,
                         breaks=c(min(results$Area.x), -1.96, 1.96, max(results$Area.x)),
                         labels=c("sig","nonsig", "sig"))

#Mapping GWR significance
areamap <- mapview(results, xcol = "coord.x", ycol = "coord.y", crs=4326, zcol = 'area_cat',
                      layer.name = 'Area') + 
  mapview(sections.sf, layer.name = 'City Sections')

crimemap <- mapview(results, xcol = "coord.x", ycol = "coord.y", crs=4326, zcol = 'crime_cat',
                   layer.name = 'Crime') + 
  mapview(sections.sf, layer.name = 'City Sections')

popmap <- mapview(results, xcol = "coord.x", ycol = "coord.y", crs=4326, zcol = 'pop_cat',
                   layer.name = 'Population') + 
  mapview(sections.sf, layer.name = 'City Sections')

agemap <- mapview(results, xcol = "coord.x", ycol = "coord.y", crs=4326, zcol = 'age_cat',
                   layer.name = 'Age') + 
  mapview(sections.sf, layer.name = 'City Sections')

roommap <- mapview(results, xcol = "coord.x", ycol = "coord.y", crs=4326, zcol = 'rooms_cat',
                   layer.name = 'Rooms') + 
  mapview(sections.sf, layer.name = 'City Sections')

resultsmap <- mapview(results_df, xcol = "Lon", ycol = "Lat", crs=4326, zcol = 'localR2',
                      layer.name = 'GWR Local R²') + 
  mapview(sections.sf, layer.name = 'City Sections')

#Maps
crimemap 
popmap
agemap
areamap
roommap

#GWR local r squared
resultsmap
