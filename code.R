# The microtubules growth model explains the shape and size of the cytoplasmic pattern
# Ling Jin 
# BISC 444 Final Project

rm(list=ls())


setwd("./Desktop/USC/BISC444/Project")
# visualize the ER image ----
library(OpenImageR)
library(EBImage)
ER = readImage('./ER.tif') # Import the .tif image
if (interactive()) 
  display(ER) # Display the ER image
threshold_ER = mean(ER) # Determine the threshold for binarization
Binarize_ER = 1*(ER>threshold_ER) # Binarize
if (interactive()) 
  display(Binarize_ER) # Display the binarized ER image

# clean up image and smooth the contour for each cell
md.im <- medianFilter(Binarize_ER, size =20) # Apply a median filter with window size of 20
display(md.im)

# display(gblur(Binarize_ER, sigma=4)) # Apply a gaussian blur with window size of 20
# contour(md.im) #Show the cotour in a coordinate system

library(dbscan)
library(ggplot2)
library(deldir)
library(ggvoronoi)
library(raster)
library(grid)
library(mmand)

# visualize the Nuclei image ---- 
NLS <- readImage('./mCherryNLS.tif')
if (interactive()) 
  display(NLS)

# binarize the image ---- 
# threshold_NLS = 0.2 #set the threshold manually 
# Binarize_NLS = 1*(NLS > threshold_NLS) 
Binarize_NLS = 1*(NLS < threshold(NLS, method="kmeans")) #set the threshold by k-means
if (interactive()) 
  display(Binarize_NLS) 
img = image(t(Binarize_NLS),col=grey(seq(0, 1))) #another way to display

# visualize the Nuclei image ---- 
coord <- which(Binarize_NLS==1, arr.ind=TRUE)
plot(coord)

# clustered nuclei pixels based on density
clust <- dbscan(coord, eps = 20, minPts = 1) #cluster pixels based on density
# plot all points in the same cluster with the same color
plot(coord[clust$cluster %in% 1:55,], col=clust$cluster[clust$cluster>0]) 

# plot all points in the same cluster with the cluster number
plot(coord[clust$cluster %in% 1:55,],type='n')
text(x=coord[,1], y=coord[,2], labels = clust$cluster)

# Find the center for each cluster ----
nuclei_coord <- c()
for(i in 1:nrow(table(clust$cluster))){
  x_mean = mean(coord[clust$cluster==i,1])
  y_mean = mean(coord[clust$cluster==i,2])
  nuclei_coord = rbind(nuclei_coord,c(x_mean, y_mean))
}

plot(x=nuclei_coord[,1],y=nuclei_coord[,2]) # plot a set of single pixels, representing the center point of each nucleus

# Test out the accuracy in visually checking the overlapping rate between the simulated nuclei and the original nuclei image ----
# The superimposing of the nuclei map with the non-transparent original nuclei image with GFP-NLS labels.
#plot(x=nuclei_coord[,2],y=nuclei_coord[,1])
#rasterImage(NLS_rotate,0,0,2048,2048,col = "transparent")

# The superimposing of the nuclei map with the original nuclei image with GFP-NLS labels.
NLS_rotate <- readImage('./mCherryNLS_rotate.tif')
ggplot(data = as.data.frame(nuclei_coord), 
       aes(x=nuclei_coord[,2],y=nuclei_coord[,1])) + 
  annotation_custom(rasterGrob(NLS_rotate, 
                               width = unit(1,'npc'), 
                               height = unit(1,'npc')), 
                    0, 2048, 0, 2048) +
  geom_point(color = "white",
             size = 5, shape = 1, alpha = 0.5)

# Generate the voronoi diagram ----
a <- voronoi_polygon(as.data.frame(nuclei_coord),x="V1",y="V2")
b <- fortify_voronoi(a)
plot(a)

# Test out the accuracy in visually checking the overlapping rate between the voronoi diagram and the original nuclei image ----
ER_raster <- raster('./ER.tif')
ER_df <- as.data.frame(ER_raster, xy = TRUE)

df <- data.frame(nuclei_coord[,1], nuclei_coord[,2])

ggplot(df, aes(nuclei_coord[,2], nuclei_coord[,1])) +
  stat_voronoi(geom = "path",size=1) +
  geom_point()+
  geom_raster(data = ER_df , aes(x = x, y = y,fill = ER),alpha =0.7,show.legend = FALSE)


