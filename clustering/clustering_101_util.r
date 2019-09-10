# Data Clustering and Segmentation Analysis
# k-medoids in Action - Part 3 Utilities
#
# Author: Jacob Cybulski
# Version: V1.2, Aug 2019
#
# Problem: Provide utility functions for clustering and their visualisation. 
# Method: Deployment, Feature Engineering, k-Medoids, PCA, Visualisation

# Use the relevant packages
require(ggplot2) # ggplot
require(cluster) # pam, silhouette
require(factoextra) # fviz_nbclust, fviz_pca_ind, fviz_pca_var, fviz_add

# Preserve plot geometry
std.width <- getOption("repr.plot.width")
std.height <- getOption("repr.plot.height")
std.par <- par()

# Change the plot geometry
set_plot_dimensions <- function(width_pc, height_pc, margin=c(2,2,2,2)) {
    par(mar = margin)
    options(repr.plot.width=std.width * width_pc, repr.plot.height = std.height * height_pc)}

# Plots data in clus clusters, using PCA and PAM cluster models
plot.pam.in.pca <- function(pca, pam.model, data, clus, comp1=1, comp2=2, title="PCA Chart", 
                            ellipse=TRUE, xlim=NULL, ylim=NULL, col=NULL, alpha=0.1) {
    # Find cluster details
    k <- nrow(pam.model$clusinfo)
    medoids.in.pca <- data.frame(scale(pam.model$medoids, pca$center, pca$scale) %*% pca$rotation)  

    # Rotate medoids to PCA coordinates
    PoV <- round(100*pca$sdev^2/sum(pca$sdev^2),digits=2)

    # Set default colour palette
    #colfunc <- colorRampPalette(c("yellow", "black"))
    #colfunc <- colorRampPalette(c("navy", "red", "orange"))
    #colfunc <- colorRampPalette(c("deepskyblue3", "royalblue4", "red3", "red1", "yellow", "green"))
    colfunc <- colorRampPalette(c("limegreen", "royalblue", "cyan3", "magenta3", "gold3", "red2"))
    if (!is.null(col)) colfunc <- colorRampPalette(col)
    cols <- colfunc(k)

    # Plot new clustered data according to the PCA
    pca.plot <- ggplot(data, aes(x=data[, comp1], y=data[, comp2])) + 
        ggtitle(title) +
        labs(x=paste("PC", comp1, " (", PoV[comp1], "%)", sep=""), 
             y=paste("PC", comp2, " (", PoV[comp2], "%)", sep="")) +
        geom_point(aes(color=cols[clus]), size=0.5, pch=20) +
        scale_colour_identity(guide="legend", breaks=cols, labels = as.factor(1:k), 
                              aesthetics = c("colour", "fill"))

    # Add legend
    pca.plot <- pca.plot +
        guides(colour = guide_legend(title="Clusters", override.aes = list(size = 5)))

    # Add elipse
    if (ellipse) pca.plot <- pca.plot +     
        stat_ellipse(aes(color=cols[clus], fill=cols[clus]),
            geom = "polygon", type = "t", linetype = 1, level=0.6, alpha = alpha, show.legend=FALSE)

    # Add medoid points
    pca.plot <- pca.plot +
        geom_point(data=medoids.in.pca, aes(x=medoids.in.pca[, comp1], y=medoids.in.pca[, comp2]),
             color="black", size=5, pch=10) +
        geom_text(data=medoids.in.pca, aes(x=medoids.in.pca[, comp1], y=medoids.in.pca[, comp2], label=paste(1:k)),
              fontface="bold", hjust=-0.7, vjust=-0.5)

    # Add axes limits
    if (!is.null(xlim)) pca.plot <- pca.plot + xlim(xlim[1],xlim[2])
    if (!is.null(ylim)) pca.plot <- pca.plot + ylim(ylim[1],ylim[2])   

    # Plot
    return(pca.plot)
}

# Plots data in clus clusters, using PCA and PAM cluster models
# Note that medoids are already in PCA coordinates
plot.clus.rpca <- function(rpca.model, clus.model, data, cluster=NULL, comp1=1, comp2=2, col=c("blue", "red"), 
                           title="PCA Chart", xlim=NULL, ylim=NULL, alpha=0.1, level=0.8) {
    col.ramp <- colorRampPalette(col)
    cols <- col.ramp(kNo)
    medoids <- data.frame(clus.model$medoids)
    if (is.null(cluster)) cluster <- data.frame(Cluster=clus.model$clustering)$Cluster
    PoV <- round(100*rpca.model$sdev^2/sum(rpca.model$sdev^2),digits=2)
    k <- nrow(medoids)
    g <- ggplot(data=data, aes(x=data[, comp1], y=data[, comp2])) +
        ggtitle(title) +
        geom_point(aes(color=cols[cluster]), size=0.3, pch=20) +
        labs(x=paste("PC", comp1, " (", PoV[comp1], "%)", sep=""), 
             y=paste("PC", comp2, " (", PoV[comp2], "%)", sep="")) +
        scale_colour_identity(guide="legend", breaks=cols, labels = as.factor(1:k), 
                              aesthetics = c("colour", "fill")) +
        guides(colour = guide_legend(title="Clusters", override.aes = list(size = 5))) +
        stat_ellipse(aes(color=cols[cluster], fill=cols[cluster]),
                geom = "polygon", type = "t", linetype = 1, level=level, alpha=alpha, show.legend=FALSE) +
        # xlim(-6, 6) + ylim(-5, 5) +
        geom_point(data=medoids, aes(x=medoids[, comp1], y=medoids[, comp2]),
             color="black", size=5, pch=10) +
        geom_text(data=medoids, aes(x=medoids[, comp1], y=medoids[, comp2], label=paste(1:k)),
              fontface="bold", hjust=-0.7, vjust=-0.5)

    # Add axes limits
    if (!is.null(xlim)) g <- g + xlim(xlim[1],xlim[2])
    if (!is.null(ylim)) g <- g + ylim(ylim[1],ylim[2])   

    return(g)
}


# Plots outliers in PC coordinates
plot.outl.in.pca <- function(pca, data, outl, comp1=1, comp2=2, title="Outlier Boundary", 
        xlim=NULL, ylim=NULL, col=c("royalblue", "red2"), alpha=0.3, level=0.8, size=1.5,
        ellipse=TRUE, cluster=NULL, add=NULL) {
    colfunc <- colorRampPalette(col)
    cols <- colfunc(2)
    
    ndata <- data[outl==1,]
    set_plot_dimensions(1, 0.7)
    if (!is.null(add)) {
        g <- add
    } else {
        g <- ggplot(data=data, aes(x=data[, comp1], y=data[, comp2])) +
            ggtitle(title)
    }
    g <- g + geom_point(aes(color=cols[outl]), size=size, pch=20) +
        labs(x=paste("PC", comp1, sep=""), y=paste("PC", comp2, sep="")) +
        scale_colour_identity(guide="legend", breaks=cols, labels = c("Within", "Outside"), 
                              aesthetics = c("colour", "fill")) +
        guides(colour = guide_legend(title="Boundary", override.aes = list(size = 5)))
    if (ellipse) g <- g + stat_ellipse(type = "t", linetype = 2, col="red", level=level)

    # Add axes limits
    if (!is.null(xlim)) g <- g + xlim(xlim[1],xlim[2])
    if (!is.null(ylim)) g <- g + ylim(ylim[1],ylim[2])   

    return(g)
}

    
