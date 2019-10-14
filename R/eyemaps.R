#   Code for generating a visual map of human cell migratin
#  TODO:
#      - convert to markup
#      - (re)capture all data from 20090PR's and injection sites from 2008PR's
#      - label injection site
#      - plot axis/titles
#      - package the code
#      - write it up? for where? Application note?
#
#  https://www.r-graph-gallery.com/27-levelplot-with-lattice.html
#  https://stackoverflow.com/questions/52239714/center-align-text-of-unequal-lengths-in-rows-in-r-dataframe-and-fill-flanking-po
#
# Process
#  1) convert all rows to single strings
#  2) "center" them by adding characters to either end
#  3) break back into columns
#  4) convert characters to specific values
#  5) levelplot
#
#  Key: 0 = zero human cells () (blue,3)
#       I = contains some number of human cells  (red,4)
#       U = unobservable(or partially) (light blue,2)
#       D = not counted/damaged (gray,1)
#       X = Filled to keep all centered (white, 0)

#Another resource if I wanted to add texture to certain areas
#  https://stackoverflow.com/questions/9415344/using-patterns-in-addition-instead-of-background-colors-in-lattice-plots


# Basic approach: Instead of trying to remap the grid, each (section, image length) coordinate is a point to locate a colored symbol. 
#    0) start with (section, image length) coordinates      
#    1) translate to grid based (x,y) location 
#    2) translate to circularized (x,y) location
#    3) get coordinates of optic nerve
#    4) get coordinates of injection site
#    5) plot points
#    6) draw and shade circle suggesting eye
#
# Eye details:
#    - samples/sections cut horizontally starting at the top of the eye  
#    - sections are only started to be counted once human cells observed
#    - yields details count data for a band across the midline of the eye
#    - generally sections are 5 um thick
#    - generally image lengths are 680 um wide
#
# Start with Background:
#    - determine extent (sections * thickness x imagelengths * image length)
#    - randomly generate points within the extent
#    - circularize the coordinates
#    - plot the random points in a "background" color
#
#
# Add annotation of the injection site.
#    Injection slides are annotated on the 2009, but it doesn't indicate
#    the precise section/image length. However, 2008PR might have it.
#    However, the information I looked at for one example didn't match up.
#
# ??QUESTIONS??:
#    what if the injection site annotation doesn't match between 2008PR and 2009PR?
#    what to do with gaps in the recorded data on 2009PR?
#



#Check a directory
#thefolder <- "/Users/kiehlt/Documents/career/NSCI/GLP/2009/data"

#get a list of files
#allfiles <- list.files(thefolder, pattern="\\.txt$")

rpe.plot.pr.rescue <- function(thefile, title, xlab="X axis", ylab="Y axis"){
	
	#Read in the row data
	
	#
	
	
}
rpe.get.data <- function(thefile, title, xlab="Image Length", ylab="Section") {
	##
	##
	##
	##
	scalefig <- 2000
	sectionthickness <- 5
	sectionsperslide <- 8
	imagelength <- 680
	
	#Read data from the given file
	thedata <-
		read.table(thefile,
			sep = "\t",
			header = TRUE
		)
	thedata[is.na(thedata)] <- "" #get rid of any dangling NA's
	
	#convert from (section,slide) to microns of thickness to determine the values on one axis
	vertical <- sectionthickness * thedata[,1] * sectionsperslide + (thedata[,2]*sectionthickness)

	#concatenate all of the imagelength annotations for each slide
	vals <- do.call(paste0, thedata[, 3:ncol(thedata)])
	#print(vals)
	
	#buffer each end of the section "string" with X's to align all of the sections anchoring them at their centers
	nr <- nchar(vals) #get the number of imagelengths recorded for each section
	mx <- max(nr)     #find the max number of image lengths recorded of all sections
	i1 <- ceiling((mx - nr) / 2) #find the difference between the longest and each of the others, then halve that value
	out <-
		ifelse(i1 > 0, paste0(strrep("X", i1), vals, strrep("X", i1)),
			  vals) #stick X's to the ends
	extended <- substr(out, 1, mx) #trim off any extra X's
	
	## Split out to separate columns
	extsplt <- str_split(extended, pattern = "")
	
	finmat <- do.call("cbind", extsplt)
	
	#  Key:
	
	#       X = Filled to keep all centered (magenta, 0)
	#       D = not counted/damaged (pink,1)
	#       U = unobservable(or partially) (white,2)
	#       0 = zero human cells () (cyan,3)
	#       I = contains some number of human cells  (navy,4)
	#       N = Indicates evidence of the injection site (circled around center of mass, 5)
	#
	#convert values to something we can level plot
	# replace X
	finmat[finmat == "0"] <- "3"
	# replace D
	finmat[finmat == "X"] <- "0"
	# replace U
	finmat[finmat == "D"] <- "1"
	# replace U
	finmat[finmat == "U"] <- "2"
	# replace I
	finmat[finmat == "I"] <- "4"
	# replace N
	finmat[finmat == "N"] <- "5"
	
	#pick our set of colors
	breaks <- (c(-1, 0, 1, 2, 3, 4))
	cols <- colorRampPalette(c("magenta", "plum1", "white", "paleturquoise1", "navy"))
	#return(levelplot(finmat, at = breaks, col.regions = cols, xlab=xlab, ylab=ylab, main=title ))
		
	#find the coordinates of each label
	
	five <- which(finmat=="5", arr.ind=TRUE)
	centroid <- cbind(five[,1], vertical[five[,2]])
	#plotfive <- cbind(five[,1], vertical[five[,2]])
	centroid[,1] <- centroid[,1] * imagelength
	centroid <- cbind(mean(centroid[,1]), mean(centroid[,2]))
	
	four <- which(finmat=="4", arr.ind = TRUE)
	plotfour <- cbind(four[,1],vertical[four[,2]])
	plotfour[,1] <- plotfour[,1] * imagelength
	
	three <- which(finmat=="3", arr.ind = TRUE)
	plotthree <- cbind(three[,1], vertical[three[,2]])
	plotthree[,1] <- plotthree[,1] * imagelength

	two <- which(finmat=="2", arr.ind = TRUE)
	plottwo <- cbind(two[,1], vertical[two[,2]])
	plottwo[,1] <- plottwo[,1] * imagelength
	
	one <- which(finmat=="1", arr.ind = TRUE)
	plotone <- cbind(one[,1], vertical[one[,2]])
	plotone[,1] <- plotone[,1] * imagelength
	
	#determine the range
	maxrange <- (ncol(thedata)-2) * imagelength
	
	#plot the background
	rpe.viz(plotscale=scalefig)

	pme <- ((plotone/maxrange) * 2) -1  #convert values to a unit grid (2x2)
	cmat <- matrix(rpe.circularize(pme[,2], pme[,1]),ncol=2) #transform values to a unit circle (radius=1)
	cmat <- jitter(cmat, amount=0.05) #jitter the points to reduce overlap
	#matplot(cmat[,1], cmat[,2], pch=19, col="navy")
	points(cmat * scalefig,pch=19, cex=0.75, col="plum1")
	
	pme <- ((plottwo/maxrange) * 2) -1
	cmat <- matrix(rpe.circularize(pme[,2], pme[,1]),ncol=2)
	cmat <- jitter(cmat, amount=0.05)
	#matplot(cmat[,1], cmat[,2], pch=19, col="navy")
	points(cmat * scalefig,pch=19, cex=0.75,col="white")
	
	pme <- ((plotthree/maxrange) * 2) -1
	cmat <- matrix(rpe.circularize(pme[,2], pme[,1]),ncol=2)
	cmat <- jitter(cmat, amount=0.05)
	#matplot(cmat[,1], cmat[,2], pch=19, col="navy")
	points(cmat * scalefig,pch=19, cex=0.75,col= rgb(175/255,238/255,238/255, alpha=0.5)) #"paleturquoise1)"
	
	#convert the values to a unit square "diameter" = 2
	pme <- ((plotfour/maxrange) * 2) -1
	cmat <- matrix(rpe.circularize(pme[,2], pme[,1]),ncol=2)
	cmat <- jitter(cmat, amount=0.05)
	#matplot(cmat[,1], cmat[,2], pch=19, col="navy")
	points(cmat * scalefig,pch=19, cex=0.75,col="navy")
	
	# Circle the Injection site
	#convert the values to a unit square "diameter" = 2
	pme <- ((centroid/maxrange) * 2) -1
	cmat <- matrix(rpe.circularize(pme[,2], pme[,1]),ncol=2)
	#matplot(cmat[,1], cmat[,2], pch=19, col="navy")
	points(cmat * scalefig,pch=13, cex=2.5, lwd=2.0, col="yellow")
	
	#shift to center on origin
	
	#circularize
	
	four <- cbind(vertical, four[,1])
	#which(finmat=="")
		
	#return(finmat)
}

#rpe.plot(paste(thefolder, allfiles[1], sep="/"), title=allfiles[1])

# Determine the extent of a prototypical eye
#    - ~3mm diameter  ("3-4")
#    - 5um section thickness
#    - 680 um image length (68?)
#    - 
#
#rpe.extent <- c()

rpe.circularize <- function(ex, ey){
	# After building the grid of data, we have to transform to layout on a circle
	#         http://squircular.blogspot.com/2015/09/mapping-circle-to-square.html      
	#         https://www.xarg.org/2017/07/how-to-map-a-square-to-a-circle/
	#         https://stackoverflow.com/questions/1621831/how-can-i-convert-coordinates-on-a-square-to-coordinates-on-a-circle
	#
	cx <- ex * sqrt(1-((ey * ey)/2))
	cy <- ey * sqrt(1-((ex * ex)/2))
	return(c(cx, cy))
}

rpe.viz <- function(plotscale=1){
	# This function opens a plot of a rough unit circle filled with randomly 
	#     placed gray points. This serves as a background to plot our actual
	#     data against.
	#
	#
	#
	
	#generate random x/y points within a 2x2 grid
	exs <- runif(5500)*2 -1
	eys <- runif(5500)*2 -1
	
	#put those x/y's in a matrix
	eucpoints <- matrix(exs, ncol=2)
	
	#map our 2x2 grid data to points on a unit circle 
	circpoints <- matrix(rpe.circularize(exs, eys), ncol=2)
	
	#scale the points if neccesary
	circpoints <- circpoints * plotscale
	
	#plot the random points
	matplot(circpoints[,1], circpoints[,2], pch=19, col="#BBBBBB", cex=2.5, 
		   ylab="Dorsal(Image Length um)", xlab = "Nasal (Section)K")
	
	#Do this again with a different gray tone to fill in that background
	exs <- runif(2500)*2 -1
	eys <- runif(2500)*2 -1
	eucpoints <- matrix(exs, ncol=2)
	circpoints <- matrix(rpe.circularize(exs, eys), ncol=2)
	#matplot(circpoints[,1], circpoints[,2], pch=19, col="#CCCCCC", cex=1.5)
     points(circpoints * plotscale, pch=19,col="#DDDDDD", cex=1.5)
}

rpe.render <- function(){
	pme <- ((plotfour/4500) * 2) -1
	cmat <- matrix(rpe.circularize(pme[,2], pme[,1]),ncol=2)
	matplot(cmat[,1], cmat[,2], pch=19, col="navy")
	points(cmat,pch=19, col="navy")
	
	#plot random point field as background
	
	#overlay observations
}
