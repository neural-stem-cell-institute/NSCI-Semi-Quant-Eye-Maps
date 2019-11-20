#   Code for generating a visual map of human cell migrating within the eye
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
# Note that the scale for each eye is independent so the maximal value (regardless of how small)
#    will always appear purple. In the case of these two the maximum value is 2, and it appears 
#    that values close to this are distributed across the eye. Note that this effect is even more
#    apparent in some of the visualizations for experiment 36 (such as figures 33, 35,36,37, 38 
#    and 43) where the maximal values range from 1 to 1.5 and similarly appear to be distributed 
#    evenly across the figures. Furthermore, note that our method of counting in these cases does 
#    not distinguish between a single cell and a contiguous line of cells comprising a single layer.
#    Each of these cases would be counted as "1" as the maximal layer depth.
#
# TODOs:
#    - make dots bigger
#    - don't cut off the threshold
#    - Try coloring categorically for photoreceptor rescue
#    - 0.0 is a color, perhaps white 
#    - scale each experiment?
#    - different scale per experiment
#    - potentially change the iteration so we iterate by experiment or eye and can group each eye's images together.



#Check a directory
#thefolder <- "/Users/kiehlt/Documents/career/NSCI/GLP/2009/data"
thefile <- "/Users/kiehlt/Documents/github-projects/RSCC-Eye-Maps/data/pr-rescue/3002-30B1L.txt"
nsci.jitteramount <- 0.05
nsci.globalmax.prrescue <- 7.4 #TODO: add check to make sure actual data isn't more than globally defined max
nsci.useglobalmax <- FALSE
#get a list of files
#allfiles <- list.files(thefolder, pattern="\\.txt$")
rgb2hex <- function(rngnb){
	sprintf('#%s',paste(as.hexmode(c(rngnb[1],rngnb[2],rngnb[3])),collapse = ''))
}

#' Shift values left to right in a vector. This is used when centering the values in a matrix.
#' 
shiftrow <- function(vec, shiftn, ignoreleading=0){
	for (i in length(vec):ignoreleading+1){
		
		if(i>shiftn+ignoreleading){
		    vec[i] <- vec[i-shiftn]
		} else {
			vec[i] <- ''
		}
	}
	
	return(vec)
	
		
}

#' Re
redist <- function(x){
	(x-min(x))/diff(range(x))
	}

valuesToLevels <- function(values, levels)
	
plotLevelPoints <- function(x, y, value){
	# could use colorRampPalette and the bias option
	btoy <- colorRamp(c("blue", "yellow"))
	print(max(x))
}

rpe.plot.pr.rescue <- function(thefile, title="Eye Map", totalslides=100, xlab="X axis", ylab="Y axis", threshold=2.0, usethreshold=FALSE){
	
	#Scale/size constants
	scalefig <- 4000  #changed from 2000, eye has diameter of 4mm
	sectionthickness <- 5
	sectionsperslide <- 8
	imagelength <- 450  #may cut this in 2/3 to handle scaling (from 680)
	dotsize <- 1.25 #size of the plotted points (cex in points)
	
	#Define set of colors for plotting the range
	btoy <- colorRamp(c("white", "yellow", "red", "purple"))
	message("before read.table")
	#Read in the eye data
	eyedata <- read.table(thefile, sep="\t", stringsAsFactors = FALSE)
	message("after read.table")
	#center the sections
	for (rnum in 1:nrow(eyedata)){
		#for each row, count the empty values
		n <- length(which(eyedata[rnum,] ==''))
		if(n>1){ #only shift if we are actually going to shift
			n <- floor(n/2)
			#print(paste("shifting row", rnum, "by", n, "positions."))
			#shift the existing values by half the number of empty values, 
			#     - ignoring the first three columns (eye id, slide, section)
			#     - Assuming empty elements are at the "right end" of the row
			eyedata[rnum,] <- shiftrow(eyedata[rnum,], shiftn=n, ignoreleading = 3) 
			}
		}
	
	#create a list of points to plot
	
	#convert from (section,slide) to microns of thickness to determine the values on one axis
	#    Assumes that the rows are slide, section, observation, observation .... etc.
	#    Ignores the first row, which are labels (eyedata[-1,2])
	collabels <- eyedata[1,] #retain the column labels in case we need them
	rowidentifiers <- eyedata[, c(1,2,3)] #retain the row labels for later
	vertical <- sectionthickness * as.numeric(eyedata[-1,2]) * sectionsperslide + (as.numeric(eyedata[-1,3])*sectionthickness)

	#adjust vertical to account for "uncounted" sections at top and bottom of eye
	vertical <- vertical + (scalefig - totalslides*8*5)/2
	
	eyedata <- eyedata[-1,c(-1,-2,-3)] #drop the row with the column labels and the columns with row identifiers
	
	imagelength <- scalefig/ncol(eyedata) * 1.1 #scale the imagelength based on the number of columns in eyedata
	
	#find locations of injection site (ignoring the first three columns)
	is <- which(eyedata[,]=="N", arr.ind=TRUE)
	#isloc <- c(mean(is[,1]), mean(is[,2]))
			 
	#find locations of optic nerve
	on <- which(eyedata[,]=="O", arr.ind=TRUE)
	#onloc <- c(mean(on[,1]), mean(on[,2]))
	
	#which rows have either injection site or optic nerve
	isonrows <- unique(c(is[,1], on[,1]))

	#remove rows with injection site and optic nerve
	eyedata <- eyedata[-isonrows,] 

	#Convert to array x,y values (get the locations for each point)
	points <- which(eyedata!='', arr.ind=TRUE)

	#get values and scale to 0..1 range
	values <- as.numeric(eyedata[points])
	
	#capture the actual maximum value before normalization and before thresholding
	maxval <- max(values) #lastedit
	
	#zero out values below threshold
	if(usethreshold){ #lastedit
		values[which(values<threshold)] <- 0.0  # NEW
	}
	
	
	#count the number below threshold, if it's the same as the total number than don't do this.
	#  Can't do this calculation if all the values are zero
message("here")
	
	if(usethreshold){   #lastedit
	if(length(which(values<threshold)) < length(values)){
		 if(nsci.useglobalmax){
		 	values <- ifelse(values >0, (values-threshold)/(nsci.globalmax.prrescue-threshold), values) #subtracting out threshold so the color range is scaled nicely starting at the threshold value
		 	# values <- values/nsci.globalmax.prrescue.threshold
		 }else{
		    values <- ifelse(values >0, (values-threshold)/(max(values-threshold)))#subtracting out threshold so the color range is scaled nicely starting at the threshold value
		 	# values <- values/max(values) 
		 }
	}}else{
	values<- values/maxval  #lastedit
}
	message(paste0("here..." , min(values), ", ", max(values)))
	#Convert x,y to grid coordinates (vertical still contains rows with is and on)
	points[,1] <- vertical[-isonrows][points[,1]]  #convert x to eye coords based on section thicknesses
     points[,2] <- points[,2]*imagelength #convert y grid to eye coords based on image lengths

     #TODO: figure out coords for injection site and optic nerve
     isloc <- c(mean(vertical[is[,1]]), mean(is[,2])*imagelength)
     onloc <- c(mean(vertical[on[,1]]), mean(on[,2])*imagelength)		 
     #isloc <- c(mean(is[,1]), mean(is[,2]))
     #onloc <- c(mean(on[,1]), mean(on[,2]))

     #convert points to unit grid
     
     #scale section dimension (dorsal/ventral) based on total sections or known proportion
     maxdim.section <- max(scalefig, max(points[,1])) # (TODO: determine specifics for each eye)This scaling max defines the spread of the sections in the final vizualization
     #scale image length dimension (nasal/temporal) based on the longest section, with the most image lengths
     maxdim.imglngth <- max(scalefig, max(points[,2]))
     
     points[,1] <- ((points[,1]/maxdim.section) *2) -1
	points[,2] <- ((points[,2]/maxdim.imglngth) *2) -1
	
	isloc[1] <- ((isloc[1]/maxdim.section) *2) -1
	isloc[2] <- ((isloc[2]/maxdim.imglngth) *2) -1
	onloc[1] <- ((onloc[1]/maxdim.section) *2) -1
	onloc[2] <- ((onloc[2]/maxdim.imglngth) *2) -1
	
	#jitter the points
	points <- jitter(points, amount=nsci.jitteramount)

	#Circularize the points
	#    negative signs rotate to proper orientation for visualization
	pmat  <- matrix(rpe.circularize(-points[,2], -points[,1]), ncol=2) #swapped 2 and 1 here

	ismat <- matrix(rpe.circularize(-isloc[2], -isloc[1]), ncol=2)

	onmat <- matrix(rpe.circularize(-onloc[2], -onloc[1]), ncol=2)



	

	
	#Determine colors for each point based on previously defined colorRamp
	pcolors <- apply(round(btoy(values)), 1, rgb2hex)
	
	#color cells below threshold - didn't work at this stage
	#    since the values in values had been normalized
	#pcolors[which(values<threshold)] <- "#ffffff"
	
	ppallette <- colorRampPalette(c("magenta", "plum1", "white", "paleturquoise1", "navy"))

	#plot the background
	rpe.viz(plotscale=1 * 2000)

	
	#plot the points
	points(pmat* 2000,pch=19, cex=dotsize, col=pcolors)
	
	#make a legend
	rpe.addLegend(maxval = maxval)
	
	#Plot the Injection Site
	points(ismat* 2000,pch=13, cex=2.5, lwd=2.0, col="blue")
	
	#Plot the optic nerve
	points(onmat* 2000,pch=13, cex=2.5, lwd=2.0, col="green")

     cols <- colorRampPalette(c("magenta", "plum1", "white", "paleturquoise1", "navy"))	
     breaks <- (c(-1, 0, 1, 2, 3, 4))
     levels(factor(eyedata[eyedata!='']))
     
     #associate colors with the levels
	return(eyedata)
	
}

rpe.get.data <- function(thefile, title, totalslides=100, xlab="Image Length", ylab="Section") {
	## Plots HLA Tracking data
	##
	##
	##
	scalefig <- 4000
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
	#    Assumes that the rows are: slide, section, observation, observation .... etc.
	vertical <- sectionthickness * as.numeric(thedata[,1]) * sectionsperslide + (as.numeric(thedata[,2])*sectionthickness)
 
	#adjust vertical to account for "uncounted" sections at top and bottom of eye ADDED
	vertical <- vertical + (scalefig - totalslides* 8*5)/2
	
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
	
	finmat <- do.call("cbind", extsplt) #finmat each column is a section, rows are imagelengths
     
	#datmat <- do.call("rbind", extsplt) #extmat, each row is a section  10/25
	
	## Update: make sure we're scaled properly in the imagelength direction	
	imagelength <- scalefig/nrow(finmat) #10/25
	#imagelength <- scalefig/ncol(datmat) #10/25
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
	# replace O
	finmat[finmat == "O"] <- "6"
	

	
	#pick our set of colors
	breaks <- (c(-1, 0, 1, 2, 3, 4))
	cols <- colorRampPalette(c("magenta", "plum1", "white", "paleturquoise1", "navy"))
	#return(levelplot(finmat, at = breaks, col.regions = cols, xlab=xlab, ylab=ylab, main=title ))
		
	#find the coordinates of each label
	
	five <- which(finmat=="5", arr.ind=TRUE) #the columns of 5 are coordinates in (imagelength)
	centroid <- cbind(five[,1], vertical[five[,2]])
	centroid[,1] <- centroid[,1] * imagelength
	centroid <- cbind(mean(centroid[,1]), mean(centroid[,2])) #find injection site centroid
	
	six <- which(finmat=="6", arr.ind=TRUE) #the columns of 6 are coordinates in (imagelength)
	oncentroid <- cbind(six[,1], vertical[six[,2]]) 
	oncentroid[,1] <- oncentroid[,1] * imagelength
	oncentroid <- cbind(mean(oncentroid[,1]), mean(oncentroid[,2])) #find optic nerve centroid
	
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
	rpe.viz(plotscale=2000) #todo: change to -1,+1 range
      
	#Note: in the below circularization steps the negative
	#    signs rotate to proper orientation for visualization
	
	#TODO: confirm that section and imagelength dimensions are being normalized to respective max values
	#scale section dimension (dorsal/ventral) based on total sections or known proportion
	#maxdim.section <- max(scalefig, max(points[,1])) # (TODO: determine specifics for each eye)This scaling max defines the spread of the sections in the final vizualization
	#scale image length dimension (nasal/temporal) based on the longest section, with the most image lengths
	#maxdim.imglngth <- max(scalefig, max(points[,2]))
	
	pme <- ((plotone/maxrange) * 2) -1  #convert values to a unit grid (2x2)
	#pme[,c(1,2)] <- pme[,c(2,1)] #reorient TODO: double check this
	
	cmat <- matrix(rpe.circularize(-pme[,1], -pme[,2]),ncol=2) #transform values to a unit circle (radius=1)
	cmat <- jitter(cmat, amount=nsci.jitteramount) #jitter the points to reduce overlap
	#matplot(cmat[,1], cmat[,2], pch=19, col="navy")
	points(cmat * scalefig/2,pch=19, cex=0.75, col="plum1")
	
	pme <- ((plottwo/maxrange) * 2) -1
	cmat <- matrix(rpe.circularize(-pme[,1], -pme[,2]),ncol=2)
	cmat <- jitter(cmat, amount=nsci.jitteramount)
	#matplot(cmat[,1], cmat[,2], pch=19, col="navy")
	points(cmat * scalefig/2,pch=19, cex=0.75,col="white")
	
	pme <- ((plotthree/maxrange) * 2) -1
	cmat <- matrix(rpe.circularize(-pme[,1], -pme[,2]),ncol=2)
	cmat <- jitter(cmat, amount=nsci.jitteramount)
	#matplot(cmat[,1], cmat[,2], pch=19, col="navy")
	points(cmat * scalefig/2,pch=19, cex=0.75,col= rgb(175/255,238/255,238/255, alpha=0.5)) #"paleturquoise1)"
	
	#convert the values to a unit square "diameter" = 2
	pme <- ((plotfour/maxrange) * 2) -1
	cmat <- matrix(rpe.circularize(-pme[,1], -pme[,2]),ncol=2)
	cmat <- jitter(cmat, amount=nsci.jitteramount)
	#matplot(cmat[,1], cmat[,2], pch=19, col="navy")
	points(cmat * scalefig/2,pch=19, cex=0.75,col="navy")
	
	# Circle the Injection site
	#convert the values to a unit square "diameter" = 2
	pme <- ((centroid/maxrange) * 2) -1
	cmat <- matrix(rpe.circularize(-pme[,1], -pme[,2]),ncol=2)
	points(cmat * scalefig/2,pch=13, cex=2.5, lwd=2.0, col="yellow")
	
	# Circle the Injection site
	#convert the values to a unit square "diameter" = 2
	pme <- ((oncentroid/maxrange) * 2) -1
	cmat <- matrix(rpe.circularize(-pme[,1], -pme[,2]),ncol=2)
	points(cmat * scalefig/2,pch=13, cex=2.5, lwd=2.0, col="green")
	
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
	# Assumes the data is already layed out on a unit square
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
	matplot(circpoints[,1], circpoints[,2], pch=19, col="#D3D3D3", cex=2.5,
		   xlab="", ylab="")
#		   ylab="Dorsal(Image Length um)", xlab = "Nasal (Section)K")
	mtext(text="Dorsal", side=3, line=0, outer=FALSE)
	mtext(text="Ventral", side=1, line=2, outer=FALSE)
	mtext(text="Nasal", side=2, line=2, outer=FALSE)
	mtext(text="Temporal", side=4, line=0, outer=FALSE)
	
	
	#Do this again with a different gray tone to fill in that background
	exs <- runif(2500)*2 -1
	eys <- runif(2500)*2 -1
	eucpoints <- matrix(exs, ncol=2)
	circpoints <- matrix(rpe.circularize(exs, eys), ncol=2)
	#matplot(circpoints[,1], circpoints[,2], pch=19, col="#CCCCCC", cex=1.5)
     points(circpoints * plotscale, pch=19,col="#E7E7E7", cex=1.5)
}

rpe.render <- function(){
	pme <- ((plotfour/4500) * 2) -1
	cmat <- matrix(rpe.circularize(pme[,2], pme[,1]),ncol=2)
	matplot(cmat[,1], cmat[,2], pch=19, col="navy")
	points(cmat,pch=19, col="navy")
	
	#plot random point field as background
	
	#overlay observations
}

rpe.addLegend <- function(maxval=4){
	if(nsci.useglobalmax){maxval <- nsci.globalmax.prrescue}
	maxval <- round(maxval, digits=2)	
	lgd_ = rep(NA, 5)
#	lgd_[c(1,10)] = c(maxval, "< threshold")
		lgd_[c(1,10)] = c(maxval, 0)
	legend(x = -2000, y = -1500,
		  legend = lgd_,
		  fill = colorRampPalette(colors = c( 'purple','red', 'yellow','white'))(10),
		  border = NA,
		  bg="white",
		  y.intersp = 0.2,
		  cex = 0.85, text.font = 1.3)
	
}

