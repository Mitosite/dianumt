args=commandArgs(T)

#reference bed input
ref_bed <- read.table('../../data/ref.numts.for.plotting.bed', header=FALSE)

#job bed input
job_bed <- read.table(args[1], header=FALSE)

#step1 coverage file input
s1_cov <- read.table(args[2], header=FALSE)

#step 2/3 coverage file input
s2_cov <- read.table(args[3], header=FALSE)

#x and y values for step 1
s1x0 <- s1_cov$V2
s1x1 <- 2:(length(s1x0)+1)
s1y1 <- s1_cov$V3
s1ymax <- max(s1y1)

#x and y values for step 1
s2x0 <- s2_cov$V2
s2x1 <- 2:(length(s2x0)+1)
s2y1 <- s2_cov$V3
s2ymax <- max(s2y1)

#List of chromosomes in user data
chrom <- unique(job_bed$V1)  #Chromosome number

#Produce titles for graph
library(stringr)

chrom_noM <- c()
graph_titles <- c()

for (i in unique(job_bed$V1)) {
if (i!="chrM") {
			#Eliminate all non-
			if (str_sub(i, 1, 4)!="chrG") { 
			chrom_noM <- c(chrom_noM, i)
			num <- str_sub(i, 4, -1)
			title <- paste("Nuclear chromosome", num, sep=" ")
			graph_titles <- c(graph_titles, title)
		}			
	}
}

library(GenomicRanges)

#Create reference genomic range
refgr <- lapply(split(ref_bed, ref_bed$V4), function(i){
    GRanges(seqnames = i$V1, ranges = IRanges(start = i$V2, 
    											end = i$V3, 
    											names = i$V4))
})

#Create job genomic range
jobgr <- lapply(split(job_bed, job_bed$V4), function(i){
    GRanges(seqnames = i$V1, ranges = IRanges(start = i$V2, 
    											end = i$V3, 
    											names = i$V4))
})

print(jobgr)
names <- names(jobgr)

sets <- c()

for (i in 1:length(names)) {
	setname <- paste("jobgr$", names[i], sep="")
	sets <- c(sets, setname)
}

#Consistent graphing parameters
#Size of data panel label
tracklabelsize <- 1.8
#Size of axis label digits
axislabelsize <- 1.5
#Number of ticks on axis
axisticks <- 5

library(karyoploteR)

#########################################
### FOR MITOCHONDRIAL CHROMOSOME PLOT ###
#########################################

mtname <- paste("chrMmitograph.png", sep="") #constructs path to correct directory
png(mtname, width=3000, height=1500, units="px")

mt<-plotKaryotype(labels.plotter = NULL, chromosomes=c('chrM'), plot.type=2)

kpAddMainTitle(mt, main='Mitochondrial chromosome', cex=4)

#Add graduated markers to ideogram
kpAddBaseNumbers(mt, 
				tick.dist = 1000,
				tick.len = 10, 
				tick.col="red", 
				cex=1.5,
                minor.tick.dist = 100, 
                minor.tick.len = 5, 
                minor.tick.col = "gray",
                )

#UPPER BOTTOM PANEL
kpDataBackground(
	mt,
	data.panel=1,
	r0=0.0, 
	r1=0.45,
	col='#d8fdff',     #Pale blue
	)

kpAxis(
	mt,
	data.panel=1,
	r0=0.0, 
	r1=0.45,
	cex=axislabelsize,
	ymax=s2ymax,
	numticks=axisticks,
	) 

kpAddLabels(mt, 
			data.panel=1, 
			labels='Step 2 coverage', 
			r0=0.0, 
			r1=0.45, 
			cex=tracklabelsize, 
			srt=90, 
			pos=3,
			label.margin=0.05,
			)

kpBars(
	mt,
	chr="chrM", 
	x0=s2x0, 
	x1=s2x1, 
	y1=s2y1, 
	ymax=s2ymax, 
	col="#0278ff", 			#Medium blue
	border=NA,
	r0=0.0, 
	r1=0.45,
	data.panel=1,  
	)

#UPPER TOP PANEL
#Panel background
kpDataBackground(mt,
				data.panel=1,
				r0=0.55, 
				r1=1.00,
				col='#d8fdff',     #Pale blue
				)
kpAxis(mt,
		data.panel=1,
		r0=0.55, 
		r1=1.00,
		cex=axislabelsize,
		ymax=s1ymax,
		numticks=axisticks,
		) 

kpAddLabels(mt, 
			data.panel=1, 
			labels='Step 1 coverage', 
			r0=0.55, 
			r1=1.0, 
			cex=tracklabelsize, 
			srt=90, 
			pos=3,
			label.margin=0.05,
			)

kpBars(
	mt,
	chr="chrM",
	x0=s1x0, 
	x1=s1x1, 
	y1=s1y1, 
	ymax=s1ymax,  
	col="#ff6302", 			#Orange
	border=NA,
	r0=0.55, 
	r1=1.0, 
	data.panel=1,  
	) 

#LOWER TOP PANEL
#Panel background
kpDataBackground(mt,
				data.panel=2,
				r0=0.0, r1=0.45,
				col='#d8fdff',     #Pale blue
				) 

kpAxis(mt,
		data.panel=2,
		data=jobgr$mitonumts,
		r0=0.0, 
		r1=0.45,
		cex=axislabelsize,
		numticks=axisticks,
		)

kpAddLabels(mt, 
			data.panel=2, 
			labels='Found numt coverage', 
			r0=0.45, 
			r1=0.0, 
			cex=tracklabelsize, 
			srt=90, 
			pos=3,
			label.margin=0.05,
			) 

#plot mitonumts coverage
kpPlotCoverage(mt,
	data.panel=2, 
	data=jobgr$mitonumts, 
	r0=0.0, 
	r1=0.45, 
	col='#ffe102',				#Yellow
	)

#LOWER BOTTOM PANEL
#Panel background
kpDataBackground(mt,
				data.panel=2,
				r0=0.55, 
				r1=1.00,
				col='#d8fdff',     #Pale blue
				) 

kpAxis(mt,
		data.panel=2,
		data=refgr$reference,
		r0=0.55, 
		r1=1.00,
		numticks=axisticks,
		cex=axislabelsize
		)

kpAddLabels(mt, 
			data.panel=2, 
			labels='Ref numt coverage', 
			r0=0.55, 
			r1=1.00,
			cex=tracklabelsize, 
			srt=90, 
			pos=3,
			label.margin=0.05,
			) 
#plot reference numt coverage
kpPlotCoverage(
	mt,
	data.panel=2, 
	data=refgr$reference, 
	r0=0.55, 
	r1=1.00,
	col='#c402ff',				#Bright purple
	)

dev.off()

####################################
### FOR NUCLEAR CHROMOSOME PLOTS ###
####################################

for (i in 1:length(chrom_noM)){
	nucname <- paste(chrom_noM[i], "mitograph.png", sep="") #constructs path to correct directory
	png(nucname, width=3000, height=1500, units="px")

	nuc<-plotKaryotype(labels.plotter = NULL, chromosomes=c(chrom_noM[i]), plot.type=2)

	kpAddMainTitle(nuc, main=graph_titles[i], cex=4)

	#Add graduated markers to ideogram
	kpAddBaseNumbers(nuc, 
					tick.dist = 10000000,
					tick.len = 10, 
					tick.col="red", 
					cex=1.5,
	                minor.tick.dist = 1000000, 
	                minor.tick.len = 5, 
	                minor.tick.col = "gray",
	                )



	#UPPER BOTTOM PANEL
	#Panel background
	kpDataBackground(nuc,
					data.panel=1,
					r0=0.0, r1=0.45,
					col='#d8fdff',     #Pale blue
					) 

	kpAxis(
		nuc,
		data.panel=1,
		data=jobgr$nuclearnumts,
		r0=0.0, 
		r1=0.45,
		cex=axislabelsize,
		numticks=axisticks,
		)

	kpAddLabels(
		nuc, 
		data.panel=1, 
		labels='Numt read density', 
		r0=0.45, 
		r1=0.0, 
		cex=tracklabelsize, 
		srt=90, 
		pos=3,
		label.margin=0.05,
		)

	kpPlotDensity(
		nuc,
		data.panel=1, 
		r0=0.0, 
		r1=0.45,
		data=jobgr$numtsreads,
		col='#ff6302',					#Orange
		)

	#UPPER TOP PANEL
	#Panel background
	kpDataBackground(nuc,
					data.panel=1,
					r0=0.55, 
					r1=1.00,
					col='#d8fdff',     #Pale blue
					)
	kpAxis(nuc,
			data.panel=1,
			r0=0.55, 
			r1=1.00,
			cex=axislabelsize,
			numticks=axisticks,
			) 

	kpAddLabels(
		nuc, 
		data.panel=1, 
		labels='Contamination density', 
		r0=0.55, 
		r1=1.0, 
		cex=tracklabelsize, 
		srt=90, 
		pos=3,
		label.margin=0.05,
		)

	kpPlotDensity(
		nuc,
		data.panel=1, 
		r0=0.55, 
		r1=1.0, 
		data=jobgr$contamination,
		col="#0278ff", 			#Medium blue
		) 

	#LOWER TOP PANEL
	#Panel background
	kpDataBackground(
		nuc,
		data.panel=2,
		r0=0.0, r1=0.45,
		col='#d8fdff',     #Pale blue
		) 

	kpAxis(nuc,
			data.panel=2,
			r0=0.0, 
			r1=0.45,
			cex=axislabelsize,
			numticks=axisticks,
			)

	kpPlotDensity(
		nuc,
		data.panel=2, 
		r0=0.0, 
		r1=0.45, 
		data=jobgr$nuclearnumts,
		col="#ffe102", 			#Yellow
		)

	kpAddLabels(
		nuc,
		labels='Found numts density', 
		data.panel=2, 
		r0=0.0, 
		r1=0.45, 
		cex=tracklabelsize, 
		srt=90, 
		pos=3,
		label.margin=0.05,
		) 

	#LOWER BOTTOM PANEL
	#Panel background
	kpDataBackground(nuc,
					data.panel=2,
					r0=0.55, 
					r1=1.00,
					col='#d8fdff',     #Pale blue
					) 

	kpAxis(nuc,
			data.panel=2,
			data=refgr$reference,
			r0=0.55, 
			r1=1.00,
			cex=axislabelsize,
			numticks=axisticks,
			)

	kpAddLabels(
		nuc, 
		data.panel=2, 
		labels='Ref numt density', 
		r0=0.55, 
		r1=1.00,
		cex=tracklabelsize, 
		srt=90, 
		pos=3,
		label.margin=0.05
		)

	kpPlotDensity(nuc,
		data.panel=2, 
		data=refgr$reference, 
		r0=0.55, 
		r1=1.00,
		col='#c402ff',				#Bright purple
		)

	dev.off()
}
