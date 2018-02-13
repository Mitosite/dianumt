# set arguments
args=commandArgs(T)
input=args[1]
output=args[2]

if (length(args) < 2) { # help printout
  
  cat("required arguments:\n")
  cat("\tinput file\n")
  cat("\toutput file name (without extension)\n")
  
} else { # rest of the program
  
  if (!file.exists(input)) { # checking that input file exists
    stop("input file does not exist")
  }
  
  # import data and assign custom column names
  mydata <- read.table(input) # read data
  colnames(mydata) <- c("chr", "position", "depth")
  
  # sort the data according to chromosomes
  goodChrOrder <- paste(c(1:22,"X","Y","MT"),sep="")
  mydata$chr <- factor(mydata$chr, levels=goodChrOrder)
  chrs=sort(unique(mydata$chr)) # list chromosomes present in the data
  
  # output pdf file
  pdf(file=paste(output,".pdf",sep=""), width=9, height=length(chrs)*4)
  
  # graphical parameters
  par(mfrow=c(length(chrs),1)) # number of rows depending on chromosomes number
  par(mar=c(3,3,3,3)) # graph margins; c(bottom, left, top, right)
  
  # plotting instructions
  for (x in chrs) { # one plot per chromosome
    plot(mydata[mydata$chr==x,2], mydata[mydata$chr==x,3], type="l", ylim=c(0,150),
         main=paste("chromosome",x))
    axis(2,at=c(0,50,100,150)) # custom axis ticks
  }
  
  dev.off()
  
}
