## Main function

                 
fCCAr <- function(peaks, bigwigs, labels, splines=10, nbins=100, ncan=5 , tf=c(), main="", bar=NULL ){    
		
	######### Parameter description ###################################
	#
	# peaks: BED file. Column 1: chr, Column 2: start, Column 3: end (REQUIRED)
	# bigwigs: a vector of characters with the path to bigwigs files (REQUIRED)
	# labels: IDs for each sample. Replicates should have the same label and BE number-ORDERED(vector of characters, REQUIRED)
	# splines: number of equidistant cubic B-splines to smooth the data and canonical variate weight functions (default, 15)
	# nbins: single ‘integer’ value denoting how many bins there should be for each window (genomation, BioC) (default, 100)
	# ncan: number of canonical components to report in the results. Cannot be higher than number of splines (default, 15)
	# tf: plot results involving only this TF or replicate (character) (default, empty vector)
	# main: title of the plot (fCCAr.pdf) (default, No title)
	# bar: In the barplot, plot only first 'bar' interaction after ranking by F-value (default, NULL)
	# 
	##################################################################

	require(genomation)
	require(fda)
	require(utils)
	require(RColorBrewer)
	require(gplots)
	require(ggplot2)
	source("multiplot.R")
	options(warn=-1)
	
	#Read genomic regions in BED format
	print("Reading peaks...")
	peaks <- readGeneric(file=peaks, chr=1, start=2, end=3) #[1:1000]  #readBed
	
    
    if (ncan > splines | ncan > length(peaks) ){ print("ncan should not be higher than the number of splines or peaks. Please lower the value of ncan.")   }

	if (ncan <= splines | ncan <= length(peaks)){
		
		print("Starting fCCAr...")

		SPLINES = splines   #
		L = nbins           #bp


		#Modify labels to detect replicates (assume replicates are ordered 1,2,...)
		d <- duplicated(labels)
		new_labels <- c()
		for (l in 1:length(labels)) {
	
			if (d[l] == FALSE) {  counter =1; 				new_labels[l]  <- paste(labels[l], 1, sep="_Rep") }
			if (d[l] == TRUE)  {  counter = counter +1;  	new_labels[l]  <- paste(labels[l], counter, sep="_Rep")  } 	
		
		}
	
	

		#Create B-spline basis (cubic, order 4) of the signals before FDA:
		bspl <- create.bspline.basis(rangeval=c(-L/2,L/2), nbasis=SPLINES, norder=4)   # Cubic B-splines
		argvalsBS <- (-L/2):(L/2)

		
		#Read bigWig Files 
		fdaData <- list()

		for (i in 1:length(bigwigs)) {
			print(paste(c("Reading bigWig file...",i,"/",length(bigwigs)),collapse="")  )

			fdamatrix <- matrix(0.0, ncol=L+1, nrow= length(peaks) ) 
			fdamatrix  <- ScoreMatrixBin(target = bigwigs[i], bin.num = L+1, windows = peaks, type="bigWig",rpm=F, strand.aware = TRUE, bin.op="max" )
			fdaData[[i]] <- Data2fd(y=t(fdamatrix), argvals= argvalsBS, basisobj=bspl) 		

			#plot(colMeans(fdamatrix), type='l', ylab="mean normalized signal", xlab="distance to TSS", main=sample_label[i])	
			#plot(fdaData[[i]], ylab="normalized signal", xlab="distance to TSS", main=sample_label[i])	

		}
		rm(fdamatrix)
		length(fdaData)



		#Prepare Data for fCCA

		x <- list()                 		# list to store output of cca.fd
		scc <- c()  						# squared Canon. Corr.
		sccM <- c() 						# MAX squared Canon. Corr.
		pair <- c() 						# labeling
		Spair <- c() 						# labeling for S
		co <- combn(x=c(new_labels), m=2)  	# all possible pairwise combinations


		# Select Sample of Interest (if any)
				
		if(length(tf)==1){
			tf_co<-c()
			for (j in 1:ncol(co)){
		
		 		if ( length( grep(patter=tf, x=co[,j]) ) >0  ) tf_co <- c(tf_co, j)
	
			}
		co <- co[,tf_co]
		}
	
	
	
		w = 1/ (1:ncan)
		Ma <-  sum(w) 	# maximum possible value for w
		S <- c()  		# weighted sums

		for (i in 1:ncol(co)) {
			print(paste(c("Performing fCCA in pair...",i,"/",ncol(co)) ,collapse="") )
		
			file1 <- which(new_labels==co[1,i])
			file2 <- which(new_labels==co[2,i])
		
			print(paste(c(new_labels[file1],"...vs...",new_labels[file2]) , collapse="") )
	
			x[[i]] = cca.fd(fdobj1=fdaData[[file1]], fdobj2=fdaData[[file2]], ncan = ncan, ccafdPar1=fdPar(bspl, 0, 0), ccafdPar2=fdPar(bspl,0, 0), centerfns=TRUE)
	
		}
    	#NCAN calculated in 'cca.fd' is always the number of splines
		#plot.cca.fd(x)


		for (i in 1:ncol(co)) {
	
			scc <- c(scc, x[[i]]$ccacorr[1:ncan]   )
			sccM <- c( sccM, rep( max (x[[i]]$ccacorr[1:ncan]), ncan) )

			file1 <- which(new_labels==co[1,i])
			file2 <- which(new_labels==co[2,i])
			pair <- c(pair,  rep( paste(new_labels[c(file1,file2)]  , collapse="_vs_") ,ncan)  )
			Spair <- 	c(Spair, paste(new_labels[c(file1,file2)], collapse="_vs_" ) )

	
			#calculate weigthed sum
			S[[i]]= sum(w* x[[i]]$ccacorr[1:ncan]    )

		}


		if ( is.null(bar) == TRUE ) {bar <- ncol(co) }  ##barplot(100* S/Ma, ylim= c(0,100)  )


		#Colormap
		#colfunc <- colorRampPalette(c("cyan","blue", "black","darkred","red" ,"green" )) 
  		#colfunc <- colorRampPalette(c("#132B43","#56B1F7"))
  		colfunc <- colorRampPalette(brewer.pal(10, "Spectral"))

  

		#ggplot2
		ggData <- data.frame(scc = scc, pair = pair, variables = rep(1:ncan, ncol(co)), sccM = sccM)	
		ggData <- transform(ggData,  pair = reorder(pair, sccM))
		p1 <- ggplot(ggData, aes(x=variables, y=scc, group=pair))  + geom_line(aes(colour = pair)) + ylim(0,1)  + theme_classic() +  theme(legend.position='none', plot.title = element_text(lineheight=.8, face="bold")) + xlab("Canonical variable (k)") + ylab( expression(Squared~canonical~correlation~(R[k] ^{2}) ) )  +   geom_point(size=0.9, shape=5, aes(colour=pair)) + scale_x_continuous(breaks=1:ncan) + scale_colour_manual(values = colfunc(ncol(co))) + ggtitle(main) # + annotate("text", size=3, x = 5, y = 0.69, label = "PHF8", colour=rev(colfunc(ncol(co)))[3] )  + annotate("text", size=3, x = 3, y = 0.55, label = "DPY30", colour=rev(colfunc(ncol(co)))[6] )      ###+ scale_y_log10()  + scale_x_log10() # + scale_y_log10() + scale_colour_hue(l=60)
		ggData <- ggData[which(ggData$variables==1), ]
		ggData$CL <- rev( colfunc(ncol(co))  )


		ggData2 <- data.frame(S = 100* (S/Ma), pair = Spair  )   
		ggData2 <- transform(ggData2,  pair = reorder(pair, S))

		ggData3 <- merge(x=ggData2, y=ggData, by='pair',  sort = FALSE )

		temp <- ggData3$CL[1:bar]
		ggData3 <- ggData3[order(-ggData3$S),]
		ggDataTXT <- ggData3

		ggData3 <- ggData3[1:bar,]
 		ggData3$CL <- temp
 		

		#print(head(ggData3))
		p2 <- qplot(pair, S, data=ggData3, geom="bar", stat="identity", fill=factor(pair),width=.6)  + coord_flip()  + ylim(0,100) + theme_classic()   +  theme(legend.position='none', text = element_text(size=10) ) + ylab("F(%)") + xlab("")  + scale_fill_manual(values = rev(ggData3$CL) ) #, fill=factor(Spair) + scale_fill_hue(l=60)


		print("Saving fCCAr.pdf...")
		pdf("fCCAr.pdf", height=6, width=3.5)	
		multiplot(p1, p2, cols=1)
		dev.off()
		
		print("Saving fCCAr.txt...")
		colnames(ggDataTXT) <- c("samples","F","squ_can_corr_k_1","k","squ_can_corr_k_1","color"	)
		write.table(x=ggDataTXT[,c(1,2,4,5,6)] , file = "fCCAr.txt", append = FALSE, quote = F, sep = "\t", row.names = F, col.names = T)
	
		print("Done...")
		
		rm(temp)
		#return(ggData3$S)
	
	}
}
