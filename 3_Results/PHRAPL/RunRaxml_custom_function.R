RunRaxml_v2 <- function(raxmlPath=paste(getwd(),"/",sep=""),raxmlVersion="raxmlHPC",
	inputPath=paste(getwd(),"/",sep=""),mutationModel,iterations,Threads,
	seed=sample(1:10000000,1),outputSeeds=FALSE,discard=FALSE){
	phylipFilesList <- list.files(inputPath,pattern="*.phylip$",full.names=FALSE)
	seed.vec <- array()
	for(numb in 1:length(phylipFilesList)){
		inputFile=phylipFilesList[numb]
		seed=sample(1:10000000,1)
		currentSeed <- seed
		seed.vec <- c(seed.vec,currentSeed)
		outputFile <- paste(inputFile,".tre",sep="")
		if(mutationModel=="file"){
			mutationModelFile <- (utils::read.table(paste(inputPath,"mutationModel.txt",sep="")))
			thisMutationModel <- as.character(mutationModelFile[[1]][numb])
		} else{
			thisMutationModel <- mutationModel
		}

 		systemCall2 <- system(paste(raxmlPath,raxmlVersion," -w ",inputPath," -s ",inputPath,inputFile," -n ",
 			outputFile," -m ",thisMutationModel," -f d -N ",iterations," -p ",
 			currentSeed, "-T ",Threads,sep=""))

		#Dicard inessential RAxML output
		if(discard==TRUE){
			if(file.exists(paste(inputPath,inputFile,".reduced",sep=""))){
				unlink(paste(inputPath,inputFile,".reduced",sep=""))
			}
			if(file.exists(paste(inputPath,"RAxML_info.",inputFile,".tre",sep=""))){
				unlink(paste(inputPath,"RAxML_info.",inputFile,".tre",sep=""))
			}
			for(run in 0:(iterations-1)){
				if(file.exists(paste(inputPath,"RAxML_log.",inputFile,".tre.RUN.",run,sep=""))){
					unlink(paste(inputPath,"RAxML_log.",inputFile,".tre.RUN.",run,sep=""))
				}
				if(file.exists(paste(inputPath,"RAxML_parsimonyTree.",inputFile,".tre.RUN.",run,sep=""))){
					unlink(paste(inputPath,"RAxML_parsimonyTree.",inputFile,".tre.RUN.",run,sep=""))
				}
				if(file.exists(paste(inputPath,"RAxML_result.",inputFile,".tre.RUN.",run,sep=""))){
					unlink(paste(inputPath,"RAxML_result.",inputFile,".tre.RUN.",run,sep=""))
				}
			}
		}
	}
	if(outputSeeds==TRUE){
		write.table(seed.vec[-1],file=paste(inputPath,"RAxML.seeds.txt",sep=""))
	}
	return(systemCall2)
}

