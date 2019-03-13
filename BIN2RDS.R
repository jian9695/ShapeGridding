rm(list = ls())
gc()

inputdir <- 'E:/Projects/Gridding/Bin/Samples/shapefile_grid_intersection_fractions/Yilong/'


files <- list.files(inputdir,pattern = 'bin', ignore.case = TRUE)
for(filename in files)
{    
  infile <-paste0(inputdir,filename)
  outfile <-  paste0(inputdir,substr(filename, 1,nchar(filename)-3),'rds')
  if(file.exists(infile) & !file.exists(outfile)){
    zz <- file(infile, "rb")
    fractionArray.num<-readBin(zz, integer(),size = 4,n = 1)# read the number of time profiles
    fractionArray.feaIds <- readBin(zz, integer(),size = 4,n = fractionArray.num)
    fractionArray.gridIds <- readBin(zz, integer(),size = 4,n = fractionArray.num)
    fractionArray.fractions <- readBin(zz, double(),size = 8,n = fractionArray.num)
    fractionArray <- data.frame(fractionArray.feaIds,fractionArray.gridIds,fractionArray.fractions)
    names(fractionArray) <- c('feaId','gridId','fraction')
    saveRDS(fractionArray, file = outfile)
    close(zz)
  }
}

