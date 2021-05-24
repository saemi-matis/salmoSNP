library("ggplot2")
islandStofnfiskurImiss<-read.table("SsaTrack_Island_Stofnfiskur_1-8_bin_missing.imiss", h=T)

p_imiss <- ggplot(islandStofnfiskurImiss, aes(x=FID, y=F_MISS)) + 
  geom_boxplot(notch=TRUE)

p_imiss
