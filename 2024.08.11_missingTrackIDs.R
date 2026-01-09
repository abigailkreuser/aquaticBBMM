## whos missing?!
library(dplyr)

doHave <- list.files('data/processedData/results/20250204bbmmRasterDataFullTracks_res4km_bmvar464.1sqkm_locerr2.5k')
head(doHave)

doHave <- stringr::str_extract_all(doHave, "\\d+")  %>%
  as.vector() %>%
  as.numeric()
head(doHave)

mvmnt <- read.csv("data/20250203_mvmnt_with_dist.csv")
wanted <- unique(mvmnt$trackID)
length(wanted)
dontHave <- wanted [!(wanted %in% doHave)]


write.csv(dontHave, "data/DONTHAVES_fullday_02.06.csv")

# 
# ## for local use
# 
# 
 wanted <- read.csv("data/DONTHAVES_fullday_02.06.csv")
# 
 wanted <- unique(wanted$x)
head(wanted)
summary(wanted)
 length(wanted)
 missingTracks <- mvmnt[mvmnt$trackID %in% wanted,]

#
