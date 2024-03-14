.libPaths("/nobackup/tmjj24/apps/R/x86_64-pc-linux-gnu-library/4.2/")

# Get arguements
args <- commandArgs(trailingOnly = TRUE)

sp_0 <- args[1]
sp_1 <- args[2]
sp_2 <- args[3]
output.dir <- args[4]

total.sp <- sum(as.numeric(c(sp_0,sp_1, sp_2)))

traits <- paste0("titia_", 1:total.sp)

write.table(cbind(traits, c(rep(0, sp_0),rep(1, sp_1),rep(2, sp_2))),
            paste0(output.dir, "/all_traits.txt"),
            quote = F, row.names = F,
            col.names = c("traits","species"), sep = '\t')
