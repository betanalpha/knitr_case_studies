L <- 50

ls <- rep(0:L, each=2)
buffered_ls <- sapply(1:length(ls), function(l) if(l %% 2 == 0) idx[l] + 0.5 else idx[l] - 0.5)

rho <- 0.25

rep_rhos <- rho**ls
write.table(data.frame(buffered_ls, rep_rhos), file = "low_rho.dat",
            row.names=FALSE, col.names=FALSE, sep=" ")

rho <- 0.95
rep_rhos <- rho**ls
write.table(data.frame(buffered_ls, rep_rhos), file = "high_rho.dat",
            row.names=FALSE, col.names=FALSE, sep=" ")

