plot_netview<- function(infile, popfile) {
	# infile=> plink dist; popfile=> metaData.
	library('netview')
	tmp<-read.table(infile)
	XHY_dist<-as.matrix(tmp)
	XHY_meta<-read.table(popfile, header=T)
	XHYOptions <- netviewOptions(selectionTitle="XHY k-Selection", nodeID="ID", nodeGroup="Group", nodeColour="Colour", communityAlgorithms=c("Walktrap", "Infomap", "Fast-Greedy"))
	graphsD3 <- netview(XHY_dist, XHY_meta, options = XHYOptions, project= unlist(strsplit(infile, '\\.'))[1], networkD3= TRUE, mst=T, save=T)

}
