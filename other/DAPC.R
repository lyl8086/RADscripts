plot_dapc<- function(infile, popfile) {
	library('adegenet')
	
	x<<-read.genepop(infile, ncode=3)
	tb1<<-read.table(popfile, header=T)
	mycol<<-as.character(unique(tb1$Colour))
	mypop<<-as.character(unique(tb1$Group))
	
	dapc1<-dapc(x, pca.select="percVar", perc.pca=80, n.da=3) #80% variance retained.
	gp<-round(dapc1$n.pca/5) #Every 5 PCs for optim.a.score.
	best<<-optim.a.score(dapc1, n.sim=100, n=gp)
	ask<-readline("Continue to plot scatter? ")
	if (ask == 'y') {
		dapc2<<-dapc(x, n.pca=best$best, n.da=3)
		scatter(dapc2, col=mycol, cstar=0, solid=1, cex=1.5, label=mypop, scree.da=F, legend=T, txt.leg=mypop)
	}
	
	ask<-readline("Continue to plot component? ")
	if (ask == 'y') {
		compoplot(dapc2, lab="", col=mycol, txt.leg=mypop, lwd=1.5,)
	}
	
	ask<-readline("Continue to plot assignment? ")
	if (ask == 'y') {
		assignplot(dapc2)
	}
}
