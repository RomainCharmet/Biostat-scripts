#chromosome number and p-values should be displayed in columns, SNPs must be sorted in chromosome:position order
c = which(is.na(a$P.value)==TRUE)
a = a[-c,]
#Remove NAs
chr = a$chr
#alternative color scheme (WARNING: random) rainbow = rgb(runif(22),runif(22),runif(22))
index = vector(length=length(unique(chr)))
v = 0
col = vector(length=length(chr))
for (i in 1:length(col))
{
        f = chr[i];
        if(f%%2!=0)
        {
	  col[i] = "orange"
	}
	else
	{
	  col[i] = "blue"
	}
        #col[i] = rainbow[f];
}
#Set dual color scheme, blue for evens and orange for odds, use colorblind friendly colors
for (i in 1:22)
{
        c = which(chr==i);
        index[i] = v+(length(c)/2);                                                                                                                                                                                       
        v = v+length(c);                                                                                                                                                                                                  
}
#Set chromosome labels positions on x-axis
logP = -log10(a$P.value)
#Create a manhattan plot
png("manhattan_plot.png",width=1000, height=1000)
par(mar=c(8,8,4,5)+0.1,mgp=c(5,2,0))
plot(logP,col=col,xaxt="n",yaxt='n',pch=16,cex=1.5,cex.lab=2.5,main="",xlab="Chromosomes",ylab="-log(P-value)")
abline(h=8.5,col="red",lwd=2)
abline(h=6,col="blue",lwd=2)
axis(1,at=index,labels=c(1:22),cex.axis=2.5)
axis(2,at=c(1:11),labels=c(1:11),cex.axis=2.5,las=1)
text()
text
dev.off()


shade = function(x1, y1, x2, y2, color = "light blue")
{
        n1 = length(x2)
        polygon(c(x1, x2[n1:1]), c(y1, y2[n1:1]), border = NA, col = color)
}
pvalue = a$P.value
nmax = length(pvalue)
stat = -log10(pvalue)
observ = sort(stat)
rang =(nmax:1)/(nmax+1)
expect = -log10(rang)

highval = sort(-log10(qbeta(0.975,1:nmax,nmax:1)))
lowval = sort(-log10(qbeta(0.025,1:nmax,nmax:1)))
#Compute inflation factor
lambda = median((a$Effect^2)/(a$StdErr^2))/0.456
#Create QQ_Plot
png("QQplot_MetaAnalysis_3studies.png", width=1000, height=1000)
par(mar=c(8,8,4,5)+0.1,mgp=c(5,2,0))
plot(x=expect,y=observ,ylim=c(0,9),xlab=expression(paste(-log[10], "expected Pvalue")),ylab=expression(paste(-log[10], "observed Pvalue")),pch=".", main="",cex.lab=2.5,cex.axis=2.5,las=1)
shade(expect,lowval,expect,highval)
points(expect,observ,pch=3,cex=0.5)
#text(x=1,y=6,labels=paste("lambda=",signif(lambda,digits=5),"",sep=""),cex=0.75)
dev.off()
