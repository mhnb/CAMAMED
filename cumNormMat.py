# Cumulative sum scaling factors.
# Calculates each column's quantile and calculates the sum up to and including that quantile.
#################################################################
# Cumulative sum scaling percentile selection
# Calculates the percentile for which to sum counts up to and scale by. Faster
# version than available in cumNormStat. Deviates from methods described in Nature Methods by
# making use of ro means for reference.
def cumNormStatFast(data):
	import os
	import numpy as np
	import pandas as pd
	import warnings
	le=0
	nr=np.size(data,0)
	nc=np.size(data,1)
	smat=np.zeros(nc,dtype=int)
	for i in range(0,nc):
		smat[i]=np.count_nonzero(data[:,i]>0)
		if smat[i]<=0:
			le=le+1;
	leng=int(max(smat))
	if le>0:
		print("Warning sample with one or zero features")
	else:
		smat2=np.empty((leng,nc))
		smat2[:]=np.NaN
		for i in range(0,nc):
			ix= np.nonzero(data[:,i]>0)
			da=data[ix,i]
			da=np.sort(da)
			#da.size
			smat2[leng-da.size:leng,i] = da
	p=np.linspace(0,1,np.size(smat2,0))
	rmat2=np.zeros((leng,nc),dtype=float)
	df = pd.DataFrame(smat2)
	rmat2=df.quantile(p)
	rmat2=rmat2.values
	smat2[np.isnan(smat2)]=0
	ref1=np.mean(smat2, axis = 1)
	nrr=np.size(rmat2,0)
	ncc=np.size(rmat2,1)
	diffr=np.zeros((nrr,ncc),dtype=float)
	for i in range(0,nc):
		diffr[:,i]=np.subtract(ref1,rmat2[:,i])
	diffr1=np.median(abs(diffr), axis = 1)
	ix= np.argmax(np.divide(abs(np.diff(diffr1)),diffr1[1:None]) > 0.1)
	x=float(ix+1)/float(len(diffr1))
	#print(x)
	if x<=0.5:
        	x=0.5
	return(x)
###############################################################
def cumNormMat(data,zi_p,g):
	import os
	import numpy as np
	import pandas as pd
	import warnings
	#os.system('clear')
	warnings.filterwarnings("ignore")
	cdp=os.getcwd()
	ixnz1= np.nonzero(np.sum(data, axis = 1)>=zi_p)
	ixnz= np.sum(data, axis = 1)>=zi_p
	ndata=data[ixnz,:]
	if g==1:
		file_p=cdp+'/files/gene_length.txt'
		gl=np.loadtxt(file_p, dtype=int)
		ngl=gl[ixnz]
		rr=np.size(ndata,0)
		for i in range(0,rr):
			a=np.divide(ndata[i,:],ngl[i])
			#print(np.sum(a))
			a=np.divide(a,0.00001)
			#print(np.sum(a))
			a=np.round(a)
			#print(np.sum(a))
			#print()
			ndata[i,:]=np.round(a)
	p=cumNormStatFast(ndata)		
	xx=ndata
	xx=xx.astype('float')
	xx[xx==0]='nan'
	df = pd.DataFrame(xx)
	qs=df.quantile(p)
	qs=qs.values
	
	newMat=np.zeros(qs.size,dtype=float)  
	xxx=np.subtract(ndata,np.finfo(float).eps)
	nr=np.size(ndata,0)
	nc=np.size(ndata,1)
	for i in range(0,nc):
		ix= np.nonzero(xxx[:,i]<=qs[i])
		newMat[i]=np.sum(xxx[ix,i])
	sl=1000.0
	newMat=np.divide(newMat,sl)
	nmat=np.zeros(ndata.shape,dtype=float)
	for i in range(0,nc):
		nmat[:,i]=np.divide(ndata[:,i],newMat[i])
	nr=np.size(data,0)
	nc=np.size(data,1)
	data=np.zeros((nr,nc),dtype=float)
	data[ixnz,:]=nmat
	return data
	