diversity.sim<-function(local=200,form=20,form.a=10,ptaph=0.2,
	pform=0.2,ploc=0.2,pmist=0.1,min.spec=1000,max.spec=2000,nbins=100,
	disp=0.25,ext=0.25)

{
	form.b<-form-form.a
	#######################
	## Simulate clade
	#######################

	sim<-simFossilTaxa(0.25,0.25,mintaxa=min.spec,maxtaxa=max.spec,mintime=nbins,maxtime=nbins,plot=F,anag.rate=0.125,prop.bifurc=0.5)





	################################################################
	## Chose clade of interest(between 0.25 and 0.75 of the total)
	################################################################


	tree<-taxa2cladogram(sim)
	
	big.enough<-F
	small.enough<-F

	upper.limit<-0.75*nrow(sim)
	lower.limit<-0.25*nrow(sim)

	while(big.enough==F | small.enough==F)
	{
		node<-sample(1:tree$Nnode,1)+length(tree$tip.label)
		taxon.list<-FindDescendants(node,tree)
		if(length(taxon.list)<upper.limit)
		{
			small.enough<-T
		}
		else(small.enough<-F)
		if(length(taxon.list)>lower.limit)
		{
			big.enough<-T
		}
		else(big.enough<-F)
	}

	taxa.of.interest<-tree$tip.label[taxon.list]





	###############################################################
	## matrix of time bins; two colums (first and last of each bin) 
	###############################################################		

	time<-matrix(nrow=nbins,ncol=2)
	time[,1]<-seq(nbins,1)	
	time[,2]<-seq(nbins-1,0)






	##################################################################
	## create two objects: names and listnames	
	## listnames is vector of all the taxon names
	## names is a list nbins long. each element is listnames repeated
	##################################################################

	sequence<-seq(1:nrow (time))
	names<-as.list(sequence)
	for (i in 1:length(names))
	{
		names[[i]]<-array(dim=nrow(sim))
	}
	for (i in 1:length(names))
	{
		names[[i]]<-rownames(sim)
	}
	listnames<-rownames(sim)




	######################################################################################
	## totdiv: actual diversity of whole clade, including taxa outside clade of interest
	######################################################################################
	totdiv<-array(dim=nrow(time),data=0)
	for (i in 1:length(totdiv))
	{
		for (j in 1:nrow(sim))
		{
			if(sim[j,3]>time[i,2] && sim[j,4]<time[i,1])
			{
				totdiv[i]<-totdiv[i]+1
			}
		}
	}




	##############################################################
	## adiv: actual diversity in each bin of the clade of interest
	##############################################################
	adiv<-array(dim=nrow(time),data=0)
	for (i in 1:length(adiv))
	{
		for (j in 1:nrow(sim))
		{
			if(sim[j,3]>time[i,2] && sim[j,4]<time[i,1])
			{
				if(rownames(sim)[j]%in%taxa.of.interest)
				{	
					adiv[i]<-adiv[i]+1
				}
			}
			else(names[[i]][j]<-0)
		}
	}



	##############################################################
	## Find the actual range of the clade of interest
	##############################################################

	first.app<-min(which(adiv!=0))

	if(adiv[nbins]==0)
	{
		last.app<-max(which(adiv!=0))
	}

	if(adiv[nbins]!=0)
	{
		last.app<-nbins
	}


	clade.range<-first.app:last.app


	adiv<-adiv[clade.range]

	





	#######################################
	## Make it so each list element of 
	## "names" contains only the names of 
	## taxa present in that bin
	#######################################

	for (i in 1:length(names))
	{
		subnames<-which(names[[i]]!=0)
		names[[i]]<-names[[i]][subnames]
	}





	#######################################
	## create formations matrix
	#######################################


	sequence<-seq(1:nrow (time))
	form.pres<-as.list(sequence)
	for (i in 1:length(form.pres))
	{
		form.pres[[i]]<-matrix(nrow=length(names[[i]]),ncol=form,data=0)
		rownames(form.pres[[i]])<-names[[i]]
	}




	######################################################################
	## find first time bin for all taxa, randomly assign start formation
	######################################################################
	first.app<-vector(length=length(listnames))
	names(first.app)<-listnames
	for(i in 1:length(first.app))
	{
		first.app[i]<-max(which(time[,1]>=sim[i,3]))
		rand.form<-sample(1:form.a,1)
		form.pres[[first.app[i]]][names(first.app[i]),rand.form]<-1
	}





	#################################
	## dispersal and local extinction
	#################################
	#################################################################
	#################################################################
	## After running code below you have form.pres: a list nbins long
	## Each list element represents each time bin
	## In each list element is a matrix
	## Each column of that matrix is a formation
	## Each row is a taxon alive in that time bin
	## A score of 1 means the taxon is present in that formation
	#################################################################
	#################################################################

	for(i in 1:length(form.pres))
	{
		for(j in 1:nrow(form.pres[[i]]))
		{
			if(i!=first.app[rownames(form.pres[[i]])[j]])
			{
				form.pres[[i]][j,]<-form.pres[[i-1]][rownames(form.pres[[i]])[j],]
			}
			if(sum(form.pres[[i]][j,])!=form.a)
			{
				rand<-runif(1)
				if(rand<disp)
				{
					rand.form<-which(form.pres[[i]][j,]==1)
					while(rand.form==which(form.pres[[i]][j,]==1) || rand.form%in% which(form.pres[[i]][j,]==1))
					{
						rand.form<-sample(1:form.a,1)
					}
					form.pres[[i]][j,rand.form]<-1
				}
			}
			if(sum(form.pres[[i]][j,])!=1)
			{
				rand<-runif(1)
				if(rand<ext)
				{
					rand.form<-sample(which(form.pres[[i]][j,]==1),1)
					form.pres[[i]][j,rand.form]<-0
				}
			}
		}
	}






	#################################
	## create matrix of localities
	#################################
	########################################################################
	########################################################################
	## After running code below you have the object called matrices
	## This works the same as form.pres, but now each column in each matrix
	## is a locality
	########################################################################
	########################################################################
	local.per.form<-local/form
	sequence<-seq(1:length(totdiv))
	matrices<-as.list(sequence)
	for (i in 1:length(matrices))
	{
		matrices[[i]]<-matrix(nrow=totdiv[i],ncol=local,data=0)
		rownames(matrices[[i]])<-names[[i]]
	}
	for(i in 1:length(matrices))
	{
		for(j in 1:nrow(form.pres[[i]]))
		{
			there<-which(form.pres[[i]][j,]==1)
			for(k in 1:length(there))
			{
				matrices[[i]][j,(local.per.form*there[k]-(local.per.form-1)):(local.per.form*there[k])]<-1
			}
		}
	}



	###############################################################
	## vectors of sampling probs per time bin 
	###############################################################
	ptaphs<-rep(ptaph[1],nrow(time))
	if(length(ptaph)>1)
	{
		ptaphs<-sample(seq(ptaph[1],ptaph[2],0.1),nrow(time),replace=T)
	}

	pforms<-rep(pform[1],nrow(time))
	if(length(pform)>1)
	{
		pforms<-sample(seq(pform[1],pform[2],0.1),nrow(time),replace=T)
	}
	
	plocs<-rep(ploc[1],nrow(time))
	if(length(ploc)>1)
	{
		plocs<-sample(seq(ploc[1],ploc[2],0.1),nrow(time),replace=T)
	}




	##############################################################################
	## Taphonomic pruning
	## Remove some taxa from some of the localities
	## A specific taxon survives in a specific locality with the probability ptaph
	##############################################################################

	for(i in 1:length(matrices))
	{
		for(j in 1:nrow(matrices[[i]]))
		{
			there<-which(matrices[[i]][j,]==1)
			for(k in 1:length(there))
			{
				rand<-runif(1)
				if(rand > ptaphs[i])
				{
					matrices[[i]][j,there[k]]<-0
				}
			}
		}
	}	



	#####################################################################################
	## Prune unpreserved formations (in practice turn all scores of all taxa in 
	## all localities in that formation to 0)
	## Formation survives into the fossil record with a score of pform
	## Also produces a vector called forms.sampled; a 1 means that formation is sampled
	####################################################################################

	forms.sampled<-matrix(nrow=nbins,ncol=form,data=1)
	locs.sampled<-matrix(nrow=nbins,ncol=local,data=1)
	
	for (i in 1:length(matrices))
	{
		for(j in 1:form)
		{
			ran<-runif(1)
			if (ran > pforms[i])
			{
				matrices[[i]][,(10*j-9):(10*j)]<-0
				forms.sampled[i,j]<-0
				locs.sampled[i,(10*j-9):(10*j)]<-0
			}
		}
	}


	###################################################################################
	## Prune unsampled localities 
	## Locality survives into the fossil record with a score of ploc
	## Also produces matrix (locs.sampled) showing which localities have been sampled
	###################################################################################

	for (i in 1:length(matrices))
	{
		for (j in 1:ncol(matrices[[i]]))
		{
			ran<-runif(1)
			if (ran > plocs[i])
			{
				matrices[[i]][,j]<-0
				locs.sampled[i,j]<-0
			}
		}
	}




	#################################################################
	## Create list of occurences
	## Code below creates objects occur and occurences, which are 
	## lists nbins long,
	## Each list element is a vector of the taxa sampled in that bin 
	##(occur) or list of occurences (occurences)
	#################################################################

	sequence<-seq(1:length(adiv))
	occur<-as.list(sequence)	
	occurences<-occur	


	for (i in 1:length(matrices))
	{
		for (j in 1:nrow(matrices[[i]]))
		{
			num.occs<-sum(matrices[[i]][j,])
			occs<-rep(names[[i]][j],num.occs)
			if (j==1)
			{
				occurences[[i]]<-occs
			}
			else (occurences[[i]]<-c(occurences[[i]],occs))
			occur[[i]]<-unique(occurences[[i]])
		}
	}



	###################################################################
	## Calculate sample in bin diversity (sibdiv) of clade of interest
	###################################################################

	sibdiv<-array(dim=length(matrices),data=0)
	for (i in 1:length(occur))
	{
		sibdiv[i]<-length(which(occur[[i]]%in%taxa.of.interest))
	}



	#############################################################################
	## Remove time bins before and after clade of interest appears in the record 
	#############################################################################


	sibdiv<-sibdiv[clade.range]


	#############################################
	## Calculate sampled age range for each taxon 
	#############################################
	obs.ranges<-matrix(nrow=length(listnames),ncol=2,data=NA)
	rownames(obs.ranges)<-listnames
	for(i in 1:length(listnames))
	{
		for(j in 1:length(occur))
		{
			if(listnames[i]%in%occur[[j]])
			{
				obs.ranges[i,2]<-time[j,2]
			}
		}
		for(j in length(occur):1)
		{
			if(listnames[i]%in%occur[[j]])
			{
				obs.ranges[i,1]<-time[j,1]
			}	
		}
	}
	obs.ranges<-obs.ranges[which(is.na(rowSums(obs.ranges))==F),]
	


	#############################################
	## Calculate range through diversity (rtdiv)
	#############################################

	rtdiv<-array(dim=length(matrices),data=0)
	for (i in 1:length(rtdiv))
	{
		for (j in 1:nrow(obs.ranges))
		{
			if(obs.ranges[j,1]>time[i,2] && obs.ranges[j,2]<time[i,1])
			{
				if(rownames(obs.ranges)[j]%in%taxa.of.interest)
				{
					rtdiv[i]<-rtdiv[i]+1
				}
			}
		}
	}

	rtdiv<-rtdiv[clade.range]


	###########################################################
	## Create phylogeny of sampled taxa from clade of interest
	###########################################################


	prunedtree<-drop.tip(tree,which(tree$tip.label%in%rownames(obs.ranges)==F))
	prunedtree<-drop.tip(prunedtree,which(prunedtree$tip.label%in%taxa.of.interest==F))



	########################################################################
	## Add errors to phylogeny
	## Random nearest node interchange with probability pmist for each node
	########################################################################

	for (i in 1:prunedtree$Nnode)
	{
		x<-runif(1)
		if (x < pmist)
		{
			prunedtree<-rNNI(prunedtree)
		}
	}




	######################################################
	## Time Calibrate tree
	######################################################

	timetree<-timePaleoPhy(prunedtree,obs.ranges[prunedtree$tip.label,],type="basic",add.term=T)



	########################################################
	## Calculate phylogenetic diversity estimate (pdediv)
	########################################################

	pdediv<-array(dim=nrow(time),data=0)
	node.dates<-dateNodes(timetree)
	branch.range<-matrix(nrow=nrow(timetree$edge),ncol=2)
	for(i in 1:nrow(timetree$edge))
	{
		branch.range[i,1]<-node.dates[timetree$edge[i,1]]
		branch.range[i,2]<-node.dates[timetree$edge[i,2]]
	}
	for (i in 1:length(pdediv))
	{
		for (j in 1:nrow(branch.range))
		{
			if(branch.range[j,1]>time[i,2] && branch.range[j,2]<time[i,1])
			{
				pdediv[i]<-pdediv[i]+1
			}
		}
	}
	pdediv<-pdediv[clade.range]



	########################################################
	## Calculate squares
	########################################################


	squares<-array(dim=nrow(time),data=0)

	for(i in 1:length(occurences))
	{
		tax<-table(occurences[[i]][which(occurences[[i]]%in%taxa.of.interest)])
		S<-length(tax)
		s1<-length(which(tax==1))
		squares[i]<-S + s1^2 * sum(tax^2) / (sum(tax)^2 - s1 * S)
	}
	squares<-squares[clade.range]
	################################################################
	## Calculate: no.form.a (formations bearing clade of interest)
	##            no.form.b (formations bearing larger clade)
	##            no.loc.a (localities bearing clade of interest)
	##            no.loc.b (localities bearing larger clade)
	################################################################

	no.form.a<-array(dim=nbins)
	no.form.b<-array(dim=nbins)
	no.loc.a<-array(dim=nbins)
	no.loc.b<-array(dim=nbins)

	for(i in 1:length(occur))
	{
		fossil.bearing.locs<-vector(length=0)
		clade.bearing.locs<-vector(length=0)
		if(length(occur[[i]]!=0))
		{
			for(j in 1:length(occur[[i]]))
			{		
				fossil.bearing.locs<-unique(c(fossil.bearing.locs,which(matrices[[i]][occur[[i]][j],]==1)))
				if(occur[[i]][j]%in%taxa.of.interest)
				{
					clade.bearing.locs<-unique(c(clade.bearing.locs,which(matrices[[i]][occur[[i]][j],]==1)))
				}	
			}
		}		
		no.loc.a[i]<-length(clade.bearing.locs)
		no.loc.b[i]<-length(fossil.bearing.locs)
		no.form.a[i]<-length(unique(trunc(clade.bearing.locs/10)))
		no.form.b[i]<-length(unique(trunc(fossil.bearing.locs/10)))
	}	


	no.form.a<-no.form.a[clade.range]
	no.form.b<-no.form.b[clade.range]
	no.loc.a<-no.loc.a[clade.range]
	no.loc.b<-no.loc.b[clade.range]


	#############################################################################
	## Calculate: no.form.c (all formations that the clade "could" have been in)
	##            no.form.d (all formations)
	##            no.loc.c (all formations that the clade "could" have been in)
	##            no.loc.d (localities bearing larger clade)
	#############################################################################

	no.form.c<-rowSums(forms.sampled[,1:10])
	no.form.d<-rowSums(forms.sampled)
	no.loc.c<-rowSums(locs.sampled[,1:100])
	no.loc.d<-rowSums(locs.sampled)

	no.form.c<-no.form.c[clade.range]
	no.form.d<-no.form.d[clade.range]
	no.loc.c<-no.loc.c[clade.range]	
	no.loc.d<-no.loc.d[clade.range]



	##############################
	##############################
	#####
	##### Residuals: 8 sorts
	#####
	##############################
	##############################


	############################################################
	### Formations bearing clade of interest
	############################################################

	res.form.a<-resid(lm(log10(sibdiv+1)~log10(no.form.a+1)))
	


	############################################################
	### Formations bearing larger clade
	############################################################

	res.form.b<-resid(lm(log10(sibdiv+1)~log10(no.form.b+1)))



	############################################################
	### Formations where clade could have been preserved
	############################################################

	res.form.c<-resid(lm(log10(sibdiv+1)~log10(no.form.c+1)))



	############################################################
	### All formations
	############################################################

	res.form.d<-resid(lm(log10(sibdiv+1)~log10(no.form.d+1)))






	############################################################
	### Localities bearing clade of interest
	############################################################

	res.loc.a<-resid(lm(log10(sibdiv+1)~log10(no.loc.a+1)))



	############################################################
	### Localities bearing larger clade
	############################################################

	res.loc.b<-resid(lm(log10(sibdiv+1)~log10(no.loc.b+1)))



	############################################################
	### Localities where clade could have been preserved
	############################################################

	res.loc.c<-resid(lm(log10(sibdiv+1)~log10(no.loc.c+1)))



	############################################################
	### All localities
	############################################################

	res.loc.d<-resid(lm(log10(sibdiv+1)~log10(no.loc.d+1)))



	##############################################################
	## Create matrix of results
	##############################################################



	results.matrix<-array(data=NA,dim=c(24,length(clade.range)))
	rownames(results.matrix)<-c("True Diversity","SIB TDE","RT TDE","PDE","Squares","RDE Form A","RDE Form B","RDE Form C","RDE Form D","RDE Loc A","RDE Loc B","RDE Loc C","RDE Loc D","Proxy Form A","Proxy Form B","Proxy Form C","Proxy Form D","Proxy Loc A","Proxy Loc B","Proxy Loc C","Proxy Loc D","pform","ploc","ptaph")




	#############################################
	#############################################
	#####
	##### Put results into a matrix
	#####
	#############################################
	#############################################
	results.matrix[1,]<-adiv
	results.matrix[2,]<-sibdiv
	results.matrix[3,]<-rtdiv
	results.matrix[4,]<-pdediv
	results.matrix[5,]<-squares
	results.matrix[6,]<-res.form.a
	results.matrix[7,]<-res.form.b
	results.matrix[8,]<-res.form.c
	results.matrix[9,]<-res.form.d
	results.matrix[10,]<-res.loc.a
	results.matrix[11,]<-res.loc.b
	results.matrix[12,]<-res.loc.c
	results.matrix[13,]<-res.loc.d
	results.matrix[14,]<-no.form.a
	results.matrix[15,]<-no.form.b
	results.matrix[16,]<-no.form.c
	results.matrix[17,]<-no.form.d
	results.matrix[18,]<-no.loc.a
	results.matrix[19,]<-no.loc.b
	results.matrix[20,]<-no.loc.c
	results.matrix[21,]<-no.loc.d
	results.matrix[22,]<-pforms[clade.range]
	results.matrix[23,]<-plocs[clade.range]
	results.matrix[24,]<-ptaphs[clade.range]

	
	return(results.matrix)
}





