
bcf.open<-function(filename,requireIndex)
	{
	assert(is.character(filename))
	assert(is.logical(requireIndex))
	.Call("RBcfFileOpen",filename,requireIndex);
	}
bcf.close<-function(fp)
	{
	.Call("RBcfFileClose",fp);
	}
bcf.hdr<-function(fp)
	{
	.Call("RBcfHeader",fp);
	}
	
bcf.hdr.nsamples<-function(fp)
	{
	.Call("RBcfNSamples",fp);
	}

bcf.hdr.samples<-function(fp)
	{
	.Call("RBcfSamples",fp);
	}
