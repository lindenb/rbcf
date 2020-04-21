
bcf.open<-function(filename)
	{
	.Call("RBcfFileOpen",filename);
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
