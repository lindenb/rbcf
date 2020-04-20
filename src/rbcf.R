PRIVATE_BWA_LIBRARIES_LOADED=FALSE;

private_load_rbcf_libraries<-function()
	{
	if(!PRIVATE_BWA_LIBRARIES_LOADED)
		{
		##load.dynamic.libraries(c("libbwa.so","librbwa.so"))
		dyn.load("htslib/libhts.so");
		dyn.load("librbcf.so");		
		}
	PRIVATE_BWA_LIBRARIES_LOADED=TRUE;
	}

bcf.open<-function(filename)
	{
	private_load_rbcf_libraries();
	.Call("RBcfFileOpen",filename);
	}
bcf.close<-function(fp)
	{
	private_load_rbcf_libraries();
	.Call("RBcfFileClose",fp);
	}

