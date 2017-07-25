def symmetry_stat():
	data = open(outfile_map_wholegen_dil_masked.bed)
	EU = data[:,2]#EU3 = model_3$frequency_EU
	AS = data[:,3]#AS3 = model_3$frequency_AS
	print (EU)
	print (AS)
	
	#return EU, AS

Symm = symmetry_stat