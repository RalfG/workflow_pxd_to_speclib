include: "get_data.snakefile"
include: "search_data.snakefile"
include: "make_speclib.snakefile"


rule targets:
	input:
		speclib/spectral_library.peprec,
		speclib/spectral_library.mgf
