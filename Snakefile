include: "get_data.smk"
include: "search_data.smk"
include: "make_speclib.smk"


rule targets:
	input:
		"speclib/spectral_library.peprec",
		"speclib/spectral_library.mgf"
