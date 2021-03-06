# iqpirquant

Quantitation of cross-links and deadend peptides based on 2-plex or multiplexed iqPIR cross-linkers


Place the following files in a common directory on Linux for running iqpirQuant.php:

iqpirQuant.php
MyIqPirMplex.php
MyStatistics.php
MyPepXML.php
MyDeadendPepXML.php
t_test_pvals.txt
anova_pvals.txt


Running iqpirQuant.php on data cross-linked with multiplexed iqPIR requires:

XLinkProphet results in pepXML format


Several parameter files (see ParameterFileInformation.pdf for details):

	iqpir.params file in format generated by iqpirQuant.php WRITEPARAMS command with specified input file and output directory
	sample_map file specifying raw files and associated samples
	sample_channel file specifying sample cross-linker stump modification and reporter masses
	sample_ratios file specifying numerator and denominator pairwise channel ratios
	
where a sample is defined as a multiplexed sample comprising 2-6 individual samples, each cross-linked with a distinct iqPIR cross-linker

In directory with source files, type iqpirQuant.php to see list of all command options.  
	(See iqpirquantRunCommands.pdf for detailed description)
	

Output files are written for each sample specified in the sample_map and all pairwise channels specified in the sample_ratios files

Foreach sample N:

	sample.N.iqpir.txt: Calculated apportionments (relative concentrations) of quantified ions among the channels and their corresponding envelope intensity differences
	
	sample.N.numids.txt: Number of identifications of cross-links in each channel
	
	sample.N.apports.txt: Calculated cross-link apportionments (mean, stdev, num contributing ions, outliers)
	
	sample.N.respairs.txt: Cross-link data aggregated at the protein residue pair level 
	
	sample.N.reshubs.txt: Cross-link data sorted by cross-linked protein residues
	
	sample.N.protpairs.txt: Cross-link data aggregated at the cross-linked protein pair level
	

For each pairwise numerator (num) and denominator (denom) channels of sample N:

	sample.N.num-denom.chroms.txt: Calculated log2ratio of quantified ions among the channels and their signal to noise and corresponding envelope intensity differences
	sample.N.num-denom.ratios.txt: Cross-link log2ratios (mean, stdev, num contributing ions, outliers, p-values)
	sample.N.num-denom.respairs.txt: Cross-link data aggregated at the protein residue pair level
	sample.N.num-denom.reshubs.txt: Cross-link data sorted by cross-linked protein residues
	sample.N.num-denom.protpairs.txt: Cross-link data aggregated at the cross-linked protein pair level
