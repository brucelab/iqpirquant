#!/usr/bin/php
<?php


include_once("MyPepXML.php"); // use this library for the insert statement
include_once("MyIqPir.php"); // use this library for the insert statement
include_once("MyIqPirMplex.php"); // use this library for the insert statement


function describeParams($multiplex = true) {
	echo "The sample map file contains ".($multiplex ? "3" : "4")." required tab-delimited columns listed below\n";
	echo "rawfile\t\tname of file without suffix (either full path, or assumed to be present in masschroq_dir of masschroq.params file) [e.g. 103119_3672_TACd8_3679_shamd0_F11_14_mango_4h_2]\n";
	if(! $multiplex) echo "fwd/rev\t\tfwd for RH/SH, rev for SH/RH [e.g. fwd]\n";
	echo "sample_numer\tunique number for a computed ratio (with light and heavy conditions) [e.g. 1 for ratio 1]\n";
	echo "biorep_number\tbiological replicate contributing to sample_number ratio (each biorep, coming from a different mixed heavy and light labeled sample, is normalized on its own to ensure median log ratio of 0) [e.g. 1]\n";
	echo "\n";
	//echo "Note that the sample_number and biorep_number columns of the sample map file are only used during the LOGRATIOS analysis step,\nwhereas the name, rawfile, fwd/rev, and alignment group columns are used by MassChroQ during the RUN step\n";
	if(! $multiplex) {
		echo "Once the LOGRATIOS step is completed for a sample map, one must create a sample def file before upload to XLinkDB using a browser\n\n";
		echo "The sample def file contains 4 required tab-delimited columns listed below, including the header line\n";
		echo "sample\t\t\tratio number corresponding to sample_number value in sample map file to which the ratios correspond [e.g. 2]\n";
		echo "experiment_condition\tunique brief name of the experiment (light) condition (this will be displayed in XLinkDB tables) [e.g. TAC]\n";
		echo "reference_condition\tbrief name of the reference (heavy) condition [e.g. sham]\n";
		echo "description\t\tbrief description of the ratio [e.g. TAC vs sham 6 biological replicates and 2 technical replicates for each bio sample]\n";
		echo "\n";
	}
	else {
		echo "Once the LOGRATIOS step is completed for a sample map, pairwise channel ratios are calculated based on the sample_ratios file specified in your iqpir.params file\n\n";
		echo "The sample_ratios file contains 7 required tab-delimited columns listed below, including the header line\n";
		echo "sample\t\t\tratio number corresponding to sample_number value in sample map file to which the ratios correspond [e.g. 2]\n";
		echo "experiment_condition\tunique brief name of the experiment (numerator) condition (this will be displayed in XLinkDB tables) [e.g. TAC]\n";
		echo "experiment_channel\tname of isobaric cross-linker for numerator [e.g. 826]\n";
		echo "reference_condition\tbrief name of the reference (denominator) condition [e.g. sham]\n";
		echo "reference_channel\tname of isobaric cross-linker for denominator [e.g. 808]\n";
		echo "description\t\tbrief description of the ratio [e.g. TAC vs sham 6 biological replicates and 2 technical replicates for each bio sample, or just a space '' if you don't want a description]\n";
		echo "normalize\t\tyes or no: whether to subtract median sample cross-link log2ratio to center distribution at 0 [e.g. yes]\n";
		echo "\n";
	}
	echo "New ratios can be computed without re-running RUN by modifying the sample_number and biorep_number columns of the sample map\n";
	echo "and re-running LOGRATIOS, then creating a new sample def file referenceing the sample_numbers for upload to XLinkDB through a browser\n\n";
}

function writeParams($print = false, $multiplex = true) {
	$file = $print ? 'php://stdin' : "iqpir.params.new";
	//echo "Here ready to write file {$file}....\n"; exit(1);
	if(! $print && file_exists($file)) {
		echo "file {$file} already exists.  Please first delete it if you want to overwrite it with the sample copy.\n";
		exit(1);
	}
	$fout = fopen($file, "w");
	fwrite($fout, "sample_map\t[full-path-to-samplemap-file]\n");
	if($multiplex) {
		fwrite($fout, "sample_ratios\t[full-path-to-samplemap-file]\n");
	}
	fwrite($fout, "email\t[your-address-to-which-notifications-are-sent]\n");
	if(! $multiplex) fwrite($fout, "heavy_modmass\t[327]\n");
	fwrite($fout, "xlinkprophetfile\t[full-path-to-xlinkprophet-pepxml]\n");
	fwrite($fout, "fdr\t[e.g. 0.01]\n");
	fwrite($fout, "filter_crit\t[composite_probability]\n");
	fwrite($fout, "iqpir_output_dir\t[full-path-directory-in-which-to-write-results]\n");
	if(! $multiplex) fwrite($fout, "normalize\t[true or false]\n");

	if($multiplex) {
		fwrite(STDOUT, "\nOPTIONAL parameter to specify sample_channels file with non-default channel stump and reporter masses (see below)\n");
		fwrite(STDOUT, "sample_channels\t[full-path-to-samplechannels-file]\n\n");
		fwrite(STDOUT, "\tOnly include this parameter if specify in sample_channels file for each sample non-default channel stump modification masses in increasing mass order, and reporter masses in decreasing mass order with format:\n");

		fwrite(STDOUT, "\tsample\tmodmass1,modmass2,...\treportermass1,reportermass2,...\n");
		fwrite(STDOUT, "\tDEFAULT VALUE:\n1\t325.127386,327.134095,329.14080505,330.13784005,332.14454905,334.15125905\t825.490275,821.47686,817.463,815.469,811.455947,807.442528\n");
		fwrite(STDOUT, "\t2\t325.127386,327.134095,329.14080505,330.13784005,332.14454905,334.15125905\t825.490275,821.47686,817.463,815.469,811.455947,807.442528\n");
		fwrite(STDOUT, "\t......\n");
	}



	fclose($fout);
	if(! $print) echo "\nFile {$file} written.  Edit the file, entering your relevant information for the analysis\n";

	if(! $multiplex) {
		echo "Sample map format: 4 tab delimited columns for raw file name, fwd/rev, sample #, sample biorep #\n";
		
		return;
	}

	$file = $print ? 'php://stdin' : "sample_ratios.new";
	if(! $print && file_exists($file)) {
		echo "file {$file} already exists.  Please first delete it if you want to overwrite it with the sample copy.\n";
		exit(1);
	}
	$fout = fopen($file, "w");
	fwrite($fout, "sample\texperiment_condition\texperiment_channel\treference_condition\treference_channel\tdescription\tnormalize\n");
	fclose($fout);
	if(! $print) echo "\nFile {$file} written.  Edit the file, entering your relevant information for the anlaysis, and change iqpir.params sample_ratios to reference this file\n";

	echo "Sample map format: 3 tab delimited columns for raw file name, sample #, sample biorep #\n";


}

function runCombination($params) {
global $self;
$prefix = $self . " " . $params;
$commands = array("{$prefix} RUN", // 0
			"{$prefix} RUN DE", // 1
			"{$prefix} LOGRATIOS", // 2
			"{$prefix} LOGRATIOS DE", // 3
			"{$prefix} LOGRATIOS INTRA", // 4
			"{$prefix} XL_DE_NORMALIZE", // 5
			"{$prefix} LOGRATIOS USE_CURRENT_NORM", // 6
			"{$prefix} LOGRATIOS DE USE_CURRENT_NORM", // 7
			"{$prefix} LOGRATIOS INTRA USE_CURRENT_NORM", // 8
			"{$prefix} DE_INTRA_PROTEINRATIOS" // 9
				);
for($k = 5; $k < 10; $k++) {
	echo $commands[$k] . "\n";
	system($commands[$k]);
}

}

function createMultiplexTwoChannelParams($twochannelparams) {
$myfile = fopen($twochannelparams, "r");
$params = array();
while($line = rtrim(fgets($myfile))) {
	$next = explode("\t", $line);
	$params[$next[0]] = $next[1];
}
fclose($myfile);
$dir = getcwd() . "/";
$new_samplemap = $dir . "sample_map_multi.txt"; //$params['sample_map'] . "multi";
$fout = fopen($new_samplemap, "w");
$myfile = fopen($params['sample_map'], "r");
$sample_orient = array();
while($line = rtrim(fgets($myfile))) {
	$next = explode("\t", $line);
	fwrite($fout, $next[0] . "\t" . $next[2] . "\t" . $next[3] . "\n");
	if(! array_key_exists($next[2], $sample_orient)) $sample_orient[$next[2]] = $next[1];
	if($next[1] != $sample_orient[$next[2]]) {
		echo "Error: cannot make multiplex params when have fwd and rev orientations for a sample ({$next[2]})\n";
		exit(1);
	} 
}
fclose($myfile);
fclose($fout);
$params['sample_map'] = $new_samplemap;
unset($params['heavy_modmass']);
$normalize = $params['normalize'] = "true";
unset($params['normalize']);
if($params['iqpir_output_dir'][strlen($params['iqpir_output_dir'])-1] != "/") $params['iqpir_output_dir'] .= "/";
$params['iqpir_output_dir'] .= "multiplex/";
$params['sample_channels'] = $dir . "sample_channels.txt";
$fout = fopen($params['sample_channels'], "w");
foreach($sample_orient as $sample => $orient) fwrite($fout, $sample . "\t" . "325.127386,327.134095\t811.455947,807.442528\n");
fclose($fout);
$params['sample_ratios'] = $dir . "sample_ratios.txt";
$fout = fopen($params['sample_ratios'], "w");
fwrite($fout, "sample\texperiment_condition\texperiment_stumpmodmass\treference_condition\treference_stumpmodmass\tdescription\tnormalize\n");
foreach($sample_orient as $sample => $orient) {
	fwrite($fout, $sample . "\t");
	if($orient == "rev") fwrite($fout, "808\t807.442528\t3812\t811.455947\t\t");
	else fwrite($fout, "812\t811.455947\t808\t807.442528\t\t");
	fwrite($fout, $normalize ? "yes" : "no");
	fwrite($fout, "\n");
}
fclose($fout);
$fout = fopen("iqpir_multi.params", "w");
$header = array("sample_map", "sample_channels", "sample_ratios");
for($k = 0; $k < count($header); $k++) {
	fwrite($fout, $header[$k] . "\t" . $params[$header[$k]] . "\n");
	unset($params[$header[$k]]);
}
foreach($params as $param => $value) fwrite($fout, $param . "\t" . $value . "\n");
fclose($fout);
echo "Multichannel params file iqpir_multi.params written\n";
}

// one filepath for each sample
function createParamFiles($filepaths, $email, $multichannel = false) {
$files = array();
$dir = getcwd() . "/";
for($k = 0; $k < count($filepaths); $k++) {
	array_push($files, glob($filepaths[$k]));
	if(count($files[count($files)-1]) == 0) {
		echo "Error: no files found in {$filepaths[$k]}\n";
		exit(1);
	}
}
$outfile = $multichannel ? "sample_map_multi.txt" : "sample_map.txt";
$fout = fopen($outfile, "w");
for($k = 0; $k < count($files); $k++) {
	for($j = 0; $j < count($files[$k]); $j++) {
		$nextpos = strrpos($files[$k][$j], "/");
		if($nextpos !==false) $files[$k][$j] = substr($files[$k][$j], $nextpos + 1);
		fwrite($fout, $files[$k][$j] . ($multichannel ? "" : "\tfwd") . "\t" . ($k+1) . "\t1\n");
	}
}
fclose($fout);
if($multichannel) {
	$fout = fopen("sample_channels.txt", "w");
	for($k = 0; $k < count($files); $k++) {

		fwrite($fout, ($k+1) . "\t325.127386,327.134095,329.14080505,330.13784005,332.14454905,334.15125905\t825.490275,821.47686,817.463,815.469,811.455947,807.442528\n");
	}
	fclose($fout);
	$fout = fopen("sample_ratios.txt", "w");
		fwrite($fout, "sample\texperiment_condition\texperiment_stumpmodmass\treference_condition\treference_stumpmodmass\tdescription\tnormalize\n");
	for($k = 0; $k < count($files); $k++) {

		fwrite($fout, ($k+1) . "\t327\t327.134095325\t325.127386\t\tno\n");
		fwrite($fout, ($k+1) . "\t329\t329.14080505\t325.127386\t\tno\n");
		fwrite($fout, ($k+1) . "\t330\t330.13784005\t325.127386\t\tno\n");
		fwrite($fout, ($k+1) . "\t332\t332.14454905\t325.127386\t\tno\n");
		fwrite($fout, ($k+1) . "\t334\t334.15125905\t325.127386\t\tno\n");
	}
	fclose($fout);
	$fout = fopen("iqpir_multi.params", "w");
	fwrite($fout, "sample_map\t{$dir}sample_map_multi.txt\n");
	fwrite($fout, "sample_channels\t{$dir}sample_channels.txt\n");
	fwrite($fout, "sample_ratios\t{$dir}sample_ratios.txt\n");
	fwrite($fout, "email\t{$email}\n");
	fwrite($fout, "xlinkprophetfile\t{$dir}iprophet-xl.pep.xml\n");
	fwrite($fout, "deadendprophetfile\t{$dir}Deadend_search/iprophet.pep.xml\n");
	fwrite($fout, "fdr\t0.01\n");
	fwrite($fout, "iqpir_output_dir\t{$dir}multiplex/\n");
	fwrite($fout, "longarm_reporters\ttrue\n");
	fclose($fout);
	echo "Multi-channel params file for ".count($files)." samples writtent to iqpir_multi.params\n";
}
else {
	$fout = fopen("iqpir.params", "w");
	fwrite($fout, "sample_map\t{$dir}sample_map.txt\n");
	fwrite($fout, "email\t{$email}\n");
	fwrite($fout, "heavy_modmass\t327\n");
	fwrite($fout, "xlinkprophetfile\t{$dir}iprophet-xl.pep.xml\n");
	fwrite($fout, "deadendprophetfile\t{$dir}Deadend_search/iprophet.pep.xml\n");
	fwrite($fout, "fdr\t0.01\n");
	fwrite($fout, "iqpir_output_dir\t{$dir}quant/\n");
	fwrite($fout, "normalize\tfalse\n");
	fwrite($fout, "include_reporters\tfalse\n");
	fclose($fout);
	echo "Two-channel params file for ".count($files)." samples writtent to iqpir.params\n";
}
}

$user = get_current_user();
$multichannel = true; //$user == "keller" || $user == "root" || $user == "";
$self = "iqpirQuant.php";
if(count($argv) < 2) {
	echo "\n  Usage: \033[1m{$self} < iqpir.params file > TASK\033[0m\n";
	echo "  TASKS: RUN (DE)                                    [computes quantitation for peptide and stump-containig fragments in each scan (DE for deadend)]\n";
	echo "         LOGRATIOS (DE)                              [after RUN: computes crosslink log ratios for all samples (DE for deadend)]\n";
	echo "         LOGRATIOS SAMPLES=x,y,... (DE)              [after RUN: computes crosslink log ratios for one or more specified samples x,y,... (DE for deadend)]\n";
	echo "         LOGRATIOS INTRA                             [after RUN: computes intra-protein crosslink log ratios for all samples]\n";
	echo "         LOGRATIOS INTRA SAMPLES=x,y,...             [after RUN: computes intra-protein crosslink log ratios for specified samples]\n";
	if($multichannel) echo "         LOGRATIOS APPORTS(=x,y,...)                 [after RUN: computes crosslink channel apportionments for all or specified samples]\n";
	echo "         UPLOAD  < sample.def/ratios file > DE       [after LOGRATIOS DE: generates tab delimited file with deadend based protein sample def ratios for upload into XLinkDB]\n";
	echo "         UPLOAD  < sample.def/ratios file > INTRA    [after LOGRATIOS INTRA: generates tab delimited file with intra-link based protein sample def ratios for upload into XLinkDB]\n";
	echo "         UPLOAD  < sample.def/ratios file > DE_INTRA [after LOGRATIOS DE and INTRA: generates tab delimited file with deadend/intra-link based protein sample def ratios for upload into XLinkDB]\n";
	echo "         NONIDENTS                                   [after LOGRATIOS: identifies samples lacking cross-link identifications]\n";
	echo "         COMPARE_DE_XL_SITES                         [after LOGRATIOS and LOGRATIOS DE: compares the number of cross-link sites with deadend quantitation]\n";
	echo "\n";
	echo "  or:    \033[1m{$self} PRINTPARAMS (2-CHANNEL)\033[0m [to print sample file contents to terminal (2-CHANNEL for non-multiplex version)]\n";
	echo "         \033[1m{$self} WRITEPARAMS (2-CHANNEL)\033[0m [to write sample iqpir.params file in current directory (2-CHANNEL for non-multiplex version)]\n";
	if($multichannel) echo "         \033[1m{$self} WRITEPARAMS < /full/path/to/sample1/pepxml,/full/path/to/sample2/pepxml,...> < email > (MULTI_CHANNEL) \033[0m [to write sample iqpir.params file in current directory]\n";
	echo "         \033[1m{$self} DESCRIBE (2-CHANNEL)\033[0m [desribes how sample.map and sample.def are formatted for analysis (2-CHANNEL for non-multiplex version)]\n";
	if($multichannel) echo "         \033[1m{$self} < iqpir.params> WRITEPARAMS_MULTI\033[0m [convert two-channel iqpir.params file to multi-channel iqpir_multi_params]\n";
	echo "\n";
	exit(1);
}

$two_channel = false;
if(count($argv) == 3 && $argv[2] == "2-CHANNEL") {
	array_pop($argv);
	$two_channel = true;
}
if(count($argv)===2) {
	if($argv[1] == "WRITEPARAMS") {
		writeParams(false, ! $two_channel);
	}
	else if($argv[1] == "PRINTPARAMS") {
		writeParams(true, ! $two_channel);
	}
	else if($argv[1] == "DESCRIBE") {
		describeParams(! $two_channel);
	}
	exit(1);
}
else if(count($argv) == 4 && $argv[1] == "WRITEPARAMS") {
	createParamFiles(explode(",", $argv[2]), $argv[3]);
	exit(1);
}
else if(count($argv) == 5 && $argv[1] == "WRITEPARAMS" && $argv[4] == "MULTI_CHANNEL") {
	createParamFiles(explode(",", $argv[2]), $argv[3], true);
	exit(1);
}

if(count($argv)===3 && $argv[2]==="COMBO") {
	$self = "" . $self;
	runCombination($argv[1]);
	exit(1);
}

if(count($argv)===3 && $argv[2]==="AMEND") {
	$iqpir = new MyIqPir($argv[1]); // multiplex
	$newparams = $iqpir->amendHomopeptideCrosslinkQuant();
	if($newparams[1] > 0) {
		$iqpir = new MyIqPir($newparams[0]); // multiplex
		$iqpir->computeCombinedCrosslinkQuant();
		echo $newparams[1] . " cross-link log2ratios amended with {$newparams[0]} ready to upload to XLinkDB\n";
	}
}

// read iqpir.params to see if multiplex or not
$multiplex = false;
$valid = array();
exec("grep ^sample_ratios {$argv[1]}", $valid);
if(count($valid) > 0) {
	$iqpir = new MyIqPirMplex($argv[1]); // multiplex
	$multiplex = true;
}
else {
	if(count($argv)===3 && $argv[2] == "WRITEPARAMS_MULTI") {
		createMultiplexTwoChannelParams($argv[1]);
		exit(1);
	}


	$iqpir = new MyIqPir($argv[1]);
}

if($argv[count($argv)-1]=="USE_CURRENT_NORM") {
	array_pop($argv);
	$iqpir->setUseCurrentNorm();
	echo "Using current normalization factors\n";
}

$deadend = false;
if($argv[count($argv)-1]=="DE") {
	$deadend = true;
	array_pop($argv);
	$iqpir->setDeadendMode();
	echo "Set deadend mode\n"; //exit(1);
}
if(count($argv)===3) {
	if($argv[2]==="RUN") {
		$iqpir->readPepXML();
			$onesample = $iqpir->getOneSampleParamsfile();
			if($onesample !== "") { 
			 	if(! file_exists($onesample)) $iqpir->generateSingleRatioParamFiles(); 
				$iqpir->generateSingeRatioIqpirFile(); 
			}
	}
	else if($argv[2]==="LOGRATIOS") {
		$iqpir->computeCombinedCrosslinkQuant();
		$onesample = $iqpir->getOneSampleParamsfile();
		if($onesample!== "" && file_exists($onesample)) {
			if($multiplex)
				$iqpir_onesample = new MyIqPirMplex($onesample); // multiplex
			else $iqpir_onesample = new MyIqPir($onesample); 
			if($deadend) $iqpir_onesample->setDeadendMode();
			$iqpir_onesample->computeCombinedCrosslinkQuant();
		}
	}
	else if($multiplex && $argv[2] == "UPLOAD") {
		$iqpir->uploadQuantitationToXlinkdbUsingSampleDef('', '');
	}
	else if($argv[2]==="NONIDENTS") {
		$iqpir->getSampleNonIdents();
	}
	else if(false && $argv[2] ==="DEADENDS") {
		$iqpir->quantifyDeadends("/net/gs/vol4/shared/brucelab/search/jdchavez/022920_iqPIR_stdmix_mango/SH_deadend_search/interact-022720_mix1_iqPIR_2h_1.pep.xml", 
		array("mod_mass_info" => "n:198.040247,K:325.127385,K:325.127385,K:327.134095"));
	}
	else if($argv[2]==="COMPOSITE") {
		$iqpir->getCompositeProteinZscores();
	}
	else if($argv[2]==="XL_DE_NORMALIZE") {
		$iqpir->normalizeCrosslinksAndDeadends();
	}
	else if($argv[2]==="DE_INTRA_PROTEINRATIOS") {
		$iqpir->combinedDeadendIntralinkProteinRatios();
	}
	else if($argv[2]==="COMPARE_DE_XL_SITES") {
		$iqpir->compareDeadendAndCrosslinkSites();
	}
}
else if(count($argv)===4) {
	if($argv[2]==="LOGRATIOS") {
		$nextpos = strpos($argv[3], "SAMPLES=");
		if($nextpos===0) {
			$samples = substr($argv[3], $nextpos + strlen("SAMPLES="));
			$samples = explode(",", $samples);
			$iqpir->computeCombinedCrosslinkQuant(array_flip($samples));
		}
		else if($argv[3] == "INTRA") {
			$iqpir->setIntraOnlyMode();
			$iqpir->computeCombinedCrosslinkQuant();
			if($multiplex)
				$iqpir = new MyIqPirMplex($argv[1]); // multiplex
			else
				$iqpir = new MyIqPir($argv[1]);
			$iqpir->combinedDeadendIntralinkProteinRatios(); 

			if(true || $multiplex) {
				$onesample = $iqpir->getOneSampleParamsfile();
				if($onesample!== "" && file_exists($onesample)) {
					if($multiplex)
						$iqpir_onesample = new MyIqPirMplex($onesample); // multiplex
					else $iqpir_onesample = new MyIqPir($onesample); 
					$iqpir_onesample->setIntraOnlyMode();
					$iqpir_onesample->computeCombinedCrosslinkQuant();
					if($multiplex)
						$iqpir_onesample = new MyIqPirMplex($onesample); // multiplex
					else $iqpir_onesample = new MyIqPir($onesample); 
					$iqpir_onesample->combinedDeadendIntralinkProteinRatios(); 
				}
			}

		}
		else if($argv[3] == "ALL") {
			$iqpir->computeCombinedCrosslinkQuant(array(), true);
		}
		else if($argv[3] == "ALLNORM") {
			$iqpir->computeCombinedCrosslinkQuant(array(), true, true);
		}
		else if($argv[3] == "APPORTS") {
			$iqpir->computeApportionments();
			if($multiplex) {
				$onesample = $iqpir->getOneSampleParamsfile();
				if($onesample!== "" && file_exists($onesample)) {
					$iqpir_onesample = new MyIqPirMplex($onesample); // multiplex
					$iqpir_onesample->computeApportionments();
				}
			}
		}
		else if(strpos($argv[3], "APPORTS=")===0) {
			$samples = explode(",", substr($argv[3], strlen("APPORTS=")));
			$iqpir->computeApportionments(array_flip($samples));

		}
		else {
			echo "Error: expected argument SAMPLE=x rather than {$argv[3]}\n";
			exit(1);
		}
	}
	else if($multiplex  && $argv[2]==="RUN") {
		$nextpos = strpos($argv[3], "SAMPLES=");
		if($nextpos===0) {
			$samples = substr($argv[3], $nextpos + strlen("SAMPLES="));
			$samples = explode(",", $samples);
			$iqpir->readPepXML(array_flip($samples));
		}
	}

	else if($argv[2] == "UPLOAD") {
		if($argv[3]=="INTRA") {
			$iqpir->setIntraOnlyMode();
			$iqpir->uploadQuantitationToXlinkdbUsingSampleDef('', '');
		}
		else if($argv[3]=="DE_INTRA") {
			$iqpir->combinedDeadendIntralinkProteinRatios(); 
			$iqpir->setDeadendIntraonlyMode();
			$iqpir->uploadQuantitationToXlinkdbUsingSampleDef('', '');
			if($multiplex) {
				$onesample = $iqpir->getOneSampleParamsfile();
				if($onesample!== "" && file_exists($onesample)) {
					$iqpir_onesample = new MyIqPirMplex($onesample); // multiplex
					$iqpir_onesample->combinedDeadendIntralinkProteinRatios(); 
					$iqpir_onesample->setDeadendIntraonlyMode();
					$iqpir_onesample->uploadQuantitationToXlinkdbUsingSampleDef('', '');
				}
			}
		}
		else if($deadend) {
			if($multiplex) $iqpir->uploadQuantitationToXlinkdbUsingSampleDef('', '');
			else $iqpir->uploadQuantitationToXlinkdbUsingSampleDef('', $argv[3]);
		}
	}
}
else if(count($argv)===5) {
	if($argv[2]==="UPLOAD") {
		if($argv[4]=="INTRA") {
			$iqpir->setIntraOnlyMode();
			$iqpir->uploadQuantitationToXlinkdbUsingSampleDef('', $argv[3]);
		}
		else if($argv[4]=="DE_INTRA") {
			$iqpir->setDeadendIntraonlyMode();
			$iqpir->uploadQuantitationToXlinkdbUsingSampleDef('', $argv[3]);
		}
	}
	else if($argv[2] == "LOGRATIOS" && $argv[3] == "INTRA") {
		$nextpos = strpos($argv[4], "SAMPLES=");
		if($nextpos===0) {
			$samples = substr($argv[4], $nextpos + strlen("SAMPLES="));
			$samples = explode(",", $samples);
			$samples = array_flip($samples);
			echo "Setting intra_only mode\n"; //exit(1);
			$iqpir->setIntraOnlyMode();
			$iqpir->computeCombinedCrosslinkQuant($samples);
			$iqpir = new MyIqPir($argv[1]);
			$iqpir->combinedDeadendIntralinkProteinRatios($samples); 
		}
	}
	else {
		echo "Error: expected argument UPLOAD rather than {$argv[2]}\n";
		exit(1);
	}
}




?>
