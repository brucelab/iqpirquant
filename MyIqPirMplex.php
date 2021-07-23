<?php

include_once("MyPepXML.php");
include_once("MyDeadendPepXML.php");
include_once("MyStatistics.php"); // use this library for the insert statement

class MyIqPirMplex
{

	private $iqpir_paramsfile = NULL;
	private $iqpir_params = NULL; // store the parameters and sample map info
	private $iqpir_ppmtolerance = 20; //25;
	private $ms2_scan = "";
	private $ms3_ppmtol_factor = 20; //20; // how much ppmtol should be increased for ms3 scans
	private $quantitation_type = "iqpir";
	private $filter_based_on_theorrelints_lt = 0.4; //0.4; //0.3; //0; // 0.2 required relative intensity agreement with theor
	private $filter_based_on_theorrelints_hv = 0.4; //0.4; //0.3; //0; // 0.2 required relative intensity agreement with theor
	private $use_optimal_ratios =  true; // 
	private $subtract_noise_from_peaks = false; // whether to reduce observed readmzxml intensities by avearge
	private $min_sig2noise = 2; //2; //3; //2; //5; // 5
	private $memory_limit = '4000M'; //'2000M' //'1000M';
	private $deadend_subdir = "deadend/";
	private $deadend = false;
	private $intra_only_subdir = "intra_only/";
	private $intra_only = false; // whether to restrict logratios to intra-protein cross-links only
	private $use_current_normalization = false;
	private $deadend_intra_subdir = "deadend_intra/";
	private $offsets = array(); // indicate the offset from the lightest mass -> array(stump mod mass, reporter mass)
	private $quantified_isotope_peaks = array(); // how many isotope peaks to quanitify in total
	private $quantified_reporter_peaks = array(); // how many isotope peaks to quanitify in total
	private $quantified_looplink_peaks = array(); // how many isotope peaks to quanitify in total
	private $modmasses_short = array(); //  integer values above the minimum mass zero offset stump
	private $modmasses_long = array(); //  integer values above the minimum mass zero offset stump
	private $sample_offsets = array(); // hashed by sample, indicate the offset from the lightest mass -> array(stump mod mass, reporter mass)
	private $sample_quantified_isotope_peaks = array(); // hashed by sample, how many isotope peaks to quanitify in total
	private $sample_quantified_reporter_peaks = array(); // hashed by sample, how many isotope peaks to quanitify in total
	private $sample_quantified_looplink_peaks = array(); // hashed by sample, how many isotope peaks to quanitify in total
	private $sample_modmasses_short = array(); //  hashed by sample, integer values above the minimum mass zero offset stump
	private $sample_modmasses_long = array(); //  hashed by sample, integer values above the minimum mass zero offset stump
	private $outputdir = "";
	private $verbose = false;
	private $reporter_apportion = false;
	private $num_stumpmods = 1;
	private $num_reporter_pks = 5; //3; //2;
	private $reporter_theor_ints = array(0.759,0.199,0.037,0.004,0.001);
	private $ratio_iontypes = array("LongarmReporter", "Fragment", "Peptide");
	private $analysis_types = array();
	private $peaklist = array();
	private $max_log2ratio = 5;
	private $max_apport_error = 0.1;
	private $max_num_apport_iters = 1000; //10000; // how many times to perturb channel apportionments toward reducing the apportionment error
	private $num_apport_restarts = 5; // number of restarted initial apportionment values to perturb through max_num_apport_iters
	private $quantify_ms2_reporters = false; // whether to spend the time to apportion them and report reuslts to sample iqpir output file file
	private $minprob_with_mincomposite = false; // whether to quantify scans with minprob less than 0.5
	private $dominant_fragid_maxpval = 0.001; // set less than 1 if check whether single frag ion is overabundant for no reason and has ratios that disagree with other frag ion ratios
	private $normalize_with_mean = false;
	private $sample_bioreps = array();
	private $min_apportionment = 0.1; //0.05; //0.025;
	private $channel_modmasses = array();
	private $default_6plex_modmasses = array("826" => "325.127386", "822" => "327.134095", "818" => "329.14080505", "816" => "330.13784005", "812" => "332.14454905", "808" => "334.15125905");
	private $default_6plex_reportermasses = array("826" => "825.490275", "822" => "821.47686", "818" => "817.463", "816" => "815.469", "812" => "811.455947", "808" => "807.442528");
	
	public function __construct ($iqpir_paramsfile) { //$samplemap, $ml_filename, $workingdir, $masschroq_outputfile, $quantitation_type="silac") {
		if(! file_exists($iqpir_paramsfile)) {
			echo "Error: iqpir params file {$iqpir_paramsfile} does not exist\n";
			exit(1);
		}
		$this->iqpir_paramsfile = $iqpir_paramsfile;
		$this->readIqPirParams();
		$this->outputdir = $this->iqpir_params['iqpir_output_dir'];
		
		if(! $this->quantify_ms2_reporters && array_key_exists("quantify_ms2_reporters", $this->iqpir_params) && 
			$this->iqpir_params["quantify_ms2_reporters"] == "yes") $this->quantify_ms2_reporters = true;
		else if($this->quantify_ms2_reporters && array_key_exists("quantify_ms2_reporters", $this->iqpir_params) && 
			$this->iqpir_params["quantify_ms2_reporters"] != "yes") $this->quantify_ms2_reporters = false;
		ini_set('memory_limit', $this->memory_limit);	
		// read all mod masses, set to offsets
		if(count($this->reporter_theor_ints) > 0) {
			if($this->num_reporter_pks < count($this->reporter_theor_ints)) {
				while(count($this->reporter_theor_ints) > $this->num_reporter_pks) array_pop($this->reporter_theor_ints);
				$tot = array_sum($this->reporter_theor_ints);
				MyStatistics::normalize($this->reporter_theor_ints, 3);
			}
			echo "Using reporter theoretical ints ".join(",", $this->reporter_theor_ints)."\n";
		}
		if(count($this->ratio_iontypes) > 0) {
			$this->ratio_iontypes = array_flip($this->ratio_iontypes);
			ksort($this->ratio_iontypes);
		}
		else {
			echo "Error 1: have no ratio iontypes with which to compute ratios\n";
			exit(1);
		}
	}
	
	
	private function setAnalysisTypes() {
		// grep from file
		if(array_key_exists("LongarmReporter", $this->ratio_iontypes)) {
			if(! array_key_exists("longarm_reporters", $this->iqpir_params) || $this->iqpir_params['longarm_reporters'] != "true") {
				unset($this->ratio_iontypes['LongarmReporter']);
			}
			else if($this->deadend) unset($this->ratio_iontypes['LongarmReporter']);
			else {		
				$file = $this->iqpir_params['xlinkprophetfile'];		

				$valid = array();
				exec("grep '<search_summary' {$file}", $valid);
				for($k = 0; $k < count($valid); $k++) {
					if(strpos($valid[$k], "SpectraST")!==false) $this->analysis_types['SpectraST'] = 1;
					else if(strpos($valid[$k], "Mango")!==false) $this->analysis_types['Mango'] = 1;
					else $this->analysis_types['ReACT'] = 1;
				}
				if(! array_key_exists("ReACT", $this->analysis_types)) unset($this->ratio_iontypes['LongarmReporter']);
			}
		}
		if(count($this->ratio_iontypes) == 0) {
			echo "Error 2: have no ratio iontypes with which to compute ratios\n";
			exit(1);
		}
	}
	
	private function setSampleChannels($sample) {
		$this->offsets = $this->sample_offsets[$sample];
		$this->quantified_isotope_peaks = $this->sample_quantified_isotope_peaks[$sample];
		$this->quantified_reporter_peaks = $this->sample_quantified_reporter_peaks[$sample];
		$this->quantified_looplink_peaks = $this->sample_quantified_looplink_peaks[$sample];
		$this->modmasses_short = $this->sample_modmasses_short[$sample];  
		$this->modmasses_long = $this->sample_modmasses_long[$sample]; 
	}
	
	private function setChannelOffsets($sample, $stumpmasses, $reportermasses) {
		if(count($stumpmasses)==0 || count($stumpmasses) != count($reportermasses)) {
			echo "Error with ".count($stumpmasses)." stumpmasses and ".count($reportermasses)." reportermasses\n";
			exit(1);
		}
		$this->sample_offsets[$sample] = array();
		$this->sample_modmasses_short[$sample] = array();
		$this->sample_modmasses_long[$sample] = array();
		$this->sample_quantified_isotope_peaks[$sample] = array();
		$this->sample_quantified_reporter_peaks[$sample] = array();
		$this->sample_quantified_looplink_peaks[$sample] = array();
		$min_index = -1;
		$min = 1000000;
		for($k = 0; $k < count($stumpmasses); $k++) {
			if($stumpmasses[$k] < $min) {
				$min = $stumpmasses[$k];
				$min_index = $k;
			}
		}
		$masstot = 0;
		for($k = 0; $k < count($stumpmasses); $k++) {
			if($masstot == 0) $masstot = 2*$stumpmasses[$k] + $reportermasses[$k];
			else if(abs(2*$stumpmasses[$k] + $reportermasses[$k] - $masstot) > 0.1) {
				echo "Error: have paired stumpmasses {$stumpmasses[$k]} and reportermasses {$reportermasses[$k]} that don't sum to the same total mass {$masstot}\n";
				exit(1);
			}
		
			$next_offset = strval(number_format($stumpmasses[$k] - $min, 0, ".", ""));
			if(array_key_exists($next_offset, $this->sample_offsets[$sample])) {
				echo "Error: have {$stumpmasses[$k]} with redundant offset {$next_offset} among ".join(",", $stumpmasses)."\n";
				exit(1);
			}
			$this->sample_offsets[$sample][$next_offset] = array($stumpmasses[$k], $reportermasses[$k]);
			// these used to tell if peptide is heavy or light, and used to convert heavy to light
			$this->sample_modmasses_short[$sample][$next_offset] = number_format($stumpmasses[$k], 0, ".", "");
			$this->sample_modmasses_long[$sample][$next_offset] = number_format($stumpmasses[$k], 2, ".", "");
			
			$this->channel_modmasses[number_format($reportermasses[$k]+0.9, 0, ".", "")] = $stumpmasses[$k];
		}		
		ksort($this->sample_offsets[$sample]);
		echo "Channel offsets: ".join(",", array_keys($this->sample_offsets[$sample]))." for sample {$sample}\n";		
		
		foreach($this->sample_offsets[$sample] as $offset => $info) {
			echo "Next offset {$offset} with values ".join(",", $info)."\n";
			for($k = 0; $k < 5; $k++) {
				$this->sample_quantified_isotope_peaks[$sample][$offset + $k] = 1;
				if($k < $this->num_reporter_pks) $this->sample_quantified_reporter_peaks[$sample][-2 * $offset + $k] = 1;
				$this->sample_quantified_looplink_peaks[$sample][2 * $offset + $k] = 1;
			}
		}
		ksort($this->sample_quantified_isotope_peaks[$sample]);
		ksort($this->sample_quantified_reporter_peaks[$sample]);
		ksort($this->sample_quantified_looplink_peaks[$sample]);
		echo "Quantified isotope peaks: ".join(",", array_keys($this->sample_quantified_isotope_peaks[$sample]))."\n";
		echo "Quantified reporter peaks: ".join(",", array_keys($this->sample_quantified_reporter_peaks[$sample]))."\n";
		echo "Quantified looplink peaks: ".join(",", array_keys($this->sample_quantified_looplink_peaks[$sample]))."\n";
		// these used with getIsotopeQuant to quanitfy only nec isotopes to compare with expected
		echo "Have ".count($this->sample_quantified_isotope_peaks[$sample])." isotope peaks to quanitify for each ion relative to {$min} monoisotopic\n";
		echo "Short Modmasses  ".join(",", $this->sample_modmasses_short[$sample])."\n"; // {$this->light_modmass_long}\n";
		echo "Long modmasses ".join(",", $this->sample_modmasses_long[$sample])."\n"; // ".join(",", $this->heavy_modmasses_long)."\n";
		foreach($this->channel_modmasses as $channel => $stump) echo "channel\t" . $channel . "\t" . $stump . "\n";
	}

	
	public function setUseCurrentNorm() {
		$this->use_current_normalization = true;
		$normfile = $this->iqpir_params['iqpir_output_dir'] . "sample_biorep_normfactors.txt";
		if(! file_exists($normfile)) {
			echo "Norm file {$normfile} does not exist, you cannot set to use current normalization\n";
			exit(1);
		}
		$myfile = fopen($normfile, "r");
		$first = true;
		while($line = rtrim(fgets($myfile))) {
			if($first) {
				$first = false;
				continue;
			}
			$next = split("\t", $line);
			$this->iqpir_params['samples'][$next[0]][$next[1] . "-" . $next[2] . "-" . $next[3]] = $next[4];
		}
		fclose($myfile);
		echo "sample\tbiorep\tnorm_add_to_log\n";
		foreach($this->iqpir_params['samples'] as $sample => $bioreps) {
			foreach($bioreps as $biorep => $norm) {
				echo $sample . "\t" . $biorep . "\t" . $norm ."\n";
			}
		}
	}
	
	public function setTraditionalParams() {
		$this->iqpir_ppmtolerance = 25;
		$this->filter_based_on_theorrelints_lt = 0.3;
		$this->filter_based_on_theorrelints_hv = 0.3;
		$this->use_optimal_ratios = false;
		$this->subtract_noise_from_peaks = true;
		$this->deadend_subdir = "deadend_trad/";
		echo "Using traditional parameter values\n";
	}

	public function setDeadendMode() {
		if($this->deadend) return; // already in deadend mode
		$this->iqpir_params['iqpir_output_dir'] .= $this->deadend_subdir;
		if($this->iqpir_params['iqpir_output_dir'][strlen($this->iqpir_params['iqpir_output_dir'])-1]!="/") $this->iqpir_params['iqpir_output_dir'] .= "/";
		if(! file_exists($this->iqpir_params['iqpir_output_dir'])) system("mkdir {$this->iqpir_params['iqpir_output_dir']}");
		$this->deadend = true;
	}

	public function setIntraOnlyMode() {
		if($this->intra_only) return;
		if($this->deadend) {
			echo "Error: intra_only mode incompatible with deadend mode\n";
			exit(1);
		} 
		$this->iqpir_params['iqpir_output_dir'] .= $this->intra_only_subdir;
		if($this->iqpir_params['iqpir_output_dir'][strlen($this->iqpir_params['iqpir_output_dir'])-1]!="/") $this->iqpir_params['iqpir_output_dir'] .= "/";
		if(! file_exists($this->iqpir_params['iqpir_output_dir'])) system("mkdir {$this->iqpir_params['iqpir_output_dir']}");
		$this->intra_only = true;
	
	}
	
	public function setDeadendIntraonlyMode() {
		$this->intra_only = true;
		$this->deadend = true;
		$this->iqpir_params['iqpir_output_dir'] .= $this->deadend_intra_subdir;
	}
	
	
private function recordMs3Offsets($modpep, $ms3scan, $mzxml, &$header, $fout = STDERR) {
	$pep1_nostumpfrags = $this->getStumpFragments($modpep, 500, 2500, false); //325.13); //exit(1);
	if($header) {
		fwrite($fout, "mzxml\tscan\tpeptide\tion\tacutal_mz\tobs_maxmz\tobs_mz_maxint\tobs_minus_actual_mz\ttheor_relintens\n");
		$header = false;
	}
	$next_mzxml = $mzxml;
	$nextpos = strrpos($next_mzxml, "/");
	if($nextpos !== false) $next_mzxml = substr($next_mzxml, $nextpos + 1);
	foreach($pep1_nostumpfrags as $key => $value) {
		$next_maxmz = $this->getIsotopeQuant($value[0], 1, $ms3scan, $mzxml, 0, "", 1, 1);
		$next_ints_fragentries = $this->getIsotopeIntensities($value[1], true);
		if($next_maxmz[1] > -1) fwrite($fout, $next_mzxml . "\t" . $ms3scan . "\t" . $modpep . "\t" . $key . "\t" . 
			join("\t", $next_maxmz)."\t" . ($next_maxmz[1] - $next_maxmz[0]) . "\t" . join(",", $next_ints_fragentries)."\n"); 
	}
}

private function getStumpFragments($peptide, $min_mass = 500, $max_mass = 2500, $include_stump = true) {
	$aa_masses = array("A" => 71.03711, "R" => 156.10111, "N" => 114.04293, "D" => 115.02694, "C" => 103.00918, "E" => 129.04259, "Q" => 128.05858,
						"G" => 57.02146, "H" => 137.05891, "I" => 113.08406, "L" => 113.08406, "K" => 128.09496, "M" => 131.04049, "F" => 147.06841,
						"P" => 97.05276, "S" => 87.03203, "T" => 101.04768, "W" => 186.07931, "Y" => 163.06333, "V" => 99.06841);
						
	// go through the peptide counting residues for both b and y ions
	$b_index = 0;
	$stump_res = -1;
	$stump_res2 = -1;
	$hydrogen_mass = 1.007825035;
	$carbon_mass = 12.000000;
	$oxygen_mass = 15.99491463;
	$proton_mass = 1.00727647;
	$nitrogen_mass = 14.003074;
	
	$stump_mod = $this->offsets[0][0];
	$modmasses = array($this->modmasses_long[0] => $stump_mod, "147.04" => 147.035385);
	
	$b_mass = 0; //$nitrogen_mass + 2 * $hydrogen_mass;
	$y_mass = $oxygen_mass + 2 * $hydrogen_mass;

	$mods = array(); // record by index
	$output = array(); // hash ion to mass
	$last = "";
	$index_first = array();
	$last_mass = 0;
	for($k = 0; $k < strlen($peptide); $k++) {
		$b_index++;
		$b_mass += $aa_masses[$peptide[$k]];
		$last_mass = $aa_masses[$peptide[$k]];
		$index_first[$b_index] = $k;
		if($k < strlen($peptide) - 1 && $peptide[$k+1]=='[') { // must get the mass and move on
			$nextpos = strpos(substr($peptide, $k+2), "]");
			if($nextpos!==false) {
				$next_mass = substr($peptide, $k+2, $nextpos);
				if(array_key_exists($next_mass, $modmasses)) {
					$next_mass = $modmasses[$next_mass]; // get the more accurate value
				}
				$b_mass += $next_mass - $last_mass;
				$mods[$b_index-1]=$next_mass - $last_mass;
				$k = $nextpos + $k + 2;
				if($stump_res == -1 && $next_mass == $stump_mod) {
					$stump_res = $b_index;
				}
				else if($stump_res2 == -1 && $next_mass == $stump_mod) {
					$stump_res2 = $b_index;
				}
			}
			else {
				echo "Error with no ] in {$peptide}\n";
				exit(1);
			}
		}
		if($b_mass < $min_mass || $b_mass > $max_mass) continue;
		if($include_stump && $stump_res > -1 && $b_index >= $stump_res) {
			$output['b' . $b_index] = array($b_mass, substr($peptide, 0, $k + 1));
			$last = 'b' . $b_index;
		}
		if(! $include_stump && ($stump_res == -1 || $b_index < $stump_res)) {
			$output['b' . $b_index] = array($b_mass, substr($peptide, 0, $k + 1));
			$last = 'b' . $b_index;
		}
	
	}
	$stump_res2 = max($stump_res, $stump_res2); // for y ions
	if($last != "") unset($output[$last]); // remove the full length b ion
	$stripped = MyPepXML::stripMods($peptide); // now do the y ions
	$y_index = 0;
	for($k = strlen($stripped)-1; $k > 0; $k--) {
		$y_index++;
		$y_mass += $aa_masses[$stripped[$k]];
		if(array_key_exists($k, $mods)) $y_mass += $mods[$k];
		if($y_mass < $min_mass || $y_mass > $max_mass) continue;
		if($include_stump && $k < $stump_res2) $output['y' . $y_index] = array($y_mass, substr($peptide,  $index_first[strlen($stripped) + 1 - $y_index]));
		else if(! $include_stump && $k >= $stump_res2) $output['y' . $y_index] = array($y_mass, substr($peptide,  $index_first[strlen($stripped) + 1 - $y_index]));
	}
	return $output;
}


private function readIqPirParams() { // iqpir.params and reference sample_map
$this->iqpir_params = array();
$myfile = fopen($this->iqpir_paramsfile, "r");
while($line=rtrim(fgets($myfile))){
	$next = split("\t", $line);
	if(count($next) > 1 && $next[1] != "") $this->iqpir_params[$next[0]] = $next[1];
}
fclose($myfile);
if(! array_key_exists('sample_map', $this->iqpir_params)) {
	echo "Error with missing sample_map in {$this->iqpir_paramsfile} file: have only ".join(",", array_keys($this->iqpir_params))."\n";
	exit(1);
}
if(! file_exists($this->iqpir_params['sample_map'])) {
	echo "Error sample_map file {$this->iqpir_params['sample_map']} doesn not exist\n";
	exit(1);
}
if(! array_key_exists('fdr', $this->iqpir_params)) {
	echo "Error: no fdr parameter set in {$this->iqpir_paramsfile}\n";
	exit(1);
}
foreach($this->iqpir_params as $key => $value) echo $key . "\t" . $value . "\n";
$myfile = fopen($this->iqpir_params['sample_map'], "r");
$this->iqpir_params['rawfiles'] = array(); // rawfile name to array of fwd/rev, sample, and biorep
$this->iqpir_params['samples'] = array(); // sample number to array of biorep numbers
while($line=rtrim(fgets($myfile))){
	$next = split("\t", $line);
	if(array_key_exists($next[0], $this->iqpir_params['rawfiles'])) {
		echo "Error: have two copies of rawfile {$next[0]} in sample map {$this->iqpir_params['sample_map']}\n";
		exit(1);
	}
	$this->iqpir_params['rawfiles'][$next[0]] = array(
		"sample" => $next[1], "biorep" => $next[2]);
	if(! array_key_exists($next[1], $this->iqpir_params['samples'])) $this->iqpir_params['samples'][$next[1]] = array();
	if(! array_key_exists($next[1], $this->sample_bioreps)) $this->sample_bioreps[$next[1]] = array();
	if(! array_key_exists($next[2], $this->sample_bioreps[$next[1]])) $this->sample_bioreps[$next[1]][$next[2]] = 1;
}
fclose($myfile);
ksort($this->iqpir_params['samples']);
if(array_key_exists("sample_channels", $this->iqpir_params)) {
	$myfile = fopen($this->iqpir_params['sample_channels'], "r");
	while($line=rtrim(fgets($myfile))){
		if(line == "") continue; // empty line
		$next = split("\t", $line);
		$this->setChannelOffsets($next[0], split(",", $next[1]), split(",", $next[2]));
	}
	fclose($myfile);
}
else { // use default
	foreach($this->iqpir_params['samples'] as $sample => $info) {
		$this->setChannelOffsets($sample, array_values($this->default_6plex_modmasses), array_values($this->default_6plex_reportermasses));
	}
}
// now make sure have offsets for each sample
foreach($this->iqpir_params['samples'] as $sample => $bioreps) {	
	if(! array_key_exists($sample, $this->sample_offsets)) {
		echo "Error: have no stumpmodmass info for sample {$sample}\n";
		echo "Please make sure your sample_channels file includes all samples in your sample_map\n";
		exit(1);
	}
}


foreach($this->iqpir_params['rawfiles'] as $file => $info) echo $file . "\t" . join(",", array_keys($info))."\t".join(",", array_values($info))."\n";
foreach($this->iqpir_params['samples'] as $sample => $bioreps) echo $sample . "\t" . join(",", array_keys($bioreps))."\n";
if($this->iqpir_params['iqpir_output_dir'][strlen($this->iqpir_params['iqpir_output_dir'])-1]!="/") $this->iqpir_params['iqpir_output_dir'] .= "/";
if(! file_exists($this->iqpir_params['iqpir_output_dir'])) system("mkdir {$this->iqpir_params['iqpir_output_dir']}");

}

public function computeRatioDeconvError($theor_relints, $obs_relints, $ratio) {
	$tot = 0;
	$z = $theor_relints[0]+$theor_relints[1]+$theor_relints[2];
	for($k = 0; $k < 5; $k++) {
		$next = $theor_relints[$k] * $ratio;
		if($k > 1) $next += $theor_relints[$k-2];
		$next /= ($z + $ratio);
		$tot += ($next - $obs_relints[$k]) * ($next - $obs_relints[$k]);
	}
	$tot /= 5;
	return sprintf("%0.2e", sqrt($tot));	
}
// to allow for multiple databases
private function setDatabaseProteinSequencesAndGenes($fasta, &$output) {
$current = "";
$init = count($output);
$myfile = fopen($fasta, "r");
while($line=rtrim(fgets($myfile))){
	if($line[0]==">") {
		$current = "";
		if(strpos($line, "rev_")===false) {
			$current = substr($line, 1);
			$current = split('\|', $current);
			if(array_key_exists($current[1], $output)) continue;
			$nextpos = strpos($current[2], " ");
			if($nextpos!==false) $current[2] = substr($current[2], 0, $nextpos);
			$output[$current[1]] = array("gene" => $current[2], "sequence" => "");
			$current = $current[1];
		}
	}
	else if($current != "") {
		$output[$current]['sequence'] .= preg_replace("/L/", "I", $line);
	}
}
fclose($myfile);
echo "Read in ".(count($output) - $init)." protein sequences from {$fasta}\n"; //exit(1);
}
private function getDatabaseProteinSequencesAndGenes($fasta) {
$current = "";
$output = array();
$myfile = fopen($fasta, "r");
while($line=rtrim(fgets($myfile))){
	if($line[0]==">") {
		$current = "";
		if(strpos($line, "rev_")===false) {
			$current = substr($line, 1);
			$current = split('\|', $current);
			$nextpos = strpos($current[2], " ");
			if($nextpos!==false) $current[2] = substr($current[2], 0, $nextpos);
			$output[$current[1]] = array("gene" => $current[2], "sequence" => "");
			$current = $current[1];
		}
	}
	else if($current != "") {
		$output[$current]['sequence'] .= preg_replace("/L/", "I", $line);
	}
}
fclose($myfile);
echo "Read in ".count($output)." protein sequences from {$fasta}\n"; //exit(1);
return $output; //array($seqs, $genes);
}

private function setScanPeaklist($mzxml, $ms3scans, $longarm_scans) {
	$this->peaklist = array($this->ms2_scan => array());
	for($k = 0; $k < count($ms3scans); $k++) {
		$this->peaklist[$ms3scans[$k]] = array();
	}
	for($k = 0; $k < count($longarm_scans); $k++) {
		$this->peaklist[$longarm_scans[$k]] = array();
	}
	foreach($this->peaklist as $scan => $info) {
		$command = "readmzXML {$mzxml} ";
		$valid = array();
		exec($command . $scan, $valid);
		$done = false;
		for($k = 0; $k < count($valid); $k++) {
			if(strpos($valid[$k], "mass")===false) continue;
			$next = preg_split('/\s+/', $valid[$k]);
			$num_cols = count($next);
			$mass_col = count($next) - 3; //1;
			$int_col = count($next) - 1; //3;
			$this->peaklist[$scan][$next[$mass_col]] = $next[$int_col];
		}
	}
}

private function getIsotopeQuant($mass, $charge, $scan, $mzxml, $noise = 0, $title = "", $num_isotopes = 3, $ppm_diff = "") {
if($charge == "") echo "Here with zero charge for mass {$mass} and scan {$scan} in {$mzxml}\n";
	$proton_mass = 1.00727647;
	//$mass += 0.984016;
	$c13_c12_massdiff = 1.003355;
	$ppmtolerance = $this->iqpir_ppmtolerance; //25; //10; //100; //100; //10;
	if($scan != $this->ms2_scan) $ppmtolerance *= $this->ms3_ppmtol_factor;
	$nextmz = ($mass + $charge * $proton_mass)/$charge;
	$pepmzs = array();
	$max_mzs = array();
	for($k = 0; $k < $num_isotopes; $k++) {
		array_push($pepmzs, $nextmz + $k * $c13_c12_massdiff / $charge);
		array_push($max_mzs, array(-1, 0)); // m/z, intensity
	}
	$pepints = array();
	for($k = 0; $k < $num_isotopes; $k++) {
		array_push($pepints, 0);
	}
	$verbose = false; 
	$command = "readmzXML {$mzxml} ";
	if($verbose) echo "readmzXML {$mzxml} {$scan}\n";
	$output_maxmz = false;
	if($ppm_diff == "") $ppm_diff = $nextmz * $ppmtolerance / 1000000;
	else $output_maxmz = true;
	if($verbose) echo "PPM {$nextmz} += {$ppm_diff}\n";
	exec($command . $scan, $valid);
	$done = false;
	for($k = 0; $k < count($valid); $k++) {
		if(strpos($valid[$k], "mass")===false) continue;
		$next = preg_split('/\s+/', $valid[$k]);

		$num_cols = count($next);
		$mass_col = count($next) - 3; //1;
		$int_col = count($next) - 1; //3;

		for($j = 0; $j < count($pepmzs); $j++) {
			if($verbose) echo "Comparing ".$next[$mass_col]." with ".($pepmzs[$j] - $ppm_diff)." and ".($pepmzs[$j] + $ppm_diff)."\n";
			if($next[$mass_col] >= $pepmzs[$j] - $ppm_diff && $next[$mass_col] <= $pepmzs[$j] + $ppm_diff) {
				if($this->subtract_noise_from_peaks) $next[$int_col] = max(0, $next[$int_col]-$noise);
				$pepints[$j] += $next[$int_col];
				if($verbose) echo "adding {$next[$int_col]} for {$pepmzs[$j]}\n";
				if($next[$int_col] > $max_mzs[$j][1]) {
					$max_mzs[$j][0] = $next[$mass_col];
					$max_mzs[$j][1] = $next[$int_col];
				}
			}
			else if($j==count($pepmzs)-1 && $next[$mass_col] > $pepmzs[$j] + $ppm_diff) { // done
					if($verbose) echo "Finished here with {$valid[$k]} since beyond ".($pepmzs[$j] + $ppm_diff)."\n";
				$k = count($valid);
			}
		} 
	}
	$sig2noise = array();
	$ave_mz_diff = 0;
	$num_mz_diffs = 0;
	for($k = 0; $k < count($pepints); $k++) {
		$pepints[$k] = number_format($pepints[$k], 0, ".", "");
		if($max_mzs[$k][0] > -1) {
			$ave_mz_diff += ($max_mzs[$k][0] - $pepmzs[$k]) * 1000000/$pepmzs[$k];
			$num_mz_diffs++;
		}
	}
	if($verbose) echo "HERE with PEPINTS: ".join(",", $pepints)." for mz ".join(",", $pepmzs)." and scan {$scan} with masstol {$ppm_diff}\n";
	if($verbose) 
		echo $title . join(",", $pepints)." with sig2noise ".join(",", $sig2noise)." for {$mass} {$charge} in scan {$scan} of {$mzxml} and mzs ".join(",", $pepmzs)."\n";
	// now compute the average mass difference
	if($output_maxmz) return array($pepmzs[0], $max_mzs[0][0], number_format($max_mzs[0][1], 0, ".", ""));
	return array($pepints, $sig2noise, $pepmzs, $num_mz_diffs > 0 ? number_format($ave_mz_diff/$num_mz_diffs, 2, ".", "") : "N/A");

}
private function  getExpectedObs_withOffsets($theor, $ch_totals) {
$output = array();
foreach($ch_totals as $offset => $value) {
	for($k = 0; $k < 5; $k++) {
		if($this->reporter_apportion && $k >= $this->num_reporter_pks) continue;
		if(! array_key_exists($offset + $k, $output)) $output[$offset + $k] = 0;
		$output[$offset + $k] += $theor[$k] * $value;
	}
}
return $output;
}
private function sortByError($a, $b) {
	if($a[1] < $b[1]) return -1;
	if($a[1] > $b[1]) return 1;
	return 0;
}

private function minimizeError_withOffsets($theor, $obs, &$ch_totals) {
if(array_sum(array_values(($ch_totals)))==0) {
	$total = 0;
	foreach($ch_totals as $offset => $value) {
		$ch_totals[$offset] = ($obs[$offset]+$obs[$offset+1])/($theor[0]+$theor[1]);
		$total += $ch_totals[$offset];
	}
	if($total > 0) {
		foreach($ch_totals as $offset => $value) $ch_totals[$offset] = number_format($ch_totals[$offset] / $total, 5, ".", "");
	}
}
$expected = $this->getExpectedObs_withOffsets($theor, $ch_totals);
$error = 0;
foreach($expected as $peak => $value) $error += (($obs[$peak] - $value)*($obs[$peak] - $value));
return sqrt($error);
}

private function perturb_withOffsets($vals, $factor = 0.9, $add_factor = 0.05) {
$output = $vals;
$indeces = array_keys($vals);
$index = rand(0, count($indeces)-1);
$type = rand(0,1) < 1 ? "add" : "mult";
$index = $indeces[$index];
if($type == "add") {
	rand(0, 1) < 1 ? $output[$index] +- $add_factor : $output[$index] += $factor;
	if($output[$index] < 0) $output[$index] = 0;


}
else {
	rand(0, 1) < 1 ? $output[$index] /= $factor : $output[$index] *= $factor;
}
$total = 0;
foreach($output as $offset => $val) $total += $val;
if($total > 0) {
	foreach($output as $offset => $val) $output[$offset] /= $total;
}
return $output;


}

private function getOptimalChannelIntensities_withOffsets($theor, $obs, $ch_totals) {  
$curr_error = $this->minimizeError_withOffsets($theor, $obs, $ch_totals);
for($k = 0; $k < $this->max_num_apport_iters; $k++) {
	$next = $this->perturb_withOffsets($ch_totals, 0.99);
	$error = $this->minimizeError_withOffsets($theor, $obs, $next);
	if($error < $curr_error) {
		$ch_totals = array();
		foreach($next as $offset => $val) $ch_totals[$offset] = number_format($val, 5, ".", "");
		$curr_error = $error;
		if($curr_error == 0) $k = $this->max_num_apport_iters;
	}
}
return array($ch_totals, number_format($curr_error, 4, ".", ""));
}

public function multiplex_withOffsets($theor, $obs) {

$results = array();
$offsets = array_keys($this->offsets);
$ch_totals = array();
$factor = 1;
if($this->reporter_apportion) $factor = -2;
else if($this->num_stumpmods == 2) $factor = 2;
foreach($this->offsets as $offset => $val) $ch_totals[$factor * $offset] = 0;
if($this->reporter_apportion) ksort($ch_totals);
array_push($results, $this->getOptimalChannelIntensities_withOffsets($theor, $obs, $ch_totals));
$num_chans = count($ch_totals);
for($k = 0; $k < $this->num_apport_restarts; $k++) {
	$ch_totals = array();
	for($j = 0; $j < count($offsets); $j++) $ch_totals[$factor * $offsets[$j]] = 0; //$ch_totals_next[$k];
	if($this->reporter_apportion) ksort($ch_totals);

	$this->minimizeError_withOffsets($theor, $obs, $ch_totals); // this to set ch_totals to equal values for all channels
	array_push($results, $this->getOptimalChannelIntensities_withOffsets($theor, $obs, $ch_totals));
}

usort($results, array($this, "sortByError"));
return $results[0];

}

private function getIsotopeQuant_withOffsets($mass, $charge, $scan, $mzxml, $noise = 0, $title = "") {
if($charge == "") echo "Here with zero charge for mass {$mass} and scan {$scan} in {$mzxml}\n";
	$proton_mass = 1.00727647;
	$c13_c12_massdiff = 1.003355;
	$ppmtolerance = $this->iqpir_ppmtolerance; //25; //10; //100; //100; //10;
	if($scan != $this->ms2_scan) $ppmtolerance *= $this->ms3_ppmtol_factor;
	$nextmz = ($mass + $charge * $proton_mass)/$charge;
	$pepmzs = array();
	$max_mzs = array();
	$pepints = array();
	$final_offset = -1;
	$quantified_peaks = $this->quantified_isotope_peaks;
	if($this->reporter_apportion) $quantified_peaks = $this->quantified_reporter_peaks;
	else if($this->num_stumpmods == 2) $quantified_peaks = $this->quantified_looplink_peaks;
	foreach($quantified_peaks as $offset => $val) {
		$pepmzs[$offset] = $nextmz + $offset * $c13_c12_massdiff / $charge;
		$max_mzs[$offset] = array(-1, 0); // m/z, intensity
		$pepints[$offset] = 0;
		$final_offset = $offset;
	}
	$verbose = false; 
	$command = "readmzXML {$mzxml} ";
	if($verbose) echo "readmzXML {$mzxml} {$scan}\n";
	$ppm_diff = $nextmz * $ppmtolerance / 1000000;
	if($verbose) echo "PPM {$nextmz} += {$ppm_diff}\n";
	exec($command . $scan, $valid);
	$done = false;
	foreach($this->peaklist[$scan] as $next_mz => $next_int) {
		if($done) continue;
		foreach($pepmzs as $offset => $mz) {
			if($next_mz >= $mz - $ppm_diff && $next_mz <= $mz + $ppm_diff) {
				if($this->subtract_noise_from_peaks) $next_int = max(0, $next_int-$noise);
				$pepints[$offset] += $next_int;
				if($verbose) echo "adding {$next[$int_col]} for {$mz}\n";
				if($next_int > $max_mzs[$offset][1]) {
					$max_mzs[$offset][0] = $next_mz;
					$max_mzs[$offset][1] = $next_int;
				}
			}
			else if($offset == $final_offset && $next_mz > $mz + $ppm_diff) { // done
					if($verbose) echo "Finished here offset {$offset} with $next_mz since beyond ".($mz + $ppm_diff)."\n";
				$done = true;
			}
		} 
	}
	$sig2noise = array();
	$ave_mz_diff = 0;
	$num_mz_diffs = 0;
	foreach($pepints as $offset => $tot) {
		$pepints[$offset] = number_format($pepints[$offset], 0, ".", "");
		if($max_mzs[$offset][0] > -1) {
			$ave_mz_diff += ($max_mzs[$offset][0] - $pepmzs[$offset]) * 1000000/$pepmzs[$offset];
			$num_mz_diffs++;
		}
	}
	$total_intensity = array_sum(array_values($pepints));
	if($total_intensity > 0) {
		foreach($pepints as $offset => $tot) $pepints[$offset] = number_format($pepints[$offset]/$total_intensity, 5, ".", "");
	}
	
	return array($pepints, $total_intensity, $pepmzs, $num_mz_diffs > 0 ? number_format($ave_mz_diff/$num_mz_diffs, 2, ".", "") : "N/A");

}



private function printRatios($light, $heavy, $title = "", $light2heavy = true) {
	$ratios = array();
	# want to get a, c, and e
	
	for($k = 0; $k < count($light);  $k++) {
		
		if(! $light2heavy && $light[$k] > 0) {
			array_push($ratios, number_format($heavy[$k]/$light[$k], 2));
		}
		else if($light2heavy && $heavy[$k] > 0) {
			array_push($ratios, number_format($light[$k]/$heavy[$k], 2));
		}
		else {
			array_push($ratios, "N/A");
		}
	}
	return $ratios;
}


private function getIsotopeIntensities($seq_or_mass, $normalize = false, $num_isotopes = 5) {
$command = "calcisotopes ";
if(is_numeric($seq_or_mass)) {
	$command .= $seq_or_mass;
}
else {
	$command .= "\"{$seq_or_mass}\"";
}
$valid = array();
exec($command, $valid);
# now parse the reuslts
$output = array();
$tot = 0; // calculate total to make ratios
for($k = 1; $k < count($valid); $k++) {
	if($num_isotopes > 0 && $k > $num_isotopes) continue;
	$next = preg_split('/\s+/', $valid[$k]); //split("\t", $valid[$k]);
	array_push($output, $next[1]);
	if($normalize) $tot += $next[1];
}
while($num_isotopes > 0 && count($output) < $num_isotopes) array_push($output, 0);
// normalize
if($normalize && $tot > 0) {
	for($k = 0; $k < count($output); $k++) {
		$output[$k] = number_format($output[$k]/$tot, 3);
	}
}
return $output;
}

private function sortBySecondThenThirdIndex($a, $b) {
	if($a[1] < $b[1]) return -1;
	if($a[1] > $b[1]) return 1;
	if($a[2] < $b[2]) return -1;
	if($a[2] > $b[2]) return 1;
	return 0;
}


private function getChannelOffset($modpep) {
	foreach($this->modmasses_short as $offset => $mass) {
		if(strpos($modpep, "K[{$mass}")!==false) return $offset;
	}
	return -1;
	echo "Error: no offset found for {$modpep} versus ".join(",", array_keys($this->modmasses_short))."\n";
	exit(1);	
}
private function convertHeavyToLight($modpep) {
	$output = $modpep;
	foreach($this->modmasses_long as $offset => $mass) {
		if($offset == 0) continue;
		$output = preg_replace("/K\[{$mass}\]/", "K[{$this->modmasses_long[0]}]", $output);
	}
	if($this->deadend) {
		foreach($this->modmasses_short as $offset => $mass) {
			if($offset == 0) continue;
			$output = preg_replace("/K\[{$mass}\]/", "K[{$this->modmasses_short[0]}]", $output);
		}
	}
	return $output;
}
private function getNumberStumpMods($modpep) {
	$output = $modpep;
	$mass = $this->modmasses_long[0];
		$output = preg_replace("/K\[{$mass}\]/", "x", $output);
	if($this->deadend) {
		$mass = $this->modmasses_short[0];
		$output = preg_replace("/K\[{$mass}\]/", "x", $output);
	}
	return substr_count($output, "x");
}

public function readPepXML($samples = array(), $spec_filter = "", $deconvolute = false) {
	$file = $this->deadend ? $this->iqpir_params['deadendprophetfile'] : $this->iqpir_params['xlinkprophetfile'];
	
	$fdr = $this->iqpir_params['fdr'];
	$level = $this->iqpir_params['filter_crit'];
	$output_dir = $this->iqpir_params['iqpir_output_dir'];
	$prot_suff = $this->deadend ? "" : "s";
	$res_suff = $this->deadend ? "" : "_pair";

	$proton_mass = 1.00727647;
	$pep_massdiff = 2.0067095;
	$pepRatio = "L/H_";
	$reporterRatio = "H/L_";
	if($deconvolute) {
		$pepRatio = "dL/H_";
		$reporterRatio = "dH/L_";
	}
	$filter_suff = $spec_filter == "" ? "" : "_" . $spec_filter;
	$filter_suff = preg_replace("/\,/", "_", $filter_suff);

	if($spec_filter != "") {
		echo "Filtering pepXML {$file} for spectra containing the specified text {$spec_filter}\n";
		$spec_filter = split(",", $spec_filter);
	}
	else $spec_filter = array();
	
	$output_index = 8;
	$file_suffix = "-iqpir.xls2"; 
	$find_residuals = false; // wheterh to look for peaks with expected h/l spacing
	$pursue_all_pep_charges = true;
	
	$report_log2ratios = true;
	$protein_seqs_and_genes = array();

	$ratio_prefix = $report_log2ratios ? "log2" : "";
foreach($this->iqpir_params['samples'] as $sample => $bioreps) {	
	
	if(count($samples) > 0 && ! array_key_exists($sample, $samples)) continue;
	$this->setSampleChannels($sample);
	$stump_modmasses = join(",", array_map(function($item) { return $item[0]; },
    		array_values($this->offsets)));
    $reporter_masses = join(",", array_map(function($item) { return $item[1]; },
    		array_values($this->offsets)));


	$chromfile = $this->iqpir_params['iqpir_output_dir'] . "sample." . $sample . ".iqpir.txt";
	$numidsfile = $this->iqpir_params['iqpir_output_dir'] . "sample." . $sample . ".numids.txt";
	$fchrom = fopen($chromfile, "w");

	if($this->deadend) fwrite($fchrom, "deadend-peptide");
	else fwrite($fchrom, "cross-link");
	fwrite($fchrom, "\txlinkdb-id\tprotein{$prot_suff}\tuniprot\tresidue{$res_suff}\t");
	fwrite($fchrom, "rawfile\tscan\tid\tion_type\tmzs\ttotal_intensity\ttheor_relints\tobs_relints\tbiorep\tnoise\t");
	fwrite($fchrom, "channel_apportionment\tapportionment_error\tchannel_stumpmodmasses\treporter_modmasses"); //.join("\t", $pepfraginfo_keys)); //."\tratio\tratio_error");
	fwrite($fchrom, "\n");

	$existing_iqpirquant = array();
	$existing_quant_parameter = "existing_" . ($this->deadend ? "deadend" : "xl"). "_iqpirfiles";
	if(array_key_exists($existing_quant_parameter, $this->iqpir_params) && $this->iqpir_params[$existing_quant_parameter] != "") {
		$existing_iqpirfile_dirs = split(",", $this->iqpir_params[$existing_quant_parameter]); //"/net/gs/vol4/shared/brucelab/search/xiaoting/Spectrast_test/interact_spectrast-xl-iqpir.xls", "/net/gs/vol4/shared/brucelab/search/xiaoting/iqPIR_mito_Tac_sham_mango_6pairs_135runs/iprophet-110920-xl-iqpir.xls");
		$existing_iqpirfiles = array();
		for($k = 0; $k < count($existing_iqpirfile_dirs); $k++) {
			if($existing_iqpirfile_dirs[$k][strlen($existing_iqpirfile_dirs[$k])-1] != "/") $existing_iqpirfile_dirs[$k] .= "/";
			$next_existing_files = glob($existing_iqpirfile_dirs[$k] . "sample.*.iqpir.txt");
			for($j = 0; $j < count($next_existing_files); $j++) array_push($existing_iqpirfiles, $next_existing_files[$j]);
		}
		if(count($existing_iqpirfiles) == 0) {
			echo "Error: no existing iqpirfiles found in directories {$this->iqpir_params[$existing_quant_parameter]}\n";
			exit(1);
		}
		$existing_iqpirquant = $this->reuseIqpirfiles($existing_iqpirfiles, $sample);
		echo "Read in existing quantitation for ".count($existing_iqpirquant)." spectral ".
			($this->deadend ? "deadend" : "cross-link")." scans\n";
	}
	if($this->deadend) {
		if(file_exists($this->iqpir_params['xlinkprophetfile'])) {
			// get the modmass_info
			$valid = array();
			exec("grep xlinker_stump_modmasses {$this->iqpir_params['xlinkprophetfile']}", $valid);
			if(count($valid) > 0) {
				$valid[0] = rtrim($valid[0]);
				$nextpos = strpos($valid[0], "xlinker_stump_modmasses=\"");
				$modmass_info = substr($valid[0], $nextpos + strlen("xlinker_stump_modmasses=\""));
				$nextpos = strpos($modmass_info, "\"");
				$modmass_info = substr($modmass_info, 0, $nextpos);
				echo "Using modmassinfo {$modmass_info}\n";
			}
		}
		if($modmass_info == "") $modmass_info = "n:198.040247,K:325.127385,K:325.127385,K:327.134095";
		$pepxml = new MyDeadendPepXML($file, $fdr, array("modmass_info" => $modmass_info));
	}
	else $pepxml = new MyPepXML($file, $fdr, $level);
	$scan = 0;
	$peps = array();
	$masses = array(); 
	$charges = array();
	$mzxml = "";
	$noisefile = "";
	$analysis_type = "";
	$seen = array();
	$spectrum = "";
	$hk2_info = array(); // by mzXML
	$bar_graphs_peptide = array();
	$bar_graphs_reporter = array();
	$bar_graphs_combined = array();	
	$ms3scans = array();
	$rawfile = "";
	$xl_prots = array();
	$xl_uniprots = array();
	$xl_residuepairs = array();
	$xlinkdb_ids = array();
	$xl_looplinks = array(); // store the first 2 columns of iqpir
	$analyze = false;
	$is_sample = false;
	$peptide_frag_exclusions = array();
	$ms3_massoffsetptr = NULL;
	$massoffset_suff = ""; 
	$massoffset_header = true;
	while($pepxml->nextLine()) {
		if($pepxml->hasStartTag("msms_run_summary")) { // end of searh result
			$mzxml = $pepxml->getTagValue("base_name") . ".mzXML";
			$noisefile = $pepxml->getTagValue("base_name") . ".noise";
			$rawfile = $pepxml->getTagValue("base_name");
			$nextpos = strrpos($rawfile, "/");
			if($nextpos !== false) {
				$rawfile = substr($rawfile, $nextpos+1);
			}
			$is_sample = array_key_exists($rawfile, $this->iqpir_params['rawfiles']) && 
				$this->iqpir_params['rawfiles'][$rawfile]['sample'] == $sample;
			if(! file_exists($mzxml)) {
				echo "Error: mzxml {$mzxml} does not exist\n";
				exit(1);
			}
			if(! file_exists($noisefile)) {
				echo "Error: noisefile {$noisefile} does not exist\n";
				exit(1);
			}
		}
		else if(! $is_sample) {
			continue;
		}
		else if($pepxml->hasStartTag("search_summary")) { // end of searh result
			$analysis_type = $pepxml->getAnalysisType();
		}
		else if($pepxml->hasStartTag("spectrum_query")) { // end of searh result
			$spectrum = $pepxml->getTagValue("spectrum");
			if($analysis_type != "ReACT") { //"Mango" || $analysis_type == "Comet" || $analysis_type == "SpectraST" || $analysis) {
				$scan = $pepxml->getTagValue("start_scan");
			}
			if($this->deadend) {
				$analyze = true;
				array_push($charges, $pepxml->getTagValue("assumed_charge"));
			}
		}
		else if($analysis_type == "ReACT" && $pepxml->hasStartTag("search_hit")) { // end of searh result
			$scan = $pepxml->getTagValue("ms2scan");
		}
		else if($analyze && $pepxml->hasStartTag("search_hit")) { 
			array_push($masses, $pepxml->getTagValue("calc_neutral_pep_mass"));
			array_push($masses, $pepxml->getTagValue("calc_neutral_pep_mass"));
		}
		else if($analyze && $pepxml->hasEndTag("search_hit")) { // end of searh result
			$analyze = false;
		}
		else if($pepxml->hasStartTag("linked_peptide")) { // end of searh result
			array_push($masses, $pepxml->getTagValue("calc_neutral_pep_mass"));
			array_push($charges, $pepxml->getTagValue("assumed_charge"));
			if($analysis_type == "ReACT") { // end of searh result
				array_push($ms3scans, $pepxml->getTagValue("ms3_scan"));
			}
			if(count($masses) > 2) {
				echo "Error have ".count($masses)." masses at {$spectrum}\n";
				exit(1);
			}
		}
		else if(count($protein_seqs_and_genes)==0 && $pepxml->hasStartTag("search_database")) { // end of searh result
			$protein_seqs_and_genes = $this->getDatabaseProteinSequencesAndGenes($pepxml->getTagValue("local_path"));
		}
		else if((! $this->deadend || $analyze) && $pepxml->hasStartTag("modification_info")) { // end of searh result
			array_push($peps, $pepxml->getTagValue("modified_peptide"));
		}
		else if($pepxml->hasEndTag("spectrum_query")) { // end of searh result
			if(count($spec_filter) > 0) { // msut see if spectrum matches any
				$found = false;
				for($k = 0; $k < count($spec_filter); $k++) {
		 			if(strpos($spectrum, $spec_filter[$k])!==false) {
		 				$found = true;
		 				$k = count($spec_filter);
		 			}
		 		}
				if(! $found) continue;
			}
			$spec_without_ch = substr($spectrum, 0, strlen($spectrum)-2); // don't need to analyze different charges of same cross-link 
			$next_unique = ($analysis_type == "ReACT" ? $rawfile . "." . $scan . "." . $scan : $spec_without_ch) . "\t" . 
				$this->convertHeavyToLight($pepxml->getPreservedOrderedModifiedPeptidePair(false, false, $analysis_type == "Looplink"));
			$channel_offset = $this->getChannelOffset($peps[0]);
			$heavy = $channel_offset > 0; //$this->isHeavy($peps[0]); //strpos($peps[0], $heavy_mod)!==false;

			if($pepxml->aboveFilterThreshold() && array_key_exists($next_unique, $seen)) {
				$next_index = $heavy ? 1 : 0;
				$seen[$next_unique][$channel_offset]++; 
			} 
			else if($pepxml->aboveFilterThreshold($this->minprob_with_mincomposite) && ! array_key_exists($next_unique, $seen)) {
				if($channel_offset < 0) {
					//echo "Channel offset -1\n"; exit(1);
					$peps = array();
					$masses = array(); 
					$charges = array();
					$next_scans = array();
					$next_noise = array();
					$ms3scans = array();
					$longarm_scans = array();
					continue; // this modified peptide has no stump mod used in this sample
				}
				$seen[$next_unique] = array(-1 => 0);
				foreach($this->offsets as $offset => $val) $seen[$next_unique][$offset] = 0;
				$noise = 0;
				$noise_valid = array();
				exec("grep -P '^".$scan."\t' {$noisefile}", $noise_valid);
				if(count($noise_valid)===1) {
					$next_noise = split("\t", $noise_valid[0]);
					$noise = max(1, $next_noise[3]);
				}
				$next_index = $heavy ? 1 : 0;
				$seen[$next_unique][$channel_offset]++; 
				if(count($existing_iqpirquant) > 0 && array_key_exists($next_unique, $existing_iqpirquant)) {
					fwrite($fchrom, join("\n", array_values($existing_iqpirquant[$next_unique])) . "\n");
					$peps = array();
					$masses = array(); 
					$charges = array();
					$next_scans = array();
					$next_noise = array();
					$ms3scans = array();
					$longarm_scans = array();
					continue;
				}
				$num_isotopes = 5;
				$dlight1 = NULL;
				$proton_mass = 1.00727647;
				$min_mz = 500;
				$max_mz = 2000;
				
							$modpeps2 = $pepxml->getModifiedPeptidePair();
							if($heavy) { // must change sequence to light version
								for($m = 0; $m < count($modpeps2); $m++) {
									$modpeps2[$m] = $this->convertHeavyToLight($modpeps2[$m]); //preg_replace("/K\[327.13\]/", "K[325.13]", $modpeps[0]);
								}
							}
							if($analysis_type == "Looplink") {
								$next_info = split("\t", $next_unique);
								$modpeps2 = array($next_info[1], $next_info[1]);
							}
							
			$peptide_pair = $pepxml->getPreservedOrderedModifiedPeptidePair(false, true);
			$next_xl = $this->deadend ? $this->convertHeavyToLight($peptide_pair[0]) : 
				$this->convertHeavyToLight($peptide_pair[0]) . "_" .  $this->convertHeavyToLight($peptide_pair[3]);
if(! array_key_exists($next_xl, $xlinkdb_ids)) {
			$xl_peps = split("_", $next_xl);

			if(! $this->deadend && array_key_exists($xl_peps[1] . "_" . $xl_peps[0], $xlinkdb_ids)) { // already have it in reverse...
				$peptide_pair = array($peptide_pair[3], $peptide_pair[4], $peptide_pair[5],
					$peptide_pair[0], $peptide_pair[1], $peptide_pair[2]);
				$next_xl = $xl_peps[1] . "_" . $xl_peps[0];		
				echo "Switching next_xl {$xl_peps[0]}_{$xl_peps[1]} for {$xl_peps[1]}_{$xl_peps[0]} since previously seen in that order\n";
			
			}
			else {
				for($z = 0; $z < count($xl_peps); $z++) {
					$xl_peps[$z] = join("_", $this->getModifiedPeptideWithoutCrosslinkmodAndKpos($xl_peps[$z]));
				}
				$xlinkdb_ids[$next_xl] = join("_", $xl_peps);
				
				
				$xl_prots[$next_xl] = $this->deadend ? $protein_seqs_and_genes[$peptide_pair[1]]['gene'] : 
					$protein_seqs_and_genes[$peptide_pair[1]]['gene'] . "-" . $protein_seqs_and_genes[$peptide_pair[4]]['gene'];
				$xl_uniprots[$next_xl] = $this->deadend ? $peptide_pair[1] : $peptide_pair[1] . "-" . $peptide_pair[4];				
				$next_pepA = preg_replace("/L/", "I", $pepxml->stripMods($peptide_pair[0]));
				if(count($protein_seqs_and_genes) > 0 && ! array_key_exists($peptide_pair[1], $protein_seqs_and_genes)) {
					echo "Error: protein A {$peptide_pair[1]} not found in protein seqs: ".join(",", array_keys($protein_seqs_and_genes))."\n";
					exit(1);
				}
				$nextposA = strpos($protein_seqs_and_genes[$peptide_pair[1]]['sequence'], $next_pepA);
				if($nextposA===false) {
					echo "ErrorA: {$next_pepA} not found in {$peptide_pair[1]} sequence {$protein_seqs_and_genes[$peptide_pair[1]]['sequence']}\n";
					exit(1);
				}
				$nextResA = $nextposA + $peptide_pair[2] + 1; //getCrosslinkProteinRes($next_pepA, $peptide_pair[2], $protein_seqs_and_genes[$peptide_pair[1]]);
				$next_pepB = preg_replace("/L/", "I", $pepxml->stripMods($peptide_pair[3]));
				if(count($protein_seqs_and_genes) > 0 && ! array_key_exists($peptide_pair[4], $protein_seqs_and_genes)) {
					echo "Error: protein B {$peptide_pair[4]} not found in protein seqs: ".join(",", array_keys($protein_seqs_and_genes))."\n";
					exit(1);
				}					
				$nextposB = strpos($protein_seqs_and_genes[$peptide_pair[4]]['sequence'], $next_pepB);
				if($nextposB===false) {
					echo "ErrorB: {$next_pepB} not found in {$peptide_pair[4]} sequence {$protein_seqs_and_genes[$peptide_pair[4]]['sequence']}\n";
					exit(1);
				}
				$nextResB = $nextposB + $peptide_pair[5] + 1; //getCrosslinkProteinRes($next_pepB, $peptide_pair[5], $protein_seqs_and_genes[$peptide_pair[4]]);
				$xl_residuepairs[$next_xl] = $this->deadend ? $nextResA : $nextResA . "-" . $nextResB;							

				if($analysis_type == "Looplink") {
					$xl_looplinks[$next_xl] = $modpeps2[0] . "_LOOPLINK\t"; 
					$next_p = split("_", $xlinkdb_ids[$next_xl]);
					$xl_looplinks[$next_xl] .= preg_replace("/K\[128.09\]/", "K", $next_p[0]) . "_" . min($next_p[1], $next_p[3]) .
						"_LOOPLINK_" . max($next_p[1], $next_p[3]);
					$next_r = split("-", $xl_residuepairs[$next_xl]);
					sort($next_r);
					$xl_residuepairs[$next_xl] = join("-", $next_r);
				}




		} // if not seen forwards

	if(true) {						
				$next_scans = ! $this->deadend && $analysis_type == "ReACT" ? $ms3scans : array($scan, $scan);
				$next_noise = ! $this->deadend && $analysis_type == "ReACT" ? array(0, 0) : array($noise, $noise);;
				if($massoffset_suff != "" && ! $this->deadend) { // && $analysis_type == "ReACT") {
					$massoffsetfile = $this->iqpir_params['iqpir_output_dir'] . $massoffset_suff;
					if(is_null($ms3_massoffsetptr)) {
						if(file_exists($massoffsetfile)) {
							$ms3_massoffsetptr = fopen($massoffsetfile, "a");
							$massoffset_header = false;
						}
						else {
							$ms3_massoffsetptr = fopen($massoffsetfile, "w");
						}
					}
					$this->recordMs3Offsets($modpeps2[0], $next_scans[0], $mzxml, $massoffset_header, $ms3_massoffsetptr);
					$this->recordMs3Offsets($modpeps2[1], $next_scans[1], $mzxml, $massoffset_header, $ms3_massoffsetptr);
				}

				$pep1_stumpfrags = $this->getStumpFragments($modpeps2[0]); //325.13); //exit(1);
				$pep2_stumpfrags = $this->deadend || ($modpeps2[0] == $modpeps2[1] && $next_scans[0] == $next_scans[1]) ? array() : $this->getStumpFragments($modpeps2[1]); //325.13); //exit(1);
				$frag_ratios = array();
				$frag_ints = array();
				$seen_fragmasses = array(); // make sure only use once
				$fragment_relints = array();
				$pep1_frag_ratios = array();
				$pep1_frag_ints = array();
				$pep2_frag_ratios = array();
				$pep2_frag_ints = array();
				$frag_fragments = array();
				$frag_mzs = array();
				$this->ms2_scan = $scan;
				$longarm_scans = array_key_exists('longarm_reporters', $this->iqpir_params) && 
					$this->iqpir_params['longarm_reporters'] == "true" && $analysis_type == "ReACT" ?
					array($ms3scans[0] + 2, $ms3scans[1] + 2) : array();
					
				$this->setScanPeaklist($mzxml, $ms3scans, $longarm_scans);
				$fragment_theor_ints = array();
				$peptide_theor_ints = array();
				$peptide_relints = array();
				$peptide_apports = array();
				$peptide_errors = array();
				$frag_apports = array();
				$frag_errors = array();
				$reporter_apport = "";
				$this->num_stumpmods = 1;				
				$samepeps = $this->deadend || $modpeps2[0]==$modpeps2[1] || $peptide_pair[0]==$peptide_pair[3];
				$first_2cols = array_key_exists($next_xl, $xl_looplinks) ? $xl_looplinks[$next_xl] : 
					$next_xl . "\t" . $xlinkdb_ids[$next_xl];
				
				$min_peptide_massdiff = 0.5; // daltons?
				$indistinct_pepmasses = abs($masses[0]-$masses[1]) < $min_peptide_massdiff;
				$all_mzs = array($scan => array(), $next_scans[0] => array(), $next_scans[1] => array()); //for this scan
				$excluded_ions = array($scan => array(), $next_scans[0] => array(), $next_scans[1] => array()); 
				// enter peptide and reporters ????
				$pep_massdiff = $channel_offset == 0 ? 0 : $this->offsets[$channel_offset][0] - $this->offsets[0][0];
				for($next_charge = 1;  $next_charge < 3; $next_charge++) {
					foreach($this->offsets as $offset => $info) {
						$next_mz = ($masses[0] - $this->offsets[$channel_offset][0] + $this->offsets[0][0] + $proton_mass)/$next_charge;
						$ppm_diff = $next_mz * $this->iqpir_ppmtolerance / 1000000;
						$all_mzs[$scan]["PepA" . $next_charge][$offset] = array($next_mz - $ppm_diff, $next_mz + $ppm_diff);
						if($masses[1] != $masses[0]) {
							$next_mz = ($masses[1] - $this->offsets[$channel_offset][0] + $this->offsets[0][0] + $proton_mass)/$next_charge;
							$ppm_diff = $next_mz * $this->iqpir_ppmtolerance / 1000000;
							$all_mzs[$scan]["PepB" . $next_charge][$offset] = array($next_mz - $ppm_diff, $next_mz + $ppm_diff);
						}
						$next_mz = $this->offsets[$offset][1] + $proton_mass;
						$ppm_diff = $next_mz * $this->iqpir_ppmtolerance / 1000000;
						$all_mzs[$scan]["Rep"][$offset] = array($next_mz - $ppm_diff, $next_mz + $ppm_diff);
					}
				}
				foreach($pep1_stumpfrags as $key => $value) {	

					if($analysis_type == "Looplink") $this->num_stumpmods = $this->getNumberStumpMods($value[1]);
					$next_mass = strval($value[0]);
					$nextpos = strpos($next_mass, ".");
					if($nextpos!==false) $next_mass = substr($next_mass, 0, $nextpos+4); // keep only first 2 digits
					if(array_key_exists($next_mass, $seen_fragmasses)) continue;
					$seen_fragmasses[$next_mass] = 1;
					
					// need function that returns both the isotope relative amounts as well as the total intensity
					$next_quant = $this->getIsotopeQuant_withOffsets($value[0], 1, $next_scans[0], $mzxml, $noise, "frag_{$key}_1");
					$next_isotope_pk_relintens = $next_quant[0];
					$next_fragintens = $next_quant[1];
					// now compute the ratios
					$next_ints_fragentries = $this->getIsotopeIntensities($value[1], true);
					if($next_ints_fragentries[0]==0) {
						echo "1. No theoreteical intensities obtained with calcisotopes for fragment ion \"{$value[1]}\": ".join(",", $next_ints_fragentries).", possibly due to unexpected peptide modification (contact Jimmy).  Will instead pass mass {$value[0]}\n";
						$next_ints_fragentries = $this->getIsotopeIntensities($value[0], true);
					}
					$next_multi = $next_fragintens == 0 ? array() : $this->multiplex_withOffsets($next_ints_fragentries, $next_isotope_pk_relintens);
					$max_error = 0.05;
					if(count($next_multi) > 0) { // && $next_multi[0][0] > 0 && $next_multi[0][2] > 0 && $next_multi[1] <= 0.25) {
						$next_id = ! $this->deadend && $analysis_type != "Looplink" &&
							$samepeps && $next_scans[0] == $next_scans[1] ? "AB:" . $key . "_{$this->num_stumpmods}"  : 
							"A:" . $key . "_1";
						$next_mzs = array();
						$zero_stump = $this->offsets[0][0];
						
						foreach($this->offsets as $offset => $info) {
							array_push($next_mzs, number_format($value[0] + $proton_mass + 
								$this->num_stumpmods * ($info[0] - $zero_stump), 7, ".", ""));
						}
						$next_mzs = join(",", $next_mzs); 
						$next_theor_relints = join(",", $next_ints_fragentries);
						$next_obs_relints = join(",", $next_isotope_pk_relintens);
						$next_channel_apport = join(",", $next_multi[0]);
						$next_apport_err = $next_multi[1];
	fwrite($fchrom, 
	
	$first_2cols . "\t" . $xl_prots[$next_xl] . "\t" . $xl_uniprots[$next_xl]. "\t". $xl_residuepairs[$next_xl]);
	fwrite($fchrom, "\t{$rawfile}\t{$next_scans[0]}\t{$next_id}\tFragment\t{$next_mzs}\t{$next_fragintens}\t{$next_theor_relints}");
	fwrite($fchrom, "\t{$next_obs_relints}\t".$this->iqpir_params['rawfiles'][$rawfile]['biorep']."\t{$next_noise[0]}\t{$next_channel_apport}\t{$next_apport_err}\t"); 
	fwrite($fchrom, "{$stump_modmasses}\t{$reporter_masses}\n");  //.join("\t", $pepfraginfo_keys)); //."\tratio\tratio_error");
					}
				} 

				foreach($pep2_stumpfrags as $key => $value) {
					$next_mass = strval($value[0]);
					$nextpos = strpos($next_mass, ".");
					if($nextpos!==false) $next_mass = substr($next_mass, 0, $nextpos+4);
					if(array_key_exists($next_mass, $seen_fragmasses)) continue;
					$seen_fragmasses[$next_mass] = 1;

					$next_quant = $this->getIsotopeQuant_withOffsets($value[0], 1, $next_scans[1], $mzxml, $noise, "frag_{$key}_1");
					$next_isotope_pk_relintens = $next_quant[0];
					$next_fragintens = $next_quant[1];
					// now compute the ratios
					$next_ints_fragentries = $this->getIsotopeIntensities($value[1], true);
					if($next_ints_fragentries[0]==0)  {
						echo "2. No theoreteical intensities obtained with calcisotopes for channel offset {channel_offset} fragment ion \"{$value[1]}\": ".join(",", $next_ints_fragentries).", possibly due to unexpected peptide modification (contact Jimmy).  Will instead pass mass {$value[0]}\n";
						$next_ints_fragentries = $this->getIsotopeIntensities($value[0], true);
					}
					$next_multi = $next_fragintens == 0 ? array() : $this->multiplex_withOffsets($next_ints_fragentries, $next_isotope_pk_relintens);
						$next_id = "B:" . $key . "_1";
						$next_mzs = array();
						$zero_stump = $this->offsets[0][0];
						foreach($this->offsets as $offset => $info) {
							array_push($next_mzs, number_format($value[0] + $proton_mass + $info[0] - $zero_stump, 7, ".", ""));
						}
						$next_mzs = join(",", $next_mzs); 
						$next_theor_relints = join(",", $next_ints_fragentries);
						$next_obs_relints = join(",", $next_isotope_pk_relintens);
						$next_channel_apport = join(",", $next_multi[0]);
						$next_apport_err = $next_multi[1];
	fwrite($fchrom, 
	
	$first_2cols . "\t" . $xl_prots[$next_xl] . "\t" . $xl_uniprots[$next_xl]. "\t". $xl_residuepairs[$next_xl]);
	fwrite($fchrom, "\t{$rawfile}\t{$next_scans[1]}\t{$next_id}\tFragment\t{$next_mzs}\t{$next_fragintens}\t{$next_theor_relints}");
	fwrite($fchrom, "\t{$next_obs_relints}\t".$this->iqpir_params['rawfiles'][$rawfile]['biorep']."\t{$next_noise[1]}\t{$next_channel_apport}\t{$next_apport_err}\t"); 
	fwrite($fchrom, "{$stump_modmasses}\t{$reporter_masses}\n");  //.join("\t", $pepfraginfo_keys)); //."\tratio\tratio_error");
					}


				} 
				if($remove_interference && count($pep2_stumpfrags) > 0 && 
					! array_key_exists($modpeps2[1], $peptide_frag_exclusions)) $peptide_frag_exclusions[$modpeps2[1]] = $excluded_ions[$next_scans[1]];
	} // if false, no fragments			


			
				
}				
				$this->num_stumpmods = $analysis_type == "Looplink" ? 2 : 1;

				$light_entries1 = $pursue_all_pep_charges ? array("1" => array(), "2" => array(), "3" => array()) : array(); // charges 1,2,3

				if($pursue_all_pep_charges && $this->num_stumpmods == 2) $light_entries1["4"] = array(); 
				
				$heavy_entries1 = array("1" => array(), "2" => array(), "3" => array()); // charges 1,2,3
				$mzs1 = array();  // from charge to array of light/heavy
				// should also compute the ratio of the additional masses at charge 1 and 2.....
				$pep_massdiff = $channel_offset == 0 ? 0 : $this->num_stumpmods * ($this->offsets[$channel_offset][0] - $this->offsets[0][0]);
				$multi_entries1 = array();
				$mutli_entries2 = array();
				foreach($light_entries1 as $key => $value) {
					$next_charge = intval($key);
					$symbol = $next_charge == $charges[0] ? "*" : "";
					$mzs1[$key] = 0;
					if($heavy) {
						if($masses[0]/$next_charge >= $min_mz && $masses[0]/$next_charge <= $max_mz) {
							$mzs1[$key] = array(($masses[0] - $pep_massdiff + $next_charge * $proton_mass)/$next_charge, ($masses[0] + $next_charge * $proton_mass)/$next_charge);
						}
					}
					else {
						if(($masses[0]+$pep_massdiff)/$next_charge >= $min_mz && ($masses[0]+$pep_massdiff)/$next_charge <= $max_mz) {
							$mzs1[$key] = array(($masses[0] + $next_charge * $proton_mass)/$next_charge, ($masses[0] + $pep_massdiff + $next_charge * $proton_mass)/$next_charge);
						}
					}
					if($mzs1[$key] == 0) continue;
					$next_quant = $this->getIsotopeQuant_withOffsets($masses[0] - $pep_massdiff, $next_charge, $scan, $mzxml, $noise, "{$symbol}pep{$key}_1 ");
					$next_isotope_pk_relintens = $next_quant[0];
					$next_pepintens = $next_quant[1];
					$multi_entries1[$key] = array("total_intensity" => $next_pepintens, "obs_intensities" => $next_isotope_pk_relintens);
					$next_ints_pepentry = $this->getIsotopeIntensities($modpeps2[0], true);
					$multi_entries1[$key]['theor_intensities'] = $next_ints_pepentry;
					if($next_ints_pepentry[0]==0)  {
						echo "2a. No theoreteical intensities obtained with calcisotopes for channel offset {$channel_offset} peptide ion {$peps[0]} \"{$modpeps2[0]}\": ".join(",", $next_ints_pepentry).", possibly due to unexpected peptide modification (contact Jimmy).  Will instead pass mass {$masses[0]}\n";
						$next_ints_pepentry = $this->getIsotopeIntensities($masses[0], true);
					}
					$multi_entries1[$key]['apportionment'] = $next_pepintens == 0 ? array() : $this->multiplex_withOffsets($next_ints_pepentry, $next_isotope_pk_relintens);
					if(count($multi_entries1[$key]['apportionment']) > 0) { // && $multi_entries1[$key]['apportionment'][0][0] > 0 && 
							$next_id = ! $this->deadend && $analysis_type != "Looplink" && $samepeps ? "AB:{$key}"  : "A:{$key}";
							$next_mzs = array();
							$zero_stump = $this->offsets[0][0];
							foreach($this->offsets as $offset => $info) {
								array_push($next_mzs, number_format(($mzs1[$key][0] * $key + 
								$this->num_stumpmods * ($info[0] - $zero_stump))/$key, 7, ".", ""));
							}
							$next_mzs = join(",", $next_mzs); 
							$next_theor_relints = join(",", $multi_entries1[$key]['theor_intensities']);
							$next_obs_relints = join(",", $multi_entries1[$key]['obs_intensities']);
							$next_channel_apport = join(",", $multi_entries1[$key]['apportionment'][0]);
							$next_apport_err = $multi_entries1[$key]['apportionment'][1];
	fwrite($fchrom, //$next_xl . "\t" . $xlinkdb_ids[$next_xl] . "\t" . 
	$first_2cols . "\t" . $xl_prots[$next_xl] . "\t" . $xl_uniprots[$next_xl]. "\t". $xl_residuepairs[$next_xl]);
	fwrite($fchrom, "\t{$rawfile}\t{$scan}\t{$next_id}\tPeptide\t{$next_mzs}\t{$next_pepintens}\t{$next_theor_relints}");
	fwrite($fchrom, "\t{$next_obs_relints}\t".$this->iqpir_params['rawfiles'][$rawfile]['biorep']."\t{$noise}\t{$next_channel_apport}\t{$next_apport_err}\t"); 
	fwrite($fchrom, "{$stump_modmasses}\t{$reporter_masses}\n");  //.join("\t", $pepfraginfo_keys)); //."\tratio\tratio_error");

						}


				} // end of peptide1
	if(! $this->deadend && ! $samepeps) {			
				if($heavy) {
					$light1ints = $this->getIsotopeQuant($masses[0] - $pep_massdiff, $charges[0], $scan, $mzxml, $noise, "light1 ");
					$heavy1ints = $this->getIsotopeQuant($masses[0], $charges[0], $scan, $mzxml, $noise, "heavy1 ");
					if($deconvolute) $dlight1 = $this->getIsotopeQuant($masses[0] - $pep_massdiff, $charges[0], $scan, $mzxml, $noise, "light1 ", $num_isotopes);
				}
				else {
					$light1ints = $this->getIsotopeQuant($masses[0], $charges[0], $scan, $mzxml, $noise, "light1 ");
					$heavy1ints = $this->getIsotopeQuant($masses[0] + $pep_massdiff , $charges[0], $scan, $mzxml, $noise, "heavy1 ");
					if($deconvolute) $dlight1 = $this->getIsotopeQuant($masses[0], $charges[0], $scan, $mzxml, $noise, "light1 ", $num_isotopes);
				}
				$ratios1 = $this->printRatios($light1ints[0], $heavy1ints[0], "Scan {$scan} Pep1 Ratios: ");
				
				$ints1 = array();


				$dlight2 = NULL;

				$light_entries2 = $pursue_all_pep_charges ? array("1" => array(), "2" => array(), "3" => array()) : array(); // charges 1,2,3
				$heavy_entries2 = array("1" => array(), "2" => array(), "3" => array()); // charges 1,2,3
				$dlight1_entries2 = array("1" => array(), "2" => array(), "3" => array()); // charges 1,2,3
				$ratios_entries2 = array("1" => array(), "2" => array(), "3" => array()); // charges 1,2,3
				$ints_entries2 = array("1" => array(), "2" => array(), "3" => array()); // charges 1,2,3
				$deconv_ratio_entries2 = array("1" => array(), "2" => array(), "3" => array()); // charges 1,2,3
				$mzs2 = array();  // from charge to array of light/heavy
				foreach($light_entries2 as $key => $value) {

					if($remove_interference && array_key_exists($modpeps2[1], $peptide_frag_exclusions) &&
						array_key_exists("PepB{$key}", $peptide_frag_exclusions[$modpeps2[0]])) {
							echo "Skipping PepB{$key} for {$modpeps2[1]} due to interference with fragment ion\n";
							continue;
					}
					$next_charge = intval($key);
					$symbol = $next_charge == $charges[1] ? "*" : "";
					$mzs2[$key] = 0;
					if($heavy) {
						if($masses[1]/$next_charge >= $min_mz && $masses[1]/$next_charge <= $max_mz) {
							$mzs2[$key] = array(($masses[1] - $pep_massdiff + $next_charge * $proton_mass)/$next_charge, ($masses[1] + $next_charge * $proton_mass)/$next_charge);
						}
					}
					else {
						if(($masses[1]+$pep_massdiff)/$next_charge >= $min_mz && ($masses[1]+$pep_massdiff)/$next_charge <= $max_mz) {
							$mzs2[$key] = array(($masses[1] + $next_charge * $proton_mass)/$next_charge, ($masses[1] + $pep_massdiff + $next_charge * $proton_mass)/$next_charge);
						}
					}
					if($mzs2[$key] == 0) continue;

					$next_quant = $this->getIsotopeQuant_withOffsets($masses[1] - $pep_massdiff, $next_charge, $scan, $mzxml, $noise, "{$symbol}pep{$key}_1 ");
					$next_isotope_pk_relintens = $next_quant[0];
					$next_pepintens = $next_quant[1];
					$multi_entries2[$key] = array("total_intensity" => $next_pepintens, "obs_intensities" => $next_isotope_pk_relintens);

					$next_ints_pepentry = $this->getIsotopeIntensities($modpeps2[1], true);
					$multi_entries2[$key]['theor_intensities'] = $next_ints_pepentry;
					if($next_ints_pepentry[0]==0)  {
						echo "2b. No theoreteical intensities obtained with calcisotopes for channel offset {$channel_offset} peptide ion {$peps[1]} \"{$modpeps2[1]}\": ".join(",", $next_ints_fragentries).", possibly due to unexpected peptide modification (contact Jimmy).  Will instead pass mass {$masses[1]}\n";
						$next_ints_pepentry = $this->getIsotopeIntensities($masses[1], true);
					}
					$multi_entries2[$key]['apportionment'] = $next_pepintens == 0 ? array() : $this->multiplex_withOffsets($next_ints_pepentry, $next_isotope_pk_relintens);
					if(count($multi_entries2[$key]['apportionment']) > 0) { // && $multi_entries2[$key]['apportionment'][0][0] > 0 && 
							$next_id = "B:{$key}";
							$next_mzs = array();
							$zero_stump = $this->offsets[0][0];
							foreach($this->offsets as $offset => $info) {
								array_push($next_mzs, number_format(($mzs2[$key][0] * $key + $info[0] - $zero_stump)/$key, 7, ".", ""));
							}
							$next_mzs = join(",", $next_mzs); 
							$next_theor_relints = join(",", $multi_entries2[$key]['theor_intensities']);
							$next_obs_relints = join(",", $multi_entries2[$key]['obs_intensities']);
							$next_channel_apport = join(",", $multi_entries2[$key]['apportionment'][0]);
							$next_apport_err = $multi_entries2[$key]['apportionment'][1];
	fwrite($fchrom, //$next_xl . "\t" . $xlinkdb_ids[$next_xl] . "\t" . 
	$first_2cols . "\t" . $xl_prots[$next_xl] . "\t" . $xl_uniprots[$next_xl]. "\t". $xl_residuepairs[$next_xl]);
	fwrite($fchrom, "\t{$rawfile}\t{$scan}\t{$next_id}\tPeptide\t{$next_mzs}\t{$next_pepintens}\t{$next_theor_relints}");
	fwrite($fchrom, "\t{$next_obs_relints}\t".$this->iqpir_params['rawfiles'][$rawfile]['biorep']."\t{$noise}\t{$next_channel_apport}\t{$next_apport_err}\t"); 
	fwrite($fchrom, "{$stump_modmasses}\t{$reporter_masses}\n");  //.join("\t", $pepfraginfo_keys)); //."\tratio\tratio_error");
	}

				} // done peptide 2
} // not deadend
$this->num_stumpmods = 1;


if($this->quantify_ms2_reporters || count($longarm_scans) > 0) {
$this->verbose = true;
$this->reporter_apportion = true;
	if($this->quantify_ms2_reporters && ( ! $remove_interference || ! array_key_exists($modpeps2[0], $peptide_frag_exclusions) ||
						array_key_exists("Rep", $peptide_frag_exclusions[$modpeps2[0]]))) {

				$next_quant = $this->getIsotopeQuant_withOffsets($this->offsets[0][1], 1, $scan, $mzxml, $noise, "R_1");
				$next_isotope_pk_relintens = $next_quant[0];
				$next_repintens = $next_quant[1];
		if(count($this->reporter_theor_ints) > 0) $next_ints_repentries = $this->reporter_theor_ints;
		else
					// now compute the ratios
				$next_ints_repentries = $this->getIsotopeIntensities($this->offsets[0][1], true, $this->num_reporter_pks);
					// use $this->offsets to dictate and output only the top
				$next_multi = $next_repintens == 0 ? array() : 
					$this->multiplex_withOffsets($next_ints_repentries, $next_isotope_pk_relintens);
				$max_error = 0.05;
				$proton_mass = 1.00727647;

				if(count($next_multi) > 0) { // && $next_multi[0][0] > 0 && $next_multi[0][2] > 0 && $next_multi[1] <= 0.25) {
					$next_id = "R:1";
					$next_mzs = array();
					$next_obs_relints =  array();
					foreach($this->offsets as $offset => $info) {
						array_push($next_mzs, $info[1] + $proton_mass);
						for($p = 0; $p < $this->num_reporter_pks; $p++) {
							array_push($next_obs_relints, $next_isotope_pk_relintens[-2 * $offset + $p]);
						}
					}
					$next_mzs = join(",", array_reverse($next_mzs)); 
					$next_theor_relints = join(",", $next_ints_repentries);
					$next_obs_relints = 
					join(",", $next_isotope_pk_relintens);
					$next_apports = array_reverse($next_multi[0]);
					krsort($next_multi[0]);
					$next_channel_apport = join(",", array_values($next_multi[0]));
					$next_apport_err = $next_multi[1];
	fwrite($fchrom, 
	$first_2cols . "\t" . $xl_prots[$next_xl] . "\t" . $xl_uniprots[$next_xl]. "\t". $xl_residuepairs[$next_xl]);
	fwrite($fchrom, "\t{$rawfile}\t{$scan}\t{$next_id}\tReporter\t{$next_mzs}\t{$next_repintens}\t{$next_theor_relints}");
	fwrite($fchrom, "\t{$next_obs_relints}\t".$this->iqpir_params['rawfiles'][$rawfile]['biorep']."\t{$noise}\t{$next_channel_apport}\t{$next_apport_err}\t"); 
	fwrite($fchrom, "{$stump_modmasses}\t{$reporter_masses}\n");  //.join("\t", $pepfraginfo_keys)); //."\tratio\tratio_error");
			}
	} // if not excluded
	else if($this->quantify_ms2_reporters) echo "Skipping Reporter for {$modpeps2[0]}-{$modpeps2[1]} due to interference with fragment ion\n";
	for($k = 0; $k < count($longarm_scans); $k++) {		
				$next_quant = $this->getIsotopeQuant_withOffsets($this->offsets[0][1], 1, $longarm_scans[$k], $mzxml, $noise, "A:L{$k}_1");
				$next_isotope_pk_relintens = $next_quant[0];
				$next_repintens = $next_quant[1];
		if(count($this->reporter_theor_ints) > 0) $next_ints_repentries = $this->reporter_theor_ints;
		else {
					// now compute the ratios
				$next_ints_repentries = $this->getIsotopeIntensities($this->offsets[0][1], true, $this->num_reporter_pks);
					// use $this->offsets to dictate and output only the top
				$next_multi = $next_repintens == 0 ? array() : 
					$this->multiplex_withOffsets($next_ints_repentries, $next_isotope_pk_relintens);
				$max_error = 0.05;
				$proton_mass = 1.00727647;

				if(count($next_multi) > 0) { // && $next_multi[0][0] > 0 && $next_multi[0][2] > 0 && $next_multi[1] <= 0.25) {
					$next_id = ($k == 0 ? "A" : "B") . ":LRep";
					$next_mzs = array();
					$next_obs_relints =  array();
					foreach($this->offsets as $offset => $info) {
						array_push($next_mzs, $info[1] + $proton_mass);
						for($p = 0; $p < $this->num_reporter_pks; $p++) {
							array_push($next_obs_relints, $next_isotope_pk_relintens[-2 * $offset + $p]);
						}
					}
					$next_mzs = join(",", array_reverse($next_mzs)); 
					$next_theor_relints = join(",", $next_ints_repentries);
					$next_obs_relints = 
					join(",", $next_isotope_pk_relintens);
					$next_apports = array_reverse($next_multi[0]);
					krsort($next_multi[0]);
					$next_channel_apport = join(",", array_values($next_multi[0]));
					$next_apport_err = $next_multi[1];
	fwrite($fchrom, //$next_xl . "\t" . $xlinkdb_ids[$next_xl] . "\t" . 
	$first_2cols . "\t"  . $xl_prots[$next_xl] . "\t" . $xl_uniprots[$next_xl]. "\t". $xl_residuepairs[$next_xl]);
	fwrite($fchrom, "\t{$rawfile}\t{$longarm_scans[$k]}\t{$next_id}\tLongarmReporter\t{$next_mzs}\t{$next_repintens}\t{$next_theor_relints}");
	fwrite($fchrom, "\t{$next_obs_relints}\t".$this->iqpir_params['rawfiles'][$rawfile]['biorep']."\t{$next_noise[$k]}\t{$next_channel_apport}\t{$next_apport_err}\t"); 
	fwrite($fchrom, "{$stump_modmasses}\t{$reporter_masses}\n");  //.join("\t", $pepfraginfo_keys)); //."\tratio\tratio_error");
			}
			
			
			
	}
	
$this->verbose = false;
$this->reporter_apportion = false;

}
			} // above thresh
			$peps = array();
			$masses = array(); 
			$charges = array();
			$next_scans = array();
			$next_noise = array();
			$ms3scans = array();
			$longarm_scans = array();
			
		} // end of spec query

	} // next sample
	//fclose($fout);
	echo "Results written to {$chromfile}\n";
	echo "\n";
	
	$fout_lh = fopen($numidsfile, "w");
	$stump_mods = array_map(function($item) { return $item[0] . "_id"; },
    		array_values($this->offsets));
	fwrite($fout_lh, "spectrum\tpepA\tproA\tkposA\tpepB\tproB\tkposB\t".join("\t", $stump_mods)."\n");
	foreach($seen as $xl => $light_heavy) {
		fwrite($fout_lh, $xl . "\t" . join("\t", array_values($light_heavy)) . "\n");
	}
	fclose($fout_lh);
	echo "Channel identification info written to {$numidsfile}\n";
	echo "\n";
} // next sample	
}

private function getWeightedAverageRatio($ratios, $ints = NULL, $outlier_remove = false, $norm_offset = 0, 
	$verbose = false, $pvalue_with_tstat = false, $no_output_wt = false, $num_dec = 2, $no_output_pval = false) {
	$default_pval = "N/A";
	if($pvalue_with_tstat) $default_pval = "\t\tN/A";
	if($outlier_remove) {
		$vals = array();
		for($k = 0;  $k < count($ratios); $k++) {
			if(is_null($ints)) {
				array_push($vals, array($ratios[$k] + $norm_offset, 1));
			}
			else {
				array_push($vals, array($ratios[$k] + $norm_offset, sqrt($ints[$k])));
			}
			if($verbose) echo "Substituting ".($ratios[$k] + $norm_offset)." for {$ratios[$k]}\n";
		}
		if($verbose) {
			$printvals = array();
			for($z = 0; $z < count($vals); $z++) {
				array_push($printvals, $vals[$z][0]);
			}
			echo "Have input values ".join(",", $printvals)." with norm offset {$norm_offset}\n";
		}
		$rejects = array();
		$total = 0;
		$total_wt = 0;
		MyStatistics::removeOutliersByIndex($vals, $rejects, 0, $verbose); // index 1
		if($verbose) {
			$printvals = array();
			for($z = 0; $z < count($vals); $z++) {
				array_push($printvals, $vals[$z][0]);
			}
			 echo "Have outliers removed ".join(",", $rejects)." and remaining vals ".join(",", $printvals)."\n";
		}
		// now get weighted average of the remainder
		$df = count($vals) - 1 ;
		// can modify here to use sqrt of intensity
		for($k = 0; $k < count($vals); $k++) {
			$total += $vals[$k][0] * $vals[$k][1];
			$total_wt += $vals[$k][1];
		}
		if($total_wt > 0 && count($vals) > 1) {
			// compute the weighted standard deviation
			$std = 0;
			$mean = $total / $total_wt;
			for($k = 0; $k < count($vals); $k++) {
				$std += $vals[$k][1] * ($vals[$k][0] - $mean) * ($vals[$k][0] - $mean);
			}
			$std *= count($vals) / ((count($vals) - 1) * $total_wt);
			$std = sqrt($std);
			$output = array(number_format($mean, $num_dec, ".", ""), count($vals));
			if(! $no_output_wt) array_push($output, $total_wt);
			array_push($output, number_format($std, $num_dec + 1, ".", ""));
			if(! $no_output_pval) array_push($output, $this->getPvalue($mean, $std, $df, $pvalue_with_tstat)); //, $pval_num));
			array_push($output, join(",", $rejects));
			return $output;
		}
		else if($total_wt > 0 && count($vals) === 1) {
			$output = array(number_format($vals[0][0], 2, ".", ""), 1);
			if(! $no_output_wt) array_push($output, $vals[0][1]);	
			array_push($output, 0);
			if(! $no_output_pval) array_push($output, $default_pval);
			array_push($output, "");
			return $output;
		}
		$output =  array("", 0);
		if(! $no_output_wt) array_push($output, 0);			
		array_push($output, 0);
		if(! $no_output_pval) array_push($output, $default_pval);
		array_push($output, "");
		return $output;
	}
	$total_wt = 0;
	$total = 0;
	for($k = 0;  $k < count($ratios); $k++) {
		if(is_null($ints)) {
			$total += $ratios[$k] + $norm_offset;
			$total_wt += 1;
		}
		else {
			$total += ($ratios[$k] + $norm_offset) * $ints[$k];
			$total_wt += $ints[$k];
		}
	}
	if($total_wt > 0 && count($ratios) > 1)  {
		$std = 0;
		$mean = $total / $total_wt;
		for($k = 0; $k < count($ratios); $k++) {
			if(is_null($ints)) {
				$std += ($ratios[$k] + $norm_offset - $mean) * ($ratios[$k] + $norm_offset - $mean);
			}
			else {
				$std += $ints[$k] * ($ratios[$k] + $norm_offset - $mean) * ($ratios[$k] + $norm_offset - $mean);
			}
		}
		$std *= count($ratios) / ((count($ratios)-1) * $total_wt);
		$std = sqrt($std);
		$output = array(number_format($mean, $num_dec, ".", ""), count($ratios), $total_wt, number_format($std, $num_dec + 1, ".", ""));
		if(! $no_output_pval) array_push($output, $this->getPvalue($mean, $std, count($ratios)-1, $pvalue_with_tstat));
		return $output;
	}
	if($total_wt > 0 && count($ratios) === 1)  {
		$output = array($ratios[0], 1, $ints[0], 0);
	}
	else $output = 	array("", 0, 0, 0);
		
	if(! $no_output_pval) array_push($output, $default_pval);
	return $output; //array("", 0, 0, 0, $default_pval);


}


public function getPvalue($mean, $std, $df, $pval_with_tstat = false) { //, $pvalue_num = "") {
	if($df < 1) return "N/A\t0\tN/A";
	$next_stats = array($mean, $std, $df + 1);
	if($next_stats[1]==-1) {
		echo "Here with ".join(",", $next_stats)."\n";
		exit(1);
	}
	$t_stat = MyStatistics::get_tstatistic($next_stats);
	$pval = MyStatistics::get_tpvalue($t_stat, $df);
	if(! $pval_with_tstat) return $pval;
	return $t_stat . "\t" . $df . "\t". $pval;
}

private function getModifiedPeptideWithoutCrosslinkmodAndKpos($modpep) {
	$copy = $modpep;
	$modmasses = array(array("K", "327.13"), array("K", "327.14"), array("K", "325.13"));
	for($k = 0; $k < count($modmasses); $k++) {
		$modpep = preg_replace( "/{$modmasses[$k][0]}\[{$modmasses[$k][1]}\]/", $modmasses[$k][0] . " ", $modpep); 	
	}
	$next = split(" ", $modpep);
	if(count($next)!==2) {
		echo "Error: have ".count($next)." crosslink modifications in {$modpep} (orig {$copy})\n";
		exit(1);
	}
	return array(join("", $next), strlen($this->stripMods($next[0]))-1);
}

private function stripMods($modpep) {
	return preg_replace( "/\[.*?\]/", "", $modpep);

}

public function getEmail() {
	return $this->iqpir_params['email'];
}


public function computeCombinedCrosslinkQuant($samples = array(), $all = false, $all_norm = false) {
// read the sample ratios file
if($all) {
	foreach($this->sample_offsets as $sample => $offsets) {
		if(count($samples)  > 0 && ! array_key_exists($sample, $samples)) continue;
		foreach($offsets as $offset => $info) {
			foreach($offsets as $offset2 => $info2) {
				if($offset == $offset2) continue;
				$this->computeTwoChannelLogratios($sample, $info[0], $info2[0], $all_norm); //exit(1);
			}
		}
	}
	return;
}

$myfile = fopen($this->iqpir_params['sample_ratios'], "r");
$first = true;
while($line=rtrim(fgets($myfile))){
	if($first) {
		$first = false;
	}
	else {
		$line = preg_replace("/\r|\n/", "", $line);
		$next = split("\t", $line);
		if(count($samples) > 0 && ! array_key_exists($next[0], $samples)) continue;
		$this->computeTwoChannelLogratios($next[0], $next[2], $next[4], $next[6] == "yes"); //exit(1);
	}
}
fclose($myfile);

}
private function overlapPeps($startPosA, $endPosA, $startPosB, $endPosB) {
	$occ = array();
	for($k = $startPosA; $k <= $endPosA; $k++) {
		$occ[$k] = 1;
	}
	for($k = $startPosB; $k <= $endPosB; $k++) {
		if(array_key_exists($k, $occ)) return true;
	}
	return false;
}


public function computeApportionments($samples = array()) {
	foreach($this->sample_offsets as $sample => $offsets) {
		if(count($samples)  > 0 && ! array_key_exists($sample, $samples)) continue;
		if(! array_key_exists($sample, $this->iqpir_params['samples'])) continue;
		$this->computeSampleApportionments($sample);
	}
}

public function computeSampleApportionments($sample) {
$file_suff = $this->deadend ? $this->deadend_subdir : "";
$file = $this->outputdir . $file_suff . "sample." . $sample . ".iqpir.txt";

if(! file_exists($file)) {
	echo "Error: {$file} does not exist\n";
	exit(1);
}
$this->setAnalysisTypes();
$file_hash = hash_file("md5", $file);
$ratio_iontypes = count($this->ratio_iontypes) > 0 ? join(",", array_keys($this->ratio_iontypes)) : "";
$this->setSampleChannels($sample);

$num_peaks = count($this->sample_quantified_isotope_peaks[$sample]);
$num_reporter_peaks = count($this->sample_quantified_reporter_peaks[$sample]);
$data = array(); // xlid to array for contributing ions with logratios
$myfile = fopen($file, "r");
$first = true;
$headers = array();
$num_index = -1;
$denom_index = -1;
$sample_biorepnorms = array();
$existing_norms = array($sample => array());
$protpair_unis = array();
$max_apport_error = $this->max_apport_error; //0.1; //0.25;
$min_apportionment = $this->min_apportionment; //0.02;
$prot_suff = $this->deadend ? "" : "s";
$res_suff = $this->deadend ? "" : "_pair";
$is_intra = array(); // keep track of which cross-links are intra-links
$modmasses = array();
$ion_stats = array();
foreach($this->ratio_iontypes as $ion => $val) $ion_stats[$ion] = array(0, 0);
while($line = rtrim(fgets($myfile))) {
	$next = split("\t", $line);
	if($first) {
		for($k = 0; $k < count($next); $k++) {
			$headers[$next[$k]] = $k;
		}
		$first = false;
	}
	else {
		if($num_index == -1) {
			$modmasses = split(",", $next[$headers["channel_stumpmodmasses"]]);
		
		}
		if($this->intra_only) {
			if(! array_key_exists($next[$headers["xlinkdb-id"]], $is_intra)) {
		
				$prots = split("-", $next[$headers["uniprot"]]);
				$xl_info = split("_", $next[$headers["xlinkdb-id"]]);
				$residues = split("-", $next[$headers["residue_pair"]]);
				$is_intra[$next[$headers["xlinkdb-id"]]] = ! $this->isInterprot($prots[0], $xl_info[0], $residues[0] - $xl_info[1], 
					$prots[1], $xl_info[2], $residues[1] - $xl_info[3]);
			}
			if(! $is_intra[$next[$headers["xlinkdb-id"]]]) {
					continue; // inter-protein
			}
		}
		
		$apports = split(",", $next[$headers["channel_apportionment"]]);
		$next_iontype =  $next[$headers["ion_type"]]; //$this->getIontype($next[$headers["id"]]);
		if(strpos($next_iontype, "Rep")!== false) {
			$next_num_peaks = $num_reporter_peaks;
		}
		else {
			$next_num_peaks = $num_peaks;
		}
		
		$excluded = "FALSE";
		if(array_key_exists($next_iontype, $this->ratio_iontypes)) $ion_stats[$next_iontype][1]++;
		if($max_apport_error > 0 && $next[$headers["apportionment_error"]] > $max_apport_error) {
			$excluded = "APPORT_ERROR";
			echo "Ignoring ion apportionment of {$next[$headers['id']]} of {$next[$headers['xlinkdb-id']]} due to error {$next[$headers['apportionment_error']]} being greater than ".
			$max_apport_error . " max value\n"; //exit(1);
		}
		else if(count($this->ratio_iontypes) > 0 && ! array_key_exists($next_iontype, $this->ratio_iontypes)) {
			$excluded = "EXCL_". strtoupper($next_iontype);
		}
		else if($this->min_sig2noise > 0 && $next[$headers["noise"]] > 0) {

			$total_intensity = $next[$headers["total_intensity"]];
			if($total_intensity < $next_num_peaks * $this->min_sig2noise * $next[$headers["noise"]]) {
				echo "Ignoring ion apportionment of {$next[$headers['id']]} of {$next[$headers['xlinkdb-id']]} due to intensity {$total_intensity} being less than ".($next_num_peaks*$this->min_sig2noise).
				" times noise {$next[$headers['noise']]} (".
				($this->min_sig2noise * $next_num_peaks*$next[$headers["noise"]]).") in scan {$next[$headers['scan']]}\n";
				$excluded = "SIG2NOISE"; // avoid the ints
			}
		}
		if(! array_key_exists($next[$headers["protein{$prot_suff}"]], $protpair_unis)) $protpair_unis[$next[$headers["protein{$prot_suff}"]]] = $next[$headers["uniprot"]];
		if(array_key_exists($next_iontype, $this->ratio_iontypes) && $excluded == "FALSE") $ion_stats[$next_iontype][0]++;
		
		if(! array_key_exists($next[0], $data)) $data[$next[0]] = array();
		$next_entry = array();
		foreach($headers as $name => $col) $next_entry[$name] = $next[$col];
		$next_entry["excluded"] = $excluded;
		array_push($data[$next[0]], $next_entry);
	}
}
fclose($myfile);
echo "Ions used:\n";
foreach($ion_stats as $ion => $info) echo $ion . "\t" . join("\t", $info). "\t" . number_format($info[0]/$info[1], 2, ".", "")."\n";
$headers['numerator_channel_index'] = -1;
$headers['denominator_channel_index'] = -1;
$headers['excluded'] = -1;

// reac in all info and set up to compute logratios for normalization and cross-link accounting
$samplefile = $this->iqpir_params['iqpir_output_dir'] . "sample." . $sample . ".apports.txt";

$chromfile = $this->iqpir_params['iqpir_output_dir'] . "sample." . $sample . ".chroms.txt";
$respairfile = $this->iqpir_params['iqpir_output_dir'] . "sample." . $sample . ".respairs.txt";
$protpairfile = $this->iqpir_params['iqpir_output_dir'] . "sample." . $sample . ".protpairs.txt";
$reshubfile = $this->iqpir_params['iqpir_output_dir'] . "sample." . $sample . ".reshubs.txt";
$numidsfile = $this->iqpir_params['iqpir_output_dir'] . "sample." . $sample . ".numids.txt";
$reporter_data = array();
$resPairMeans = array();
$resHubMeans = array();
$protPairMeans = array();
$header_cols = array("mean_apport", "std_apport", "num_apport", "outliers_apport");

$fsample = fopen($samplefile, "w");
if($this->deadend) fwrite($fsample, "deadend-peptide");
else fwrite($fsample, "cross-link");
fwrite($fsample, "\txlinkdb-id\tprotein{$prot_suff}\tresidue{$res_suff}\tmean_apport\tstd_apport\tnum_apport\toutliers_apport");
fwrite($fsample, "\t{$file_hash}\t{$ratio_iontypes}");
fwrite($fsample, "\n");
$actual_ratios = array();
$first_line = true;
foreach($data as $xl => $info) {
	$apports = array();
	$output = array();
	for($k = 0; $k < count($modmasses); $k++) {
		array_push($apports, array());
	}
	for($k = 0; $k < count($info); $k++) {
		if($info[$k]['excluded']!= "FALSE") continue;
		$next_apports = split(",", $info[$k]['channel_apportionment']);
		for($j = 0; $j < count($modmasses); $j++) array_push($apports[$j], $next_apports[$j]);
	}
	$count = 100;
	for($j = 0; $j < count($modmasses); $j++) {
		if(count($apports[$j]) == 0) continue;
		$nextPeptideAndFragInfo2 = $this->getWeightedAverageRatio($apports[$j], NULL, true, 0, false, false, true, 5);
		if(count($output) == 0) {
			for($i = 0; $i < count($nextPeptideAndFragInfo2); $i++) array_push($output, array());
		}
		for($i = 0; $i < count($nextPeptideAndFragInfo2); $i++) {
			array_push($output[$i], $nextPeptideAndFragInfo2[$i]);
		}
		$count = min($count, $nextPeptideAndFragInfo2[1]);
	}
if(count($output) == 0) continue;
	fwrite($fsample, $xl . "\t" . $info[0]['xlinkdb-id'] . "\t" . $info[0]["protein{$prot_suff}"] . "\t" . $info[0]["residue{$res_suff}"]);
	MyStatistics::normalize($output[0]);
	for($j = 0; $j < count($output); $j++) fwrite($fsample, "\t(" . join(",", $output[$j]) . ")");
	fwrite($fsample, "\n");
	if($count === 0) continue; // nothing more to do
	$nextProts = split("-", $info[0]["protein{$prot_suff}"]);
	$nextres = split("-", $info[0]["residue{$res_suff}"]); //echo "here with {$xl_residuepairs[$xl]}\n"; exit(1);
	$zero_res = $nextres[0] == 0 || (! $this->deadend && $nextres[1] == 0);
	if(! $this->deadend) {
		if(! array_key_exists($nextProts[0] . "\t" . $nextres[0], $resHubMeans)) {
			$resHubMeans[$nextProts[0] . "\t" . $nextres[0]] = array();
			for($j = 0; $j < count($modmasses); $j++) array_push($resHubMeans[$nextProts[0] . "\t" . $nextres[0]], array(array($nextProts[1], $nextres[1], $output[0][$j])));
		}
		else for($j = 0; $j < count($modmasses); $j++) array_push($resHubMeans[$nextProts[0] . "\t" . $nextres[0]][$j], array($nextProts[1], $nextres[1], $output[0][$j]));
		if(! array_key_exists($nextProts[1] . "\t" . $nextres[1], $resHubMeans)) {
			$resHubMeans[$nextProts[1] . "\t" . $nextres[1]] = array();
			for($j = 0; $j < count($modmasses); $j++) array_push($resHubMeans[$nextProts[1] . "\t" . $nextres[1]], array(array($nextProts[0], $nextres[0], $output[0][$j])));
		}
		else for($j = 0; $j < count($modmasses); $j++) array_push($resHubMeans[$nextProts[1] . "\t" . $nextres[1]][$j], array($nextProts[0], $nextres[0], $output[0][$j]));
	}		
	sort($nextProts);
	$nextProts = join("-", $nextProts);
	$nextres = $this->deadend ? $nextres[0] : $nextres[1] . "-" . $nextres[0];
	$indexA = 0;
	$indexB = 2;
	$next_prots = $info[0]["protein{$prot_suff}"];
	$next_res = $info[0]["residue{$res_suff}"];
	if(! $this->deadend && $info[0]["protein{$prot_suff}"]!=$nextProts && array_key_exists($nextProts . "\t" . $nextres, $resPairMeans)) {
		$next_prots  = $nextProts;
		$next_res = $nextres;
		$indexA = 2;
		$indexB = 0; // peptides are in reverse order
	}
	else if(! array_key_exists($info[0]["protein{$prot_suff}"] . "\t" . $info[0]["residue{$res_suff}"], $resPairMeans)) {
		$resPairMeans[$info[0]["protein{$prot_suff}"] . "\t" . $info[0]["residue{$res_suff}"]] = array();
		$resPairPepAs[$info[0]["protein{$prot_suff}"] . "\t" . $info[0]["residue{$res_suff}"]] = array();
		$resPairPepBs[$info[0]["protein{$prot_suff}"] . "\t" . $info[0]["residue{$res_suff}"]] = array();
		for($j = 0; $j < count($modmasses); $j++) {
			array_push($resPairMeans[$info[0]["protein{$prot_suff}"] . "\t" . $info[0]["residue{$res_suff}"]], array());
			array_push($resPairPepAs[$info[0]["protein{$prot_suff}"] . "\t" . $info[0]["residue{$res_suff}"]], array());
			array_push($resPairPepBs[$info[0]["protein{$prot_suff}"] . "\t" . $info[0]["residue{$res_suff}"]], array());
		}
	}
	$next_id = $next_prots . "\t" . $next_res;

	for($j = 0; $j < count($modmasses); $j++) {
		array_push($resPairMeans[$next_id][$j], $output[0][$j]);
	}
	$nextPeps = split("_", $info[0]['xlinkdb-id']);
	$resPairPepAs[$next_id][$nextPeps[$indexA]] = 1;
	if(! $this->deadend) $resPairPepBs[$next_id][$nextPeps[$indexB]] = 1;
	if(! $zero_res) { // make sure not res 0 meaning peptide not found in protein sequence
		if(! array_key_exists($nextProts, $protPairMeans)) {
			$protPairMeans[$nextProts] = array();
			for($j = 0; $j < count($modmasses); $j++) array_push($protPairMeans[$nextProts], array());
		}
		for($j = 0; $j < count($modmasses); $j++) {
			array_push($protPairMeans[$nextProts][$j], $output[0][$j]);
		}
	}
}
fclose($fsample);
echo $samplefile ." written with ".count($protPairMeans)." protpairs to process\n";
	$have_respairs = false;
	foreach($resPairMeans as $id => $channels) {
			if(! $have_respairs) {
				$fres = fopen($respairfile, "w");
				fwrite($fres, "protein{$prot_suff}\tresidue{$res_suff}\tmean_logratio\tstdev_logratio\tnumreps\tpepA_seqs");
				if(! $this->deadend) fwrite($fres, "\tpepB_seqs");
				fwrite($fres, "\t{$file_hash}\t{$ratio_iontypes}\n");
				$have_respairs = true;
			}
			$output = array();
			for($jj= 0; $jj < count($modmasses); $jj++) {
				$next_stats = MyStatistics::stats($channels[$jj], 5, 6);
				if(count($output) == 0) {
					for($i = 0; $i < count($next_stats); $i++) array_push($output, array());
				}
				for($i = 0; $i < count($next_stats); $i++) {
					array_push($output[$i], $next_stats[$i]);
				}
			}
		ksort($resPairPepAs[$id]);
		if(! $this->deadend) ksort($resPairPepBs[$id]);
		fwrite($fres, $id); 
		MyStatistics::normalize($output[0]);
		for($j = 0; $j < count($output); $j++) fwrite($fres, "\t(" . join(",", $output[$j]) . ")");
		fwrite($fres, "\t" . join(",", array_keys($resPairPepAs[$id])));
		if(! $this->deadend) fwrite($fres, "\t".join(",", array_keys($resPairPepBs[$id])));
		fwrite($fres, "\n");
	}
	if($have_respairs) {
		echo "Files {$samplefile} and {$respairfile} written for sample {$sample}\n";
		fclose($fres);
	}
	else echo "Files {$samplefile} written for sample {$sample}\n"; //exit(1);
	$first = true;
	$prot_written = false;
	$remove_protpair_outliers = true;
	foreach($protPairMeans as $id => $channels) {
		if(! array_key_exists($id, $protpair_unis)) continue;
			if($first) {
				$fprot = fopen($protpairfile, "w");
				fwrite($fprot, "protein{$prot_suff}\tuniprot{$res_suff}\tmean_logratio\tstdev_logratio\tnumreps");
				if($remove_protpair_outliers) fwrite($fprot, "\toutliers");
				fwrite($fprot, "\t{$file_hash}\t{$ratio_iontypes}\n");
				$first = false;
				$prot_written = true;
			}
			$output = array();
			for($jj = 0; $jj < count($modmasses); $jj++) {
				if(! $remove_protpair_outliers) $nextstats = MyStatistics::stats($channels[$jj], 5, 6);
				else $nextstats = $this->getWeightedAverageRatio($channels[$jj], NULL, true, 0, false, false, true, 5);
				if(count($output) == 0) {
					for($i = 0; $i < count($nextstats); $i++) array_push($output, array());
				}
				for($i = 0; $i < count($nextstats); $i++) {
					array_push($output[$i], $nextstats[$i]);
				}
			}
			fwrite($fprot, $id . "\t" . $protpair_unis[$id]);
			// normalize first array
			MyStatistics::normalize($output[0]);
			for($j = 0; $j < count($output); $j++) fwrite($fprot, "\t(" . join(",", $output[$j]) . ")");
			fwrite($fprot, "\n");
	}//exit(1);
	if($prot_written) {
		echo "Files {$protpairfile} written for sample {$sample}\n";
		fclose($fprot);
	}
	$first = true;
	uksort($resHubMeans, array($this, "sortByFirstAndSecondIndex"));
	$hub_written = false;
	foreach($resHubMeans as $id => $channels) {
		if($first) {
			$fres = fopen($reshubfile, "w");
			fwrite($fres, "protein\tresidue\tmean_logratio\tstdev_logratio\tnumreps\tpartner_protein\tpartner_residue\tmean_logratio\t{$file_hash}\t{$ratio_iontypes}\n");
			$first = false;
			$hub_written = true;
		}
			// get the stats
		$output = array();
		for($jj = 0; $jj < count($modmasses); $jj++) {
			usort($channels[$jj], array($this, "sortByThirdIndex")); //exit(1);
			$hubratios = array();
			for($k = 0; $k < count($channels[$jj]); $k++) array_push($hubratios, $channels[$jj][$k][2]);

			$nextstats = MyStatistics::stats($hubratios, 2, 3);
			if(count($output) == 0) {
				for($i = 0; $i < count($nextstats); $i++) array_push($output, array());
			}
			for($i = 0; $i < count($nextstats); $i++) {
				array_push($output[$i], $nextstats[$i]);
			}
		}
		
		fwrite($fres, $id);
		MyStatistics::normalize($output[0]);
		for($j = 0; $j < count($output); $j++) fwrite($fres, "\t(" . join(",", $output[$j]) . ")");
				
		fwrite($fres, "\n");
	}
	if($hub_written) {
		echo "File {$reshubfile} written for sample {$sample}\n";
		fclose($fres);
	}
	
}
	
private function isInterprot($proA, $pepA, $startPosA, $proB, $pepB, $startPosB) {
	if($pepB == "LOOPLINK") return false;
	if($proA!=$proB) return true;
	if(strpos($pepA, "[")!==false) $pepA = MyPepXML::stripMods($pepA);
	if(strpos($pepB, "[")!==false) $pepB = MyPepXML::stripMods($pepB);
	return $this->overlapPeps($startPosA, $startPosA + strlen($pepA) - 1, $startPosB, $startPosB + strlen($pepB) - 1);
}


public function computeTwoChannelLogratios($sample, $numerator_stumpmodmass, $denominator_stumpmodmass, $normalize = false) {
$file_suff = $this->deadend ? $this->deadend_subdir : "";
$file = $this->outputdir . $file_suff . "sample." . $sample . ".iqpir.txt";
$normfile = $this->outputdir . $file_suff . "sample_biorep_normfactors.txt";
$normfile = $this->iqpir_params['iqpir_output_dir'] . "sample_biorep_normfactors.txt";
if(! file_exists($file)) {
	echo "Error: {$file} does not exist\n";
	exit(1);
}
$this->setAnalysisTypes();
$file_hash = hash_file("md5", $file);
$ratio_iontypes = count($this->ratio_iontypes) > 0 ? join(",", array_keys($this->ratio_iontypes)) : "";
$add_zscore2ratiofile = true;
$this->setSampleChannels($sample);
$num_offset = -1;
$denom_offset = -1;
foreach($this->offsets as $offset => $info) {
	if($num_offset == -1 && ($info[0] == $numerator_stumpmodmass || 
		(false && array_key_exists($numerator_stumpmodmass, $this->default_6plex_modmasses) && 
		$info[0] == $this->default_6plex_modmasses[$numerator_stumpmodmass])
		|| (array_key_exists($numerator_stumpmodmass, $this->channel_modmasses) && 
			$info[0] == $this->channel_modmasses[$numerator_stumpmodmass])
		
		)) $num_offset = $offset;
	if($denom_offset == -1 && ($info[0] == $denominator_stumpmodmass || 
		(false && array_key_exists($denominator_stumpmodmass, $this->default_6plex_modmasses) && 
		$info[0] == $this->default_6plex_modmasses[$denominator_stumpmodmass])
		|| (array_key_exists($denominator_stumpmodmass, $this->channel_modmasses) && 
			$info[0] == $this->channel_modmasses[$denominator_stumpmodmass])
		)) $denom_offset = $offset;

}
if($num_offset == -1 || $denom_offset == -1) {
	echo "Error no offset for num {$num_offset} or denom {$denom_offset} for num {$numerator_stumpmodmass} and denom {$denominator_stumpmodmass} modmasses\n";
	exit(1);
}
// set up which obs ints to use for sig2noise
$index = 0;
$use_obs_peaks = array();
foreach($this->sample_quantified_isotope_peaks[$sample] as $pk_ind => $val) {
	if($pk_ind >= $num_offset && $pk_ind <= $num_offset + 4) $use_obs_peaks[$index] = 1;
	if($pk_ind >= $denom_offset && $pk_ind <= $denom_offset + 4) $use_obs_peaks[$index] = 1;
	$index++;
}
$index = 0;
$use_obs_reporter_peaks = array();
foreach($this->sample_quantified_reporter_peaks[$sample] as $pk_ind => $val) {
	if($pk_ind >= -2 * $num_offset && $pk_ind <= -2 * $num_offset + 4) $use_obs_reporter_peaks[$index] = 1;
	if($pk_ind >= -2 * $denom_offset && $pk_ind <= -2 * $denom_offset + 4) $use_obs_reporter_peaks[$index] = 1;
	$index++;
}
$data = array(); // xlid to array for contributing ions with logratios
$myfile = fopen($file, "r");
$first = true;
$headers = array();
$num_index = -1;
$denom_index = -1;
$sample_biorepnorms = array();
$existing_norms = array($sample => array());
$protpair_unis = array();
$max_apport_error = 0.1; //0.25;
$include_reporters = false; 
$min_apportionment = $this->min_apportionment; //0.025;
$prot_suff = $this->deadend ? "" : "s";
$res_suff = $this->deadend ? "" : "_pair";
$is_intra = array(); // keep track of which cross-links are intra-links
while($line = rtrim(fgets($myfile))) {
	$next = split("\t", $line);
	if($first) {
		for($k = 0; $k < count($next); $k++) {
			$headers[$next[$k]] = $k;
		}
		$first = false;
	}
	else {
		if($num_index == -1) {
			$modmasses = split(",", $next[$headers["channel_stumpmodmasses"]]);
			for($k = 0; $k < count($modmasses); $k++) {
				if($modmasses[$k] == $numerator_stumpmodmass || 
					(false && array_key_exists($numerator_stumpmodmass, $this->default_6plex_modmasses) && 
					$modmasses[$k] == $this->default_6plex_modmasses[$numerator_stumpmodmass])
					|| (array_key_exists($numerator_stumpmodmass, $this->channel_modmasses) && 
					$modmasses[$k] == $this->channel_modmasses[$numerator_stumpmodmass])
					
					) $num_index = $k;
				else if($modmasses[$k] == $denominator_stumpmodmass || 
					(false && array_key_exists($denominator_stumpmodmass, $this->default_6plex_modmasses) && 
					$modmasses[$k] == $this->default_6plex_modmasses[$denominator_stumpmodmass])
					|| (array_key_exists($denominator_stumpmodmass, $this->channel_modmasses) && 
						$modmasses[$k] == $this->channel_modmasses[$denominator_stumpmodmass])
					) $denom_index = $k;
			}
		}
		if($this->intra_only) {
			if(! array_key_exists($next[$headers["xlinkdb-id"]], $is_intra)) {
		
				$prots = split("-", $next[$headers["uniprot"]]);
				$xl_info = split("_", $next[$headers["xlinkdb-id"]]);
				$residues = split("-", $next[$headers["residue_pair"]]);
				$is_intra[$next[$headers["xlinkdb-id"]]] = ! $this->isInterprot($prots[0], $xl_info[0], $residues[0] - $xl_info[1], 
					$prots[1], $xl_info[2], $residues[1] - $xl_info[3]);
			}
			if(! $is_intra[$next[$headers["xlinkdb-id"]]]) {
					continue; // inter-protein
			}
		}
		
		
		
		
		$apports = split(",", $next[$headers["channel_apportionment"]]);
		$next_iontype =  $next[$headers["ion_type"]]; //$this->getIontype($next[$headers["id"]]);
		if(strpos($next_iontype, "Rep")!== false) {
			$next_use_obs_peaks = $use_obs_reporter_peaks;
		}
		else {
			$next_use_obs_peaks = $use_obs_peaks;
		}
		$excluded = "FALSE";
		if($apports[$num_index] == 0 && $apports[$denom_index] == 0) continue; 
		
		if($min_apportionment > 0 && $apports[$num_index] < $min_apportionment && $apports[$denom_index] < $min_apportionment) {
			$excluded = "ZERO_APPORT";
		}
		else if($max_apport_error > 0 && $next[$headers["apportionment_error"]] > $max_apport_error) {
			$excluded = "APPORT_ERROR";
		}
		else if(count($this->ratio_iontypes) > 0 && ! array_key_exists($next_iontype, $this->ratio_iontypes)) {
			$excluded = "EXCL_". strtoupper($next_iontype);
		}
		else if($this->min_sig2noise > 0 && $next[$headers["noise"]] > 0) {

			$total_intensity = 0;
			$pks = split(",", $next[$headers["obs_relints"]]);
			foreach($next_use_obs_peaks as $pk_ind => $val) {
				$total_intensity += $pks[$pk_ind];
			}
			$next_num_peaks = count($next_use_obs_peaks);
			$total_intensity = number_format($total_intensity * $next[$headers["total_intensity"]], 0, ".", "");
			if($total_intensity < $next_num_peaks * $this->min_sig2noise * $next[$headers["noise"]]) {
				echo "Ignoring ion ratio of {$next[$headers['id']]} of {$next[$headers['xlinkdb-id']]} due to intensity {$total_intensity} being less than ".($next_num_peaks*$this->min_sig2noise).
				" times noise {$next[$headers['noise']]} (".
				($this->min_sig2noise * $next_num_peaks*$next[$headers["noise"]]).") in scan {$next[$headers['scan']]}\n";
				$excluded = "SIG2NOISE"; // avoid the ints
			}
		}
		$next_ratio = 0;
		if($apports[$num_index] == 0) $next_ratio = -1 * $this->max_log2ratio;
		else if($apports[$denom_index] == 0) $next_ratio = $this->max_log2ratio;
		else $next_ratio = number_format(log($apports[$num_index]/$apports[$denom_index])/log(2), 2, ".", "");
		$next_orig_ratio = $next_ratio;
		$next_ratio = min($this->max_log2ratio, $next_ratio);
		$next_ratio = max(-1 * $this->max_log2ratio, $next_ratio);
		if(! array_key_exists($next[$headers["protein{$prot_suff}"]], $protpair_unis)) $protpair_unis[$next[$headers["protein{$prot_suff}"]]] = $next[$headers["uniprot"]];
		if(! $this->deadend) {
			$next_p = split("-", $next[$headers["protein{$prot_suff}"]]);
			$next_u = split("-", $next[$headers["uniprot"]]);
			if(count($next_p) == 2 && count($next_u) == 2 && 
				! array_key_exists($next_p[1] . "-" . $next_p[0], $protpair_unis)) 
					$protpair_unis[$next_p[1] . "-" . $next_p[0]] = $next_u[1] . "-" . $next_u[0];
		}

		if($normalize) {
			$next_norm_id = $next[$headers["biorep"]] . "-" . $numerator_stumpmodmass . "-" . $denominator_stumpmodmass;
			if(! array_key_exists($next_norm_id, $sample_biorepnorms)) {
				$sample_biorepnorms[$next_norm_id] = array();
				$existing_norms[$sample][$next_norm_id] = array();
			}
			array_push($sample_biorepnorms[$next_norm_id], $next_ratio);
		}
		
		if(! array_key_exists($next[0], $data)) $data[$next[0]] = array();
		$next_entry = array();
		foreach($headers as $name => $col) $next_entry[$name] = $next[$col];
		$next_entry["logratio"] = $next_ratio;
		$next_entry["excluded"] = $excluded;
		$next_entry["numerator_channel_index"] = $num_index + 1;
		$next_entry["denominator_channel_index"] = $denom_index + 1;
		$next_entry['orig_logratio'] = $next_orig_ratio;
		array_push($data[$next[0]], $next_entry);
	}
}
fclose($myfile);
$headers['numerator_channel_index'] = -1;
$headers['denominator_channel_index'] = -1;
$headers['excluded'] = -1;
if($normalize && $this->use_current_normalization) {
	$this->setUseCurrentNorm();
	$existing_norms = $this->iqpir_params['samples'];
}
if($normalize && ! $this->use_current_normalization) {
    foreach($sample_biorepnorms as $biorep_num_denom => $ratios) {
			echo "biorep {$biorep_num_denom} with ".count($ratios)." ratios\n";
			
		if($this->normalize_with_mean) $this->iqpir_params['samples'][$sample][$biorep_num_denom] = number_format(-1 * array_sum($ratios)/count($ratios), 2, ".", "");
		else $this->iqpir_params['samples'][$sample][$biorep_num_denom] = number_format(-1 * MyStatistics::median($ratios), 2, ".", "");
		$existing_norms[$sample][$biorep_num_denom] = $this->iqpir_params['samples'][$sample][$biorep_num_denom];
	}
	if(file_exists($normfile)) {
		$first = true;
		$myfile = fopen($normfile, "r");
		while($line = rtrim(fgets($myfile))) {
			if($first) {
				$first = false;
				continue;
			}
			$next = split("\t", $line);
			if($next[0] != $sample || ! array_key_exists($next[1] . "-" . $next[2] . "-" . $next[3], $existing_norms[$next[0]])) {
				if(! array_key_exists($next[0], $existing_norms)) $existing_nroms[$next[0]] = array();
				 $existing_norms[$next[0]][$next[1] . "-" . $next[2] . "-" . $next[3]] = $next[4];
				echo "Adding pre-existing normalization additive factor {$next[4]} for sample {$next[0]} biorep ".$next[1] . "-" . $next[2] . "-" . $next[3]."\n";
			}
		}
		fclose($myfile);
	}
	$fnorm = fopen($normfile, "w");
	fwrite($fnorm, "sample\tbiorep\tnumerator_stumpmodmass\tdenominator_stumpmodmass\tnorm_add_to_log\n");
		foreach($existing_norms as $norm_sample => $bioreps) {
			foreach($bioreps as $biorep => $norm) {
				$biorep = preg_replace("/\-/", "\t", $biorep);
				fwrite($fnorm, $norm_sample . "\t" . $biorep . "\t" . $norm ."\n");
			}
		}
		fclose($fnorm);	
		echo "Normalization file {$normfile} written for all samples: ".join(",", array_keys($this->iqpir_params['samples']))."\n";
} // if normalize


if($this->dominant_fragid_maxpval < 1) {
	foreach($data as $xl => $info) {
		if(strpos($xl, "_LOOPLINK")===false) continue;
		$fragids = array();
		$all_fragids = array();
		$fragids_B = array();
		$all_fragids_B = array();
		for($k = 0; $k < count($info); $k++) {
			if($info[$k]['ion_type'] != "Fragment") continue;
			if(strpos($info[$k]['id'], "A")!==false) { // A or AB
				$all_fragids[$info[$k]['id']] = 1;
				if($info[$k]['excluded']!= "FALSE") continue;
				if(! array_key_exists($info[$k]['id'], $fragids)) $fragids[$info[$k]['id']] = 0;
				$fragids[$info[$k]['id']]++;
			}
			else { // B
				$all_fragids_B[$info[$k]['id']] = 1;
				if($info[$k]['excluded']!= "FALSE") continue;
				if(! array_key_exists($info[$k]['id'], $fragids_B)) $fragids_B[$info[$k]['id']] = 0;
				$fragids_B[$info[$k]['id']]++;
			}
		} // next ion
		$next_tot = array_sum(array_values($fragids));
		if($next_tot >= 10 && count($all_fragids) >= 5) {
			arsort($fragids);
			$next_ions = array_keys($fragids);
			if($fragids[$next_ions[0]]/$next_tot >= 0.75) {
				$next_pep = split("_", $xl);
				$next_pep = MyPepXML::stripmods($next_pep[0]);
				$next_index = split(":", $next_ions[0]);
				$next_index = split("_", $next_index[1]);
				$next_index = substr($next_index[0], 1);
				if(strpos($next_ions[0], "y") !== false) {
					if($next_pep[strlen($next_pep) - $next_index] == "P") continue;
				}
				else if(strpos($next_ions[0], "b") !== false) {
					if($next_index < strlen($next_pep) && $next_pep[$next_index] == "P") 
						continue;
				}
				// still here, must mark the excluded
				// see if ratios of dom ion are different from other non-false
				$dom_ratios = array();
				$other_ratios = array();
				for($k = 0; $k < count($info); $k++) {
					if($info[$k]['excluded']!= "FALSE") continue;
					$norm = $normalize &&
						array_key_exists($info[$k]['biorep'] . "-" . $numerator_stumpmodmass . "-" . $denominator_stumpmodmass, $this->iqpir_params['samples'][$sample]) ? 
						$this->iqpir_params['samples'][$sample][$info[$k]['biorep'] . "-" . $numerator_stumpmodmass . "-" . $denominator_stumpmodmass] : 0;
					$next_ratio = $info[$k]['logratio'] + $norm;
					$next_ratio = min($this->max_log2ratio, $next_ratio);
					$next_ratio = max(-1 * $this->max_log2ratio, $next_ratio);


					if($info[$k]['ion_type'] == "Fragment" && $info[$k]['id']== $next_ions[0]) array_push($dom_ratios, $next_ratio);
					else array_push($other_ratios, $next_ratio);
				}
				$dom_stats = MyStatistics::stats($dom_ratios, 2, 3);
				if(count($other_ratios) > 0) {
					if(count($other_ratios) == 1) {
					$dom_stats[1] = max($dom_stats[1], 0.01);
						$next_z = ($other_ratios[0]-$dom_stats[0])/$dom_stats[1];
							if(abs($next_z) < 1.2) continue;
					}
					else {
				
						$other_stats = MyStatistics::stats($other_ratios, 2, 3);	
						$next_t = MyStatistics::oneWayAnovaAlt($dom_stats, $other_stats, "anova_pvals.txt");
						if($next_t[2] >= $this->dominant_fragid_maxpval) continue;
					}
				}	
				for($k = 0; $k < count($info); $k++) {
					if($info[$k]['ion_type'] != "Fragment") continue;
					if($info[$k]['excluded']!= "FALSE") continue;
					if($info[$k]['id']!= $next_ions[0]) continue;
					$data[$xl][$k]['excluded'] = "DOMINANT_FRAGID";
					//echo "Assigned this {$info[$k]['id']} to {$info[$k]['excluded']} for {$info[$k]['rawfile']} {$info[$k]['scan']}\n"; exit(1);
				}
			} // dominant id
		} // enough total ions to analyze		
		$next_tot = array_sum(array_values($fragids_B));
		if($next_tot >= 10 && count($all_fragids_B) >= 5) {
			arsort($fragids_B);
			$next_ions = array_keys($fragids_B);
			if($fragids_B[$next_ions[0]]/$next_tot >= 0.75) {
				$next_pep = split("_", $xl);
				$next_pep = MyPepXML::stripmods($next_pep[1]);
				$next_index = split(":", $next_ions[0]);
				$next_index = split("_", $next_index[1]);
				$next_index = substr($next_index[0], 1);
				if(strpos($next_ions[0], "y") !== false) {
					if($next_pep[strlen($next_pep) - $next_index] == "P") continue;
				}
				else if(strpos($next_ions[0], "b") !== false) {
					if($next_index < strlen($next_pep) && $next_pep[$next_index] == "P") 
						continue;
				}
				for($k = 0; $k < count($info); $k++) {
					if($info[$k]['ion_type'] != "Fragment") continue;
					if($info[$k]['excluded']!= "FALSE") continue;
					if($info[$k]['id']!= $next_ions[0]) continue;
					$data[$xl][$k]['excluded'] = "DOMINANT_FRAGID";
					//echo "Assigned this {$info[$k]['id']} to {$info[$k]['excluded']} for {$info[$k]['rawfile']} {$info[$k]['scan']}\n";
				}
			} // dominant id
		} // enough total ions to analyze		
	} // next xl

} // if analyze used fragment ids

// reac in all info and set up to compute logratios for normalization and cross-link accounting
$suffix = ".{$numerator_stumpmodmass}-{$denominator_stumpmodmass}";
$samplefile = $this->iqpir_params['iqpir_output_dir'] . "sample." . $sample . $suffix . ".ratios.txt";

$numerator_reportermass = "";
$denominator_reportermass = "";
foreach($this->sample_offsets[$sample] as $offset => $info) {
	if($info[0] == $numerator_stumpmodmass || 
		(array_key_exists($numerator_stumpmodmass, $this->default_6plex_modmasses) && 
		$info[0] == $this->default_6plex_modmasses[$numerator_stumpmodmass])) $numerator_reportermass = $info[1];
	else if($info[0] == $denominator_stumpmodmass || 
		(array_key_exists($denominator_stumpmodmass, $this->default_6plex_modmasses) && 
		$info[0] == $this->default_6plex_modmasses[$denominator_stumpmodmass])) $denominator_reportermass = $info[1];
}

$chromfile = $this->iqpir_params['iqpir_output_dir'] . "sample." . $sample . $suffix .".chroms.txt";
$respairfile = $this->iqpir_params['iqpir_output_dir'] . "sample." . $sample . $suffix . ".respairs.txt";
$protpairfile = $this->iqpir_params['iqpir_output_dir'] . "sample." . $sample . $suffix . ".protpairs.txt";
$reshubfile = $this->iqpir_params['iqpir_output_dir'] . "sample." . $sample . $suffix . ".reshubs.txt";
$numidsfile = $this->iqpir_params['iqpir_output_dir'] . "sample." . $sample . $suffix . ".numids.txt";
$reporter_data = array();
$resPairMeans = array();
$resHubMeans = array();
$protPairMeans = array();

$fsample = fopen($samplefile, "w");
if($this->deadend) fwrite($fsample, "deadend-peptide");
else fwrite($fsample, "cross-link");
fwrite($fsample, "\txlinkdb-id\tprotein{$prot_suff}\tresidue{$res_suff}\tactual_ratio\tlogratio_mean\tlogratio_num\tlogratio_std\tlogratio_tstat\tlogratio_df\tlogatio_pval\tlogratio_outliers");
fwrite($fsample, "\tnumerator_stumpmodmass\tdenominator_stumpmodmass");
if(! $add_zscore2ratiofile) fwrite($fsample, "\t{$file_hash}\t{$ratio_iontypes}");
fwrite($fsample, "\n");
$fchrom = fopen($chromfile, "w");
$normratios = array();


if($this->deadend) fwrite($fchrom, "deadend-peptide");
$actual_ratios = array();
$first_line = true;
foreach($data as $xl => $info) {
	$ratios = array();
	$orig_ratios = array(); // for std if ratio over
	$optimal_ratios = array();
	for($k = 0; $k < count($info); $k++) {
		if($info[$k]['excluded']!= "FALSE") continue;
		$norm = $normalize &&
			array_key_exists($info[$k]['biorep'] . "-" . $numerator_stumpmodmass . "-" . $denominator_stumpmodmass, $this->iqpir_params['samples'][$sample]) ? 
				$this->iqpir_params['samples'][$sample][$info[$k]['biorep'] . "-" . $numerator_stumpmodmass . "-" . $denominator_stumpmodmass] : 0;
		$next_ratio = $info[$k]['logratio'] + $norm;
		$next_ratio = min($this->max_log2ratio, $next_ratio);
		$next_ratio = max(-1 * $this->max_log2ratio, $next_ratio);
		array_push($ratios, $next_ratio);
		array_push($orig_ratios, $info[$k]['orig_logratio']);
	}
	if($include_reporters && array_key_exists($xl, $reporter_data)) {
		for($r = 0; $r < count($reporter_data[$xl]); $r++) {
			$next_ratio = $reporter_data[$xl][$r][0] + $norm;
			$next_ratio = min($this->max_log2ratio, $next_ratio);
			$next_ratio = max(-1 * $this->max_log2ratio, $next_ratio);
			array_push($ratios, $next_ratio);
		}
	}
		
		
	$nextPeptideAndFragInfo2 = $this->getWeightedAverageRatio($ratios, NULL, true, 0, false, true, true); //, 2, false, count($ion_ids));
	$nextPeptideAndFragInfo2_orig = $this->getWeightedAverageRatio($orig_ratios, NULL, true, 0, false, true, true); //, 2, false, count($ion_ids));


	$next_actual = array_key_exists($info[0]["protein{$prot_suff}"], $actual_ratios) ? $actual_ratios[$info['proteins']] : 0;
	$outliers = $nextPeptideAndFragInfo2[count($nextPeptideAndFragInfo2)-1] == "" ? array() : split(",", $nextPeptideAndFragInfo2[count($nextPeptideAndFragInfo2)-1]);
$outliers = array_flip($outliers);
$outliers = array_keys($outliers);
	for($k = 0; $k < count($outliers); $k++) {
		$logratios = array();
		for($j = 0; $j < count($info); $j++) {
			if($info[$j]['excluded']!= "FALSE") continue;
			$norm = $normalize &&
				array_key_exists($info[$j]['biorep'] . "-" . $numerator_stumpmodmass . "-" . $denominator_stumpmodmass, $this->iqpir_params['samples'][$sample]) ? 
					$this->iqpir_params['samples'][$sample][$info[$j]['biorep'] . "-" . $numerator_stumpmodmass . "-" . $denominator_stumpmodmass] : 0;


			$next_val = $info[$j]['logratio']; //$outliers[$k]; // - $this->iqpir_params['samples'][$sample][$info[$j]['biorep'] . "-" . $numerator_stumpmodmass . "-" . $denominator_stumpmodmass];
			$next_val = min($this->max_log2ratio, $next_val);
			$next_val = max(-1 * $this->max_log2ratio, $next_val);
			$next_ratio = $info[$j]['logratio'] + $norm; //$this->iqpir_params['samples'][$sample][$info[$j]['biorep'] . "-" . $numerator_stumpmodmass . "-" . $denominator_stumpmodmass][0];
			$next_ratio = min($this->max_log2ratio, $next_ratio);
			$next_ratio = max(-1 * $this->max_log2ratio, $next_ratio);
			array_push($logratios, $next_ratio);
			//echo "Have {$info[$j]['logratio']} and {$norm}\n";
if(abs($outliers[$k] - $next_ratio) < 0.01) {
				$info[$j]['excluded'] = "OUTLIER";
				$found = true;
			}
		}
		if(! $found && $include_reporters && array_key_exists($xl, $reporter_data)) {
			for($r = 0; $r < count($reporter_data[$xl]); $r++) {
				if(abs($outliers[$k] - $this->iqpir_params['samples'][$sample][$reporter_data[$xl][$r][3]][0] 
					- $reporter_data[$xl][$r][0]) < 0.01) {
					array_push($reporter_data[$xl][$r], "OUTLIER");
					$found = true;
				}
				array_push($logratios, $reporter_data[$xl][$r][0] + $this->iqpir_params['samples'][$sample][$reporter_data[$xl][$r][3]][0]);
			}
		}
		if(! $found) {
			sort($logratios);
			echo "Error: could not find outlier normalized value {$outliers[$k]} among logratios of {$xl}: ".join(",", $logratios)."\n";
			exit(1);
		}
	}
	$next_actual = array_key_exists($info[0]["protein{$prot_suff}"], $actual_ratios) ? $actual_ratios[$info[0]["protein{$prot_suff}"]] : "";
	$num_RH_SH_ids[$xl]['log2ratio'] = $nextPeptideAndFragInfo2[0];
	$nextPeptideAndFragInfo2[2] = max($nextPeptideAndFragInfo2[2], $nextPeptideAndFragInfo2_orig[2]);

	fwrite($fsample, $xl . "\t" . $info[0]['xlinkdb-id'] . "\t" . $info[0]["protein{$prot_suff}"] . "\t" . $info[0]["residue{$res_suff}"] . "\t" . $next_actual . "\t" . join("\t", $nextPeptideAndFragInfo2));
	fwrite($fsample, "\t" . (array_key_exists($numerator_stumpmodmass, $this->default_6plex_modmasses) ? 
		$this->default_6plex_modmasses[$numerator_stumpmodmass] : $numerator_stumpmodmass) . "\t" .
		(array_key_exists($denominator_stumpmodmass, $this->default_6plex_modmasses) ? 
		$this->default_6plex_modmasses[$denominator_stumpmodmass] : $denominator_stumpmodmass));

	fwrite($fsample, "\n");
	if($normalize && count($this->sample_bioreps[$sample]) == 1 && $nextPeptideAndFragInfo2[1] > 1) 
		array_push($normratios, $nextPeptideAndFragInfo2[0] - $norm);
	if($nextPeptideAndFragInfo2[1] === 0) continue; // nothing more to do
	$nextProts = split("-", $info[0]["protein{$prot_suff}"]);
	$nextres = split("-", $info[0]["residue{$res_suff}"]); //echo "here with {$xl_residuepairs[$xl]}\n"; exit(1);
	$zero_res = $nextres[0] == 0 || (! $this->deadend && $nextres[1] == 0);
	if(! $this->deadend) {
		if(! array_key_exists($nextProts[0] . "\t" . $nextres[0], $resHubMeans)) $resHubMeans[$nextProts[0] . "\t" . $nextres[0]] = array();
		array_push($resHubMeans[$nextProts[0] . "\t" . $nextres[0]], array($nextProts[1], $nextres[1], $nextPeptideAndFragInfo2[0]));
		if(! array_key_exists($nextProts[1] . "\t" . $nextres[1], $resHubMeans)) $resHubMeans[$nextProts[1] . "\t" . $nextres[1]] = array();
		if($nextProts[0]!=$nextProts[1] || $nextres[0] != $nextres[1]) array_push($resHubMeans[$nextProts[1] . "\t" . $nextres[1]], 
			array($nextProts[0], $nextres[0], $nextPeptideAndFragInfo2[0]));
	}		
	sort($nextProts);
	$nextProts = join("-", $nextProts);
	$nextres = $this->deadend ? $nextres[0] : $nextres[1] . "-" . $nextres[0];
	$indexA = 0;
	$indexB = 2;
	$next_prots = $info[0]["protein{$prot_suff}"];
	$next_res = $info[0]["residue{$res_suff}"];
	if(! $this->deadend && $info[0]["protein{$prot_suff}"]!=$nextProts && array_key_exists($nextProts . "\t" . $nextres, $resPairMeans)) {
		$next_prots  = $nextProts;
		$next_res = $nextres;
		$indexA = 2;
		$indexB = 0; // peptides are in reverse order
	}
	else if(! array_key_exists($info[0]["protein{$prot_suff}"] . "\t" . $info[0]["residue{$res_suff}"], $resPairMeans)) {
		$resPairMeans[$info[0]["protein{$prot_suff}"] . "\t" . $info[0]["residue{$res_suff}"]] = array();
		$resPairPepAs[$info[0]["protein{$prot_suff}"] . "\t" . $info[0]["residue{$res_suff}"]] = array();
		$resPairPepBs[$info[0]["protein{$prot_suff}"] . "\t" . $info[0]["residue{$res_suff}"]] = array();
	}
	$next_id = $next_prots . "\t" . $next_res;
	array_push($resPairMeans[$next_id], $nextPeptideAndFragInfo2[0]);
	$nextPeps = split("_", $info[0]['xlinkdb-id']);
	$resPairPepAs[$next_id][$nextPeps[$indexA]] = 1;
	if(! $this->deadend) $resPairPepBs[$next_id][$nextPeps[$indexB]] = 1;

	if(! $zero_res) { // && $nextPeptideAndFragInfo2[1] > 1) { //$nextres[0] != 0 && $nextres[1] != 0) { // make sure not res 0 meaning peptide not found in protein sequence
		if(! array_key_exists($nextProts, $protPairMeans)) $protPairMeans[$nextProts] = array();
		array_push($protPairMeans[$nextProts], $nextPeptideAndFragInfo2[0]);
	}
	for($k = 0; $k < count($info); $k++) {
		$norm = 
			array_key_exists($info[$k]['biorep'] . "-" . $numerator_stumpmodmass . "-" . $denominator_stumpmodmass, $this->iqpir_params['samples'][$sample]) ? 
				$this->iqpir_params['samples'][$sample][$info[$k]['biorep'] . "-" . $numerator_stumpmodmass . "-" . $denominator_stumpmodmass] : 0;
		$next_th = split(",", $info[$k]['theor_relints']);
		$next_ob = split(",", $info[$k]['obs_relints']);
		$current_ratio = number_format(pow(2, $info[$k]['logratio']), 2, ".", "");
		$norm_logratio = $info[$k]['logratio']+$norm;
		$norm_logratio = min($this->max_log2ratio, $norm_logratio);
		$norm_logratio = max(-1 * $this->max_log2ratio, $norm_logratio);
		$norm_ratio = number_format(pow(2, $norm_logratio), 2, ".", "");
		if($first_line) {
			fwrite($fchrom, join("\t", array_keys($info[$k])));
			fwrite($fchrom, "\tnorm_logratio\tnorm_ratio");
			fwrite($fchrom, "\t{$file_hash}\t{$ratio_iontypes}\n");
			$first_line = false;

		}
		fwrite($fchrom,
			join("\t", array_values($info[$k]))); //.
		fwrite($fchrom, "\t" . $norm_logratio . "\t" . $norm_ratio);
		fwrite($fchrom, "\n");//exit(1);
	}
}
fclose($fchrom);
fclose($fsample);
echo $samplefile ." and {$chromfile} written\n";
if($normalize && count($normratios) >= 500) {
	$next_norm_ratio2 = number_format(-1 * MyStatistics::median($normratios), 2, ".", "");
	$lone_biorep = array_keys($this->sample_bioreps[$sample]);
	$next_norm_id = $lone_biorep[0] . "-" . $numerator_stumpmodmass . "-" . $denominator_stumpmodmass;
	if($this->normalize_with_mean) $this->iqpir_params['samples'][$sample][$next_norm_id] = number_format(-1 * array_sum($normratios)/count($normratios), 2, ".", "");
	else $this->iqpir_params['samples'][$sample][$next_norm_id] = number_format(-1 * MyStatistics::median($normratios), 2, ".", "");
	$existing_norms[$sample][$next_norm_id] = $this->iqpir_params['samples'][$sample][$next_norm_id];
	$fnorm = fopen($normfile, "w");
	fwrite($fnorm, "sample\tbiorep\tnumerator_stumpmodmass\tdenominator_stumpmodmass\tnorm_add_to_log\n");
	foreach($existing_norms as $norm_sample => $bioreps) {
		foreach($bioreps as $biorep => $norm) {
			$biorep = preg_replace("/\-/", "\t", $biorep);
			fwrite($fnorm, $norm_sample . "\t" . $biorep . "\t" . $norm ."\n");
		}
	}
	fclose($fnorm);	
	echo "Normalization file {$normfile} written for all samples: ".join(",", array_keys($this->iqpir_params['samples']))."\n";
	
	
}
	$have_respairs = false;
	foreach($resPairMeans as $id => $ratios) {
		if(count($ratios) < 2) continue;
		if(! $have_respairs) {
			$fres = fopen($respairfile, "w");
			fwrite($fres, "protein{$prot_suff}\tresidue{$res_suff}\tmean_logratio\tstdev_logratio\tnumreps\t95conf\tpepA_seqs");
			if(! $this->deadend) fwrite($fres, "\tpepB_seqs");
			fwrite($fres, "\t{$file_hash}\t{$ratio_iontypes}\n");
			$have_respairs = true;
		}
		$nextstats = MyStatistics::stats($ratios, 2, 3);
		ksort($resPairPepAs[$id]);
		if(! $this->deadend) ksort($resPairPepBs[$id]);
		fwrite($fres, $id . "\t" . join("\t", $nextstats)."\t" . 
			number_format(1.96 * $nextstats[1] / sqrt($nextstats[2]), 3)."\t" .
			join(",", array_keys($resPairPepAs[$id])));
		if(! $this->deadend) fwrite($fres, "\t".join(",", array_keys($resPairPepBs[$id])));
		fwrite($fres, "\n");
	}
	if($have_respairs) {
		echo "Files {$samplefile}, {$chromfile}, and {$respairfile} written for sample {$sample}\n";
		fclose($fres);
	}
	else echo "Files {$samplefile} and {$chromfile} written for sample {$sample}\n"; //exit(1);
	$first = true;
	$prot_written = false;
	$remove_protpair_outliers = true;
	foreach($protPairMeans as $id => $ratios) {
		if(count($ratios) < 2) {
			if($add_zscore2ratiofile) unset($protPairMeans[$id]);
			continue;
		}
		if($first) {
			$fprot = fopen($protpairfile, "w");
			fwrite($fprot, "protein{$prot_suff}\tuniprot{$res_suff}\tmean_logratio\tstdev_logratio\tnumreps\t95conf\ttstat\tdf\tpvalue");
			if($remove_protpair_outliers) fwrite($fprot, "\toutliers");
			fwrite($fprot, "\t{$file_hash}\t{$ratio_iontypes}\n");
			$first = false;
			$prot_written = true;
		}
		if(! $remove_protpair_outliers) {
			$nextstats = MyStatistics::stats($ratios, 2, 3);
			fwrite($fprot, $id . "\t" . $protpair_unis[$id] . "\t" . join("\t", $nextstats)."\t" . 
				number_format(1.96 * $nextstats[1] / sqrt($nextstats[2]), 3)."\t".
				$this->getPvalue($nextstats[0], $nextstats[1], count($ratios)-1,true).
				"\n");
		}
		else {
			$nextstats = $this->getWeightedAverageRatio($ratios, NULL, true, 0, false, true, true);
			fwrite($fprot, $id . "\t" . $protpair_unis[$id] . "\t" . number_format($nextstats[0], 2, ".", "") . "\t" .
				number_format($nextstats[2], 3, ".", "") . "\t" . $nextstats[1] . "\t" . 
				number_format(1.96 * $nextstats[2] / sqrt($nextstats[1]), 3)."\t" .
				$nextstats[3] . "\t" . $nextstats[4] . "\n");
		}
		if($add_zscore2ratiofile) $protPairMeans[$id] = $nextstats;
	}
	if($prot_written) {
		echo "Files {$protpairfile} written for sample {$sample}\n";
		fclose($fprot);
	}
	if($add_zscore2ratiofile && file_exists($samplefile)) {
		$fout = fopen($samplefile . ".tmp", "w");
		$myfile = fopen($samplefile, "r");
		$first = true;
		while($line = fgets($myfile)) {
			$line = preg_replace( "/\r|\n/", "", $line);
			if($first) {
				fwrite($fout, $line . "\tprotpair_z\t{$file_hash}\t{$ratio_iontypes}\n");
				$first = false;
			}
			else {
				$next = split("\t", $line);
				$nextProts = split("-", $next[2]);
				$next_z = "";
				if($next[5] !=="") {
					if(array_key_exists($next[2], $protPairMeans)) {
						if($protPairMeans[$next[2]][1] > 0) $next_z = number_format(($next[5] - $protPairMeans[$next[2]][0])/$protPairMeans[$next[2]][1], 2, ".", "");
					}
					else if(! $this->deadend && array_key_exists($nextProts[1] . "-" . $nextProts[0], $protPairMeans)) {
						if($protPairMeans[$nextProts[1] . "-" . $nextProts[0]][1] > 0) 
						$next_z = number_format(($next[5] - $protPairMeans[$nextProts[1] . "-" . $nextProts[0]][0])/
						$protPairMeans[$nextProts[1] . "-" . $nextProts[0]][1], 2, ".", "");
					}
				}
				fwrite($fout, $line . "\t" . $next_z . "\n");
			}
		}
		fclose($fout);//exit(1);
		rename($samplefile . ".tmp", $samplefile);
	} // if add z score
	$first = true;
	uksort($resHubMeans, array($this, "sortByFirstAndSecondIndex"));
	$hub_written = false;
	foreach($resHubMeans as $id => $info) {
		if($first) {
			$fres = fopen($reshubfile, "w");
			fwrite($fres, "protein\tresidue\tmean_logratio\tstdev_logratio\tnumreps\t95conf\tpartner_protein\tpartner_residue\tmean_logratio\t{$file_hash}\t{$ratio_iontypes}\n");
			$first = false;
			$hub_written = true;
			continue;
		}
		usort($info, array($this, "sortByThirdIndex")); //exit(1);
		// get the stats
		$hubratios = array();
		for($k = 0; $k < count($info); $k++) array_push($hubratios, $info[$k][2]);
		$nextstats = MyStatistics::stats($hubratios, 2, 3);
		$nextinsert = join("\t", $nextstats)."\t" . 
			number_format(1.96 * $nextstats[1] / sqrt($nextstats[2]), 3);
		
		for($k = 0; $k < count($info); $k++) {
			fwrite($fres, $id . "\t" . $nextinsert . "\t" . join("\t", $info[$k])."\n");
		}
	}
	if($hub_written) {
		echo "File {$reshubfile} written for sample {$sample}\n";
		fclose($fres);
	}
	
	
	if($normalize && ! $this->use_current_normalization && count($normratios) >= 500) {
		$this->use_current_normalization = true;
		$this->computeTwoChannelLogratios($sample, $numerator_stumpmodmass, $denominator_stumpmodmass, $normalize);
		$this->use_current_normalization = false;
	}
		
	
	
}



private function sortByFirstAndSecondIndex($a, $b) {
	if(strcmp($a[0], $b[0]) < 0) return -1;
	if(strcmp($a[0], $b[0]) > 0) return 1;
	if($a[1] < $b[1]) return -1;
	if($a[1] > $b[1]) return -1;
	return 0;
}
private function sortByThirdIndex($a, $b) {
	if($a[2] < $b[2]) return -1;
	if($a[2] > $b[2]) return 1;
	return 0;
}

	
	public function getSampleArray() {
		$output = array();
		foreach($this->iqpir_params['samples'] as $key => $value) {
			$samplefile = $this->iqpir_params['iqpir_output_dir'] . "sample." . $key . ".ratios.txt";
			if(! file_exists($samplefile)) continue;
			$output[$key] = 1;
		}
		return $output;
	}


	public function getCompositeProteinZscores() {
		if($this->deadend) {
			echo "Error: cannot run getcomposite in deadend mode\n";
			exit(1);
		}
		foreach($this->iqpir_params['samples'] as $key => $value) {
			$samplefile = $this->iqpir_params['iqpir_output_dir'] . "sample." . $key . ".ratios.txt";
			$deadendprotpairfile = $this->iqpir_params['iqpir_output_dir'] . $this->deadend_subdir . "sample." . $key . ".protpairs.txt";
			$intraprotpairfile = $this->iqpir_params['iqpir_output_dir'] . $this->intra_only_subdir . "sample." . $key . ".protpairs.txt";
			$outfile = $this->iqpir_params['iqpir_output_dir'] . "sample." . $key . ".composite.txt";
			$deadendprots = array();
			$intraprots = array();
			if(! file_exists($samplefile)) {
				echo "Error: {$samplefile} does not exist\n";
				exit(1);
			}
			if(! file_exists($deadendprotpairfile)) {
				echo "Error: {$deadendprotpairfile} does not exist\n";
				exit(1);
			}
			if(! file_exists($intraprotpairfile)) {
				echo "Error: {$intraprotpairfile} does not exist\n";
				exit(1);
			}
			$myfile = fopen($deadendprotpairfile, "r");
			$first = true;
			while($line = rtrim(fgets($myfile))) {
				if($first) {
					$first = false;
					continue;
				}
				$next = split("\t", $line);
				if($next[4] < 3) continue; // num reps
				$deadends[$next[0]] = array($next[2], $next[3], $next[4]);
				
			}
			fclose($myfile);
			echo "Read in information for ".count($deadends)." from {$deadendprotpairfile}\n";
			$myfile = fopen($intraprotpairfile, "r");
			$first = true;
			while($line = rtrim(fgets($myfile))) {
				if($first) {
					$first = false;
					continue;
				}
				$next = split("\t", $line);
				if($next[4] < 3) continue; // num reps
				$nextprots = split("-", $next[0]);
				$intraprots[$nextprots[0]] = array($next[2], $next[3], $next[4]);
				
			}
			fclose($myfile);
			echo "Read in information for ".count($intraprots)." from {$intraprotpairfile}\n";
			$myfile = fopen($samplefile, "r");
			$fout = fopen($outfile, "w");
			$first = true;
			while($line = rtrim(fgets($myfile))) {
				if($first) {
					$first = false;
					fwrite($fout, $line . "\tdeadend_proA_z\tdeadend_proB_z\tintra_proA_z\tintra_proB_z\n");
					continue;
				}
				$next = split("\t", $line);
				if($next[6] < 1) continue; // num reps
				$nextprots = split("-", $next[2]);
				$deadend_zs = array("", "");
				$intra_zs = array("", "");
				for($k = 0; $k < 2; $k++) {
					if(array_key_exists($nextprots[$k], $deadends) && $deadends[$nextprots[$k]][1] > 0) $deadend_zs[$k] = number_format(($next[5] - $deadends[$nextprots[$k]][0])/$deadends[$nextprots[$k]][1], 2, ".", "");
					if(array_key_exists($nextprots[$k], $intraprots) && $intraprots[$nextprots[$k]][1] > 0) $intra_zs[$k] = number_format(($next[5] - $intraprots[$nextprots[$k]][0])/$intraprots[$nextprots[$k]][1], 2, ".", "");
				}
				fwrite($fout, $line . "\t" . join("\t", $deadend_zs)."\t" . join("\t", $intra_zs)."\n");
			}
			fclose($myfile);
			fclose($fout);
			echo "Composite information written to {$outfile}\n"; //exit(1);
		} // next sample
	
	}

	public function combinedDeadendIntralinkProteinRatios($samples = array()) {
		$myfile = fopen($this->iqpir_params['sample_ratios'], "r");
		$first = true;
		while($line=rtrim(fgets($myfile))){
			if($first) {
				$first = false;
			}
			else {
				$line = preg_replace("/\r|\n/", "", $line);
				$next = split("\t", $line);
				if(count($samples) > 0 && ! array_key_exists($next[0], $samples)) continue;
				$this->combinedTwoChannelDeadendIntralinkProteinRatios($next[0], $next[2], $next[4]); //exit(1);
			}
		}
		fclose($myfile);
	}
	
	public function combinedTwoChannelDeadendIntralinkProteinRatios($sample, $numerator_stumpmodmass, $denominator_stumpmodmass) {
		$suffix = ".{$numerator_stumpmodmass}-{$denominator_stumpmodmass}";
		$output_dir = $this->outputdir . $this->deadend_intra_subdir;
		$de_ratiofile = $this->outputdir . $this->deadend_subdir . "sample." . $sample . $suffix .".ratios.txt";
		$de_chromfile = $this->outputdir . $this->deadend_subdir . "sample." . $sample . $suffix .".chroms.txt";
		$intra_ratiofile = $this->outputdir . $this->intra_only_subdir . "sample." . $sample . $suffix .".ratios.txt";
		$intra_chromfile = $this->outputdir . $this->intra_only_subdir . "sample." . $sample . $suffix .".chroms.txt";
			if(! file_exists($de_ratiofile)) {
				echo "Error: ratiofile {$de_ratiofile} not found\n";
				echo "Please run LOGRAGTIOS DE\n";
				exit(1);
			}
			if(! file_exists($intra_ratiofile)) {
				echo "Error: ratiofile {$intra_ratiofile} not found\n";
				echo "Please run LOGRAGTIOS INTRA\n";
				exit(1);
			}
			$proteinratios = array(); // biorep to array
			$protpair_unis = array();
			$myfile = fopen($de_ratiofile, "r");
			$first = true;
			$file_hash_ion_types = "";
			while($line = rtrim(fgets($myfile))) {
				$next = split("\t", $line);
				if($first) {
					$first = false;
					if($file_hash_ion_types == "") $file_hash_ion_types = $next[15] . "\t" . $next[16];
					continue;
				}
				if($next[6] == 0) continue; // num reps
				if(! array_key_exists($next[2], $proteinratios)) $proteinratios[$next[2]] = array();
				array_push($proteinratios[$next[2]], $next[5]);
			}
			fclose($myfile);
			$myfile = fopen($intra_ratiofile, "r");
			$first = true;
			while($line = rtrim(fgets($myfile))) {
				if($first) {
					$first = false;
					continue;
				}
				$next = split("\t", $line);
				if($next[6] == 0) continue; // num reps
				$nextprots = split("-", $next[2]);
				if(! array_key_exists($nextprots[0], $proteinratios)) $proteinratios[$nextprots[0]] = array();
				array_push($proteinratios[$nextprots[0]], $next[5]);
			}
			fclose($myfile);
			$myfile = fopen($de_chromfile, "r");
			$first = true;
			while($line = rtrim(fgets($myfile))) {
				if($first) {
					$first = false;
					continue;
				}
				$next = split("\t", $line);
				$protpair_unis[$next[2]] = $next[3];
			}
			fclose($myfile);
			$myfile = fopen($intra_chromfile, "r");
			$first = true;
			while($line = rtrim(fgets($myfile))) {
				if($first) {
					$first = false;
					continue;
				}
				$next = split("\t", $line);
				$nextprots = split("-", $next[2]);
				$next_unis = split("-", $next[3]);

				$protpair_unis[$nextprots[0]] = $next_unis[0];
				$protpair_unis[$nextprots[1]] = $next_unis[1];
			}
			fclose($myfile);
			echo "Read in ".count($protpair_unis)." unis for ".count($proteinratios)." protein ratios for sample {$sample}.{$suffix}\n"; //exit(1);
			
			$add_zscore2ratiofile = false;
			$prot_written = false;
			$remove_protpair_outliers = true;
			$first = true;
			foreach($proteinratios as $id => $ratios) {
				if(count($ratios) < 2) {
					if($add_zscore2ratiofile) unset($protPairMeans[$id]);
					continue;
				}
				$protpairfile = $output_dir . "sample." . $sample . $suffix . ".protpairs.txt";
				if(! file_exists($output_dir)) system("mkdir ".$output_dir);
				if($first) {
					$fprot = fopen($protpairfile, "w");
					fwrite($fprot, "protein\tuniprot\tmean_logratio\tstdev_logratio\tnumreps\t95conf\ttstat\tdf\tpvalue");
					if($remove_protpair_outliers) fwrite($fprot, "\toutliers");
					fwrite($fprot, "{$file_hash_ion_types}\n");
					$first = false;
					$prot_written = true;
				}
				if(! $remove_protpair_outliers) {
					$nextstats = MyStatistics::stats($ratios, 2, 3);
					fwrite($fprot, $id . "\t" . $protpair_unis[$id] . "\t" . join("\t", $nextstats)."\t" . 
						number_format(1.96 * $nextstats[1] / sqrt($nextstats[2]), 3)."\t".
						$this->getPvalue($nextstats[0], $nextstats[1], count($ratios)-1,true).
						"\n");
				}
				else {
					$nextstats = $this->getWeightedAverageRatio($ratios, NULL, true, 0, false, true, true);
					fwrite($fprot, $id . "\t" . $protpair_unis[$id] . "\t" . number_format($nextstats[0], 2, ".", "") . "\t" .
						number_format($nextstats[2], 3, ".", "") . "\t" . $nextstats[1] . "\t".
						number_format(1.96 * $nextstats[2] / sqrt($nextstats[1]), 3).
						"\t" . $nextstats[3] . "\t" . $nextstats[4] . "\n");
				}
				if($add_zscore2ratiofile) $protPairMeans[$id] = $nextstats;
			}
			fclose($fprot);
			echo "Combined protein pair information written to file {$protpairfile}\n";
	
	}

	public function normalizeCrosslinksAndDeadends() {
		if(! $this->iqpir_params['normalize'] ) return;
		if($this->deadend || $this->intra_only) {
			echo "Error: cannot normalizeCrosslinksAndDeadends in deadend or intra_only mode\n";
			exit(1);
		}
		foreach($this->iqpir_params['samples'] as $sample => $value) {
			$chromfile = $this->iqpir_params['iqpir_output_dir'] . "sample." . $sample . ".chroms.txt";
			if(! file_exists($chromfile)) {
				echo "Error: chromfile {$chromfile} not found\n";
				exit(1);
			}
			$de_chromfile = $this->iqpir_params['iqpir_output_dir'] . $this->deadend_subdir . "sample." . $sample . ".chroms.txt";
			if(! file_exists($de_chromfile)) {
				echo "Error: chromfile {$de_chromfile} not found\n";
				exit(1);
			}
			$unnorm_logratios = array(); // biorep to array
			$myfile = fopen($chromfile, "r");
			$first = true;
			while($line = rtrim(fgets($myfile))) {
				if($first) {
					$first = false;
					continue;
				}
				$next = split("\t", $line);
				if($next[10] != "FALSE") continue; // num reps
				if(! array_key_exists($next[13], $unnorm_logratios)) $unnorm_logratios[$next[13]] = array();
				array_push($unnorm_logratios[$next[13]], $next[7]);
			}
			fclose($myfile);
			$myfile = fopen($de_chromfile, "r");
			$first = true;
			while($line = rtrim(fgets($myfile))) {
				if($first) {
					$first = false;
					continue;
				}
				$next = split("\t", $line);
				if($next[10] != "FALSE") continue; // num reps
				if(! array_key_exists($next[13], $unnorm_logratios)) $unnorm_logratios[$next[13]] = array();
				array_push($unnorm_logratios[$next[13]], $next[7]);
			}
			fclose($myfile);

        	foreach($unnorm_logratios as $biorep => $ratios) {
				$this->iqpir_params['samples'][$sample][$biorep] = array(number_format(-1 * MyStatistics::median($ratios), 2, ".", ""));
				array_push($this->iqpir_params['samples'][$sample][$biorep], $this->iqpir_params['samples'][$sample][$biorep][0]);
				echo "Have {$this->iqpir_params['samples'][$sample][$biorep][0]} norm for {$sample} sample and biorep {$biorep}\n";
				echo "Have optimal {$this->iqpir_params['samples'][$sample][$biorep][1]} norm for {$sample} sample and biorep {$biorep}\n";
			}

		} // next sammple
		if(! $this->use_current_normalization) {
			$normfile = $this->iqpir_params['iqpir_output_dir'] . "sample_biorep_normfactors.txt";
			$fnorm = fopen($normfile, "w");
			fwrite($fnorm, "sample\tbiorep\tnorm_add_to_log\n");
			foreach($this->iqpir_params['samples'] as $sample => $bioreps) {
				foreach($bioreps as $biorep => $norm) {
					fwrite($fnorm, $sample . "\t" . $biorep . "\t" . number_format($norm[0], 3, ".", "")."\n");
				}
			}
			fclose($fnorm);	
			echo "Normalization file {$normfile} written for all samples: ".join(",", array_keys($this->iqpir_params['samples']))."\n";
			if(file_exists($this->iqpir_params['iqpir_output_dir'] . $this->deadend_subdir . "sample_biorep_normfactors.txt")) system("cp ".$this->iqpir_params['iqpir_output_dir'] . "sample_biorep_normfactors.txt " . $this->iqpir_params['iqpir_output_dir'] . $this->deadend_subdir);
			if(file_exists($this->iqpir_params['iqpir_output_dir'] . $this->intra_only_subdir . "sample_biorep_normfactors.txt")) system("cp ".$this->iqpir_params['iqpir_output_dir'] . "sample_biorep_normfactors.txt " . $this->iqpir_params['iqpir_output_dir'] . $this->intra_only_subdir);
		}
	
	}

	public function compareDeadendAndCrosslinkSites() {
		$myfile = fopen($this->iqpir_params['sample_ratios'], "r");
		$first = true;
		$tally = array();
		while($line=rtrim(fgets($myfile))){
			if($first) {
				$first = false;
			}
			else {
				$line = preg_replace("/\r|\n/", "", $line);
				$next = split("\t", $line);
				if(count($samples) > 0 && ! array_key_exists($next[0], $samples)) continue;
				$sample = $next[0];
				$suffix = ".{$next[2]}-{$next[4]}";
				$de_ratiofile = $this->outputdir . $this->deadend_subdir . "sample." . $sample . $suffix .".ratios.txt";
				$intra_ratiofile = $this->outputdir . $this->intra_only_subdir . "sample." . $sample . $suffix .".ratios.txt";
				if(! file_exists($de_ratiofile)) {
					echo "Error: ratiofile {$de_ratiofile} not found\n";
					exit(1);
				}
				if(! file_exists($xl_ratiofile)) {
					echo "Error: ratiofile {$xl_ratiofile} not found\n";
					exit(1);
				}
				$deadendsites = array(); // biorep to array
				$myfile_de = fopen($de_ratiofile, "r");
				$first_de = true;
				while($line_de = rtrim(fgets($myfile_de))) {
					if($first_de) {
						$first_de = false;
						continue;
					}
					$next_de = split("\t", $line_de);
					if($next_de[6] == 0) continue; // num reps
					$deadendsites[$next_de[2] . "-" . $next_de[3]] = 1; // protein-res
				}
				fclose($myfile_de);
				$myfile_xl = fopen($xl_ratiofile, "r");
				$first_xl = true;
				$tot = 0;
				$tally[$sample] = array("tot_des" => count($deadendsites), "tot_xls" => 0, "xl_0" => 0, "xl_1" => 0, "xl_2" => 0, "xl_2hd" => 0); 
				while($line_xl = rtrim(fgets($myfile_xl))) {
					if($first_xl) {
						$first_xl = false;
						continue;
					}
					$next_xl = split("\t", $line);
					if($next_xl[6] == 0) continue; // num reps
					$prots = split("-", $next_xl[2]);
					$sites = split("-", $next_xl[3]);
					$num = 0;
					if(array_key_exists($prots[0] . "-" . $sites[0], $deadendsites)) $num++;
					if(array_key_exists($prots[1] . "-" . $sites[1], $deadendsites)) $num++;
					$is_homodimer = $num == 2 && $prots[0]==$prots[1] && $sites[0] == $sites[1] ? "hd" : "";
					$tally[$sample]["xl_" . $num . $is_homodimer]++;
					$tally[$sample]["tot_xls"]++;
				}
				fclose($myfile_xl);
				if($tally[$sample]["tot_xls"] > 0) {
					$tally[$sample]["xl_0"] = number_format($tally[$sample]["xl_0"] / $tally[$sample]["tot_xls"], 3, ".", "");
					$tally[$sample]["xl_1"] = number_format($tally[$sample]["xl_1"] / $tally[$sample]["tot_xls"], 3, ".", "");
					$tally[$sample]["xl_2"] = number_format($tally[$sample]["xl_2"] / $tally[$sample]["tot_xls"], 3, ".", "");
					$tally[$sample]["xl_2hd"] = number_format($tally[$sample]["xl_2hd"] / $tally[$sample]["tot_xls"], 3, ".", "");
				}
			} // next sample
			$first_tally = true;
			ksort($tally);
			foreach($tally as $sample => $stats) {
				if($first_tally) {
					echo "sample\t".join("\t", array_keys($stats))."\n";
					$first_tally = false;
				}
				echo $sample . "\t" . join("\t", array_values($stats))."\n";
			}
		} // next sample and ratio num denom
		fclose($myfile);
	}
	// will have to adapt this to sample.iqpir.txt file format from iqpir
	private function reuseIqpirfiles($files, $sample) {
		$output = array(); // hashed by spectrum scan pepA proA modposA pepB proB modposB
		if(count($files) == 0) return $output;
		for($k = 0; $k < count($files); $k++) {
			$myfile = fopen($files[$k], "r");
			$first = true;
			while($line = rtrim(fgets($myfile))) {
				if($first) {
					$first = false;
				}
				else {
					$next = split("\t", $line);
					if(! array_key_exists($next[5], $this->iqpir_params['rawfiles']) ||
						$this->iqpir_params['rawfiles'][$next[5]]['sample'] != $sample) continue;
					if($this->iqpir_params['rawfiles'][$next[5]]['biorep'] != $next[13]) {
						$next[13] = $this->iqpir_params['rawfiles'][$next[5]]['biorep'];
						$line = join("\t", $next); // set up with new biorep
					}
					$next_xl = split("_", $next[1]);
					$next_peps = split("_", $next[0]);
					$next_prots = split("-", $next[3]);
					if($this->deadend) {
						$next_xl = array_merge($next_xl, $next_xl);
						$next_peps = array_merge($next_peps, $next_peps);
						$next_prots = array_merge($next_prots, $next_prots);
					}
					$next_scan = $next[6];
					while(strlen($next_scan) < 5) $next_scan = "0" . $next_scan;
					
					$next_unique = $next[5] . "." . $next_scan . "." . $next_scan . "\t" .
					$next_peps[0] . "\t" . $next_prots[0] . "\t" . $next_xl[1] . "\t" . $next_peps[1] . "\t" .
					$next_prots[1] . "\t" . $next_xl[3];
					$next_id = $next[7];
					if(! array_key_exists($next_unique, $output)) $output[$next_unique] = array();
					if(! array_key_exists($next_id, $output[$next_unique])) $output[$next_unique][$next_id] = $line;
				}
			}
			fclose($myfile);
		}
		return $output;

	}


}


?>
