<?php

include_once("MyPepXML.php");
include_once("MyDeadendPepXML.php");
include_once("MyStatistics.php"); // use this library for the insert statement

class MyIqPir
{

	private $iqpir_paramsfile = NULL;
	private $scan_ratiofile = NULL;
	private $scan_ratiofile_lh = NULL;
	private $iqpir_params = NULL; // store the parameters and sample map info
	private $iqpir_ppmtolerance = 20; //20; //20; //25;
	private $ms2_scan = ""; // so can distinghish ms3 from ms2
	private $ms3_ppmtol_factor = 20;
	private $quantitation_type = "iqpir";
	private $filter_based_on_theorrelints_lt = 0.4; //0.4; //0.3; //0; // 0.2 required relative intensity agreement with theor
	private $filter_based_on_theorrelints_hv = 0.4; //0.4; //0.3; //0; // 0.2 required relative intensity agreement with theor
	private $use_optimal_ratios =  true; // 
	private $max_optimal_error = 0.1; //1; //0.1; //0.1;
	private $subtract_noise_from_peaks = false; // whether to reduce observed readmzxml intensities by avearge
	private $min_sig2noise = 5; // 5
	private $memory_limit = '4000M'; //'2000M' //'1000M';
	private $deadend_subdir = "deadend/";
	private $deadend = false;
	private $intra_only_subdir = "intra_only/";
	private $intra_only = false; // whether to restrict logratios to intra-protein cross-links only
	private $use_current_normalization = false;
	private $deadend_intra_subdir = "deadend_intra/";
	private $peaklist = array();
	private $normalization_method = "median"; // default, or mode


	public function __construct ($iqpir_paramsfile) { //$samplemap, $ml_filename, $workingdir, $masschroq_outputfile, $quantitation_type="silac") {
		if(! file_exists($iqpir_paramsfile)) {
			echo "Error: iqpir params file {$iqpir_paramsfile} does not exist\n";
			exit(1);
		}
		$this->iqpir_paramsfile = $iqpir_paramsfile;
		$this->readIqPirParams();
		$this->scan_ratiofile = $this->getScanratioFilename();
		$this->scan_ratiofile_lh = $this->getScanratioFilename(false, true);
		ini_set('memory_limit', $this->memory_limit);	

		if(array_key_exists("iqpir_ppmtolerance", $this->iqpir_params) && $this->iqpir_params["iqpir_ppmtolerance"] != "")
			$this->iqpir_ppmtolerance = $this->iqpir_params["iqpir_ppmtolerance"];
		if(array_key_exists("ms3_ppmtol_factor", $this->iqpir_params) && $this->iqpir_params["ms3_ppmtol_factor"] != "")
			$this->iqpir_ppmtolerance = $this->iqpir_params["ms3_ppmtol_factor"];
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
			$next = explode("\t", $line);
			$this->iqpir_params['samples'][$next[0]][$next[1]] = array($next[2], $next[2]);
		}
		fclose($myfile);
		echo "sample\tbiorep\tnorm_add_to_log\n";
		foreach($this->iqpir_params['samples'] as $sample => $bioreps) {
			foreach($bioreps as $biorep => $norm) {
				echo $sample . "\t" . $biorep . "\t" . number_format($norm[0], 3, ".", "")."\n";
			}
		}
	}
	
	public function setDeadendMode() {
		if($this->deadend) return; // already in deadend mode
		$this->iqpir_params['iqpir_output_dir'] .= $this->deadend_subdir;
		if($this->iqpir_params['iqpir_output_dir'][strlen($this->iqpir_params['iqpir_output_dir'])-1]!="/") $this->iqpir_params['iqpir_output_dir'] .= "/";
		if(! file_exists($this->iqpir_params['iqpir_output_dir'])) system("mkdir {$this->iqpir_params['iqpir_output_dir']}");
		$this->scan_ratiofile = $this->getScanratioFilename(true);
		$this->scan_ratiofile_lh = $this->getScanratioFilename(true, true);
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

public function generateSingeRatioIqpirFile() {
	$file_suff = $this->deadend ? $this->deadend_subdir : "";

	$output_dir = $this->iqpir_params['iqpir_output_dir'];
	if($output_dir[strlen($output_dir)-1]!= "/") $output_dir .= "/";
	$output_dir .= "combined_onesample/" . $file_suff;
	if(! file_exists($output_dir)) system("mkdir {$output_dir}");
	$current_dir = getcwd() . "/";
	chdir($output_dir);
	system("ln -s {$this->scan_ratiofile} .");
	chdir($current_dir);
	$files = glob($output_dir."/*xls");

	if(count($files) == 1) echo "One sample {$files[0]} written ready for LOGRATIOS step with {$this->iqpir_params['onesample_iqpirparamsfile']}\n";
	else echo "Check to see if one sample iqpir file was sybolically linked to {$output_dir}\n";
}

public function generateSingleRatioParamFiles() {
	if(! array_key_exists("onesample_iqpirparamsfile", $this->iqpir_params) || $this->iqpir_params['onesample_iqpirparamsfile'] == "") {
		echo "Error: you must specify onesample_iqpirparamsfile parameter in your params file {$this->iqpir_paramsfile}\n";
		exit(1);
	}
	if(stripos($this->iqpir_params['onesample_iqpirparamsfile'], "iqpir")===false) {
		echo "Error: onesample_iqpirparamsfile must have iqpir in its name to be recognized, {$this->iqpir_params['onesample_iqpirparamsfile']} does not\n";
		echo "Please change your parameter value accordingly\n";
		exit(1);
	}
	$dir = "";
	$dir = getcwd() . "/";
	$nextpos = strrpos($this->iqpir_paramsfile, "/");
	if($nextpos !== false) {
		if($this->iqpir_paramsfile[0] == "/")
			$dir = substr($this->iqpir_paramsfile, 0, $nextpos) . "/";
		else $dir .= substr($this->iqpir_paramsfile, 0, $nextpos) . "/";
	}
	
	$onesample_iqpirparamsfile = $this->iqpir_params['onesample_iqpirparamsfile'];
	if($onesample_iqpirparamsfile[0] != "/") $onesample_iqpirparamsfile = $dir . $onesample_iqpirparamsfile;
	if(substr($onesample_iqpirparamsfile, strlen($onesample_iqpirparamsfile)-7) != ".params") $onesample_iqpirparamsfile .= ".params";
	if(file_exists($onesample_iqpirparamsfile)) return; // nothing to do
	$onesample_sample_map = $dir . "onesample_map.txt";
	$onesample_sample_def = $dir . "onesample_def.txt";
	$output_dir = $this->iqpir_params['iqpir_output_dir'];
	if($output_dir[strlen($output_dir)-1] != "/") $output_dir .= "/";
	$output_dir .= "combined_onesample/";
	
	$exp_cond = "combined_";
	$ref_cond = "";
	$descr = "";
	if(file_exists($dir . "sample_def.txt")) {
		$valid = array();
		exec("grep ^1 ".$dir . "sample_def.txt", $valid);
		for($k = 0; $k < count($valid); $k++) {
			$valid[$k] = rtrim($valid[$k]);
			$next = explode("\t", $valid[$k]);
			if($next[0] == 1) {
				$exp_cond .= $next[1];
				$ref_cond = $next[2];
				$descr = $next[3];
				$k = count($valid);
			}
		}
	}
	// write the new sample_map
	$fmap = fopen($onesample_sample_map, "w");
	$biorep = 0;
	$current = "";
	foreach($this->iqpir_params['rawfiles'] as $rawfile => $info) { // sample, biorep
		$next = $info['sample'] . "-" . $info['biorep'];
		if($next != $current) {
			$biorep++;
			$current = $next;
		}
		fwrite($fmap, $rawfile . "\t1\t" . $biorep . "\t" . $info['orientation'] . "\n");
	
	}
	fclose($fmap);
	echo "One sample sample_map file {$onesample_sample_map} written\n";


	$fdef = fopen($onesample_sample_def, "w");
	fwrite($fdef, "sample\texp_cond\tref_cond\tdescription\n");
	fwrite($fdef, "1\t{$exp_cond}\t{$ref_cond}\t{$descr}\n");
	fclose($fdef);
	echo "One sample sample_def file {$onesample_sample_def} written\n";
	// write it
	$fparams = fopen($onesample_iqpirparamsfile, "w");
	fwrite($fparams, "sample_map\t{$onesample_sample_map}\n");
	fwrite($fparams, "email\t{$this->iqpir_params['email']}\n");
	fwrite($fparams, "xlinkprophetfile\t{$this->iqpir_params['xlinkprophetfile']}\n");
	if(array_key_exists("deadendprophetfile", $this->iqpir_params)) fwrite($fparams, "deadendprophetfile\t{$this->iqpir_params['deadendprophetfile']}\n");
	fwrite($fparams, "fdr\t{$this->iqpir_params['fdr']}\n");
	fwrite($fparams, "filter_crit\t{$this->iqpir_params['filter_crit']}\n");
	fwrite($fparams, "iqpir_output_dir\t{$output_dir}\n");
	fwrite($fparams, "normalize\t{$this->iqpir_params['normalize']}\n");

	if(array_key_exists('normalization_method', $this->iqpir_params))  
		fwrite($fparams, "normalization_method\t{$this->iqpir_params['normalization_method']}\n");

	fclose($fparams);
	echo "One sample iqpir file {$onesample_iqpirparamsfile} written with output directory {$output_dir}\n";
}

private function getStumpFragments($peptide, $stump_mod, $min_mass = 500, $max_mass = 2500) {
	$aa_masses = array("A" => 71.03711, "R" => 156.10111, "N" => 114.04293, "D" => 115.02694, "C" => 103.00919, "E" => 129.04259, "Q" => 128.05858,
						"G" => 57.02146, "H" => 137.05891, "I" => 113.08406, "L" => 113.08406, "K" => 128.09496, "M" => 131.04049, "F" => 147.06841,
						"P" => 97.05276, "S" => 87.03203, "T" => 101.04768, "W" => 186.07931, "Y" => 163.06333, "V" => 99.06841);
						
	// go through the peptide counting residues for both b and y ions
	$b_index = 0;
	$stump_res = -1;
	$hydrogen_mass = 1.007825035;
	$carbon_mass = 12.000000;
	$oxygen_mass = 15.99491463;
	$proton_mass = 1.00727647;
	$nitrogen_mass = 14.003074;
	$modmasses = array("325.13" => 325.127386, "147.04" => 147.035385);
	if(array_key_exists(strval($stump_mod), $modmasses)) $stump_mod = $modmasses[strval($stump_mod)];
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
			}
			else {
				echo "Error with no ] in {$peptide}\n";
				exit(1);
			}
		}
		if($b_mass < $min_mass || $b_mass > $max_mass) continue;
		if($stump_res > -1 && $b_index > $stump_res) {
			$output['b' . $b_index] = array($b_mass, substr($peptide, 0, $k + 1));
			$last = 'b' . $b_index;
		}
	
	}
	if($last != "") unset($output[$last]); // remove the full length b ion
	$stripped = MyPepXML::stripMods($peptide); // now do the y ions
	$y_index = 0;
	for($k = strlen($stripped)-1; $k > 0; $k--) {
		$y_index++;
		$y_mass += $aa_masses[$stripped[$k]];
		if(array_key_exists($k, $mods)) $y_mass += $mods[$k];
		if($y_mass < $min_mass || $y_mass > $max_mass) continue;
		if($k < $stump_res) $output['y' . $y_index] = array($y_mass, substr($peptide,  $index_first[strlen($stripped) + 1 - $y_index]));
	}
	return $output;
}

public function getOneSampleParamsfile() {
	if(array_key_exists('onesample_iqpirparamsfile', $this->iqpir_params)) {
		$onesample_iqpirparamsfile = $this->iqpir_params['onesample_iqpirparamsfile'];
		if(substr($onesample_iqpirparamsfile, strlen($onesample_iqpirparamsfile)-7) != ".params") $onesample_iqpirparamsfile .= ".params";
		return $this->iqpir_params['onesample_iqpirparamsfile'];
	}
	return "";
}

private function readIqPirParams() { // iqpir.params and reference sample_map
$this->iqpir_params = array();
$myfile = fopen($this->iqpir_paramsfile, "r");
while($line=rtrim(fgets($myfile))){
	$next = explode("\t", $line);
	if(count($next) > 1 && $next[1] != "") $this->iqpir_params[$next[0]] = $next[1];
}
fclose($myfile);
if(! array_key_exists('sample_map', $this->iqpir_params) || ! file_exists($this->iqpir_params['sample_map'])) {
	echo "Error with missing sample_map in {$this->iqpir_paramsfile} or directory {$this->iqpir_params['sample_map']}\n";
	if(array_key_exists('email', $this->iqpir_params) && $this->iqpir_params['email'] != "")
		mail($this->iqpir_params['email'], "Error with missing sample_map in {$this->iqpir_paramsfile} or directory {$this->iqpir_params['sample_map']}", "Error with missing sample_map in {$this->iqpir_paramsfile} or directory {$this->iqpir_params['sample_map']}\r\n",
			"From: brucelabgenomesciences@gmail.com\r\n");
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
	$next = explode("\t", $line);
	$this->iqpir_params['rawfiles'][$next[0]] = array("orientation" => $next[1], 
		"sample" => $next[2], "biorep" => $next[3]);
	if(! array_key_exists($next[2], $this->iqpir_params['samples'])) $this->iqpir_params['samples'][$next[2]] = array();
	$this->iqpir_params['samples'][$next[2]][$next[3]] = array(0, 0);
}
fclose($myfile);
ksort($this->iqpir_params['samples']);
if(array_key_exists('normalization_method', $this->iqpir_params) && $this->iqpir_params['normalization_method'] == "mode")  {
	echo "Setting normalization method to ".$this->iqpir_params['normalization_method']."\n";
	$this->normalization_method = $this->iqpir_params['normalization_method'];
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

public function getOptimalRatio($theor_relints, $obs_relints) {
if($obs_relints[0] == 0 || $obs_relints[1] == 0 || $obs_relints[2] == 0 || $obs_relints[3] == 0) return "";

if($this->filter_based_on_theorrelints_lt > 0 && $obs_relints[0] > 0) {
	$next = $obs_relints[1]/$obs_relints[0];
	$next_theor = $theor_relints[1]/$theor_relints[0];
	if(abs($next - $next_theor) > $this->filter_based_on_theorrelints_lt) return "";
}
if($this->filter_based_on_theorrelints_hv > 0 && $obs_relints[3] > 0) {
	$next = $obs_relints[4]/$obs_relints[3];
	$next_theor = $theor_relints[2]/$theor_relints[1];
	if(abs($next - $next_theor) > $this->filter_based_on_theorrelints_hv) return "";
}

	$s = array();
	$m = array();
	$z = $theor_relints[0]+$theor_relints[1]+$theor_relints[2];
	for($k = 0; $k < 5; $k++) {
		array_push($s, $theor_relints[$k] - $obs_relints[$k]);
		if($k > 1) array_push($m, $obs_relints[$k]*$z - $theor_relints[$k-2]);
		else array_push($m, $obs_relints[$k]*$z);
	}
	$a = 0;
	$b = 0;
	$c = 0;
	for($k = 0; $k < 5; $k++) {
		$a += $m[$k] * $m[$k];
		$b += $s[$k] * $m[$k];
		$c += $s[$k] * $s[$k];
	}
	if($z == 0 || $c+$b == 0) {
		return "";
	}
	$ratio = ($a + $b*$z)/($z*$c + $b);
	if($ratio < 0) return 0;
	if(! $this->isPosSecondErrorDeriv($a, $b, $c, $z, $ratio)) {
		echo "Negative second derivative error for {$ratio}\t{$a}\t{$b}\t{$c}\t{$z}\n"; 
		return "";
	}
	if($this->max_optimal_error > 0 && 
		$this->computeRatioDeconvError($theor_relints, $obs_relints, $ratio) > $this->max_optimal_error) return "";
	return number_format($ratio, 2, ".", "");
}
private function isPosSecondErrorDeriv($a, $b, $c, $z, $ratio) {
	return $ratio * (-4 * $z * $c - 4 * $b) + 2 * $z * $z * $c + 8 * $z * $b + 6 * $a > 0;
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
			$current = explode('|', $current);
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
echo "Read in ".count($output)." protein sequences from {$fasta}\n"; 
return $output; 
}


private function getIsotopeQuant($mass, $charge, $scan, $mzxml, $noise = 0, $title = "", $num_isotopes = 3) {
	$proton_mass = 1.00727647;
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
	$ppm_diff = $nextmz * $ppmtolerance / 1000000;
	$done = false;
	foreach($this->peaklist[$scan] as $next_mz => $next_int) {
		if($done) continue;
		foreach($pepmzs as $offset => $mz) {
			for($j = 0; $j < count($pepmzs); $j++) {
				if($next_mz >= $pepmzs[$j] - $ppm_diff && $next_mz <= $pepmzs[$j] + $ppm_diff) {
					if($this->subtract_noise_from_peaks) $next_int = max(0, $next_int-$noise);
					$pepints[$j] += $next_int;
					if($verbose) echo "adding {$next_int} for {$pepmzs[$j]}\n";
					if($next_int > $max_mzs[$j][1]) {
						$max_mzs[$j][0] = $next_mz;
						$max_mzs[$j][1] = $next_int;
					}
				}
				else if($j==count($pepmzs)-1 && $next_mz > $pepmzs[$j] + $ppm_diff) { // done
					$done = true;
				}
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
	return array($pepints, $sig2noise, $pepmzs, $num_mz_diffs > 0 ? number_format($ave_mz_diff/$num_mz_diffs, 2, ".", "") : "N/A");

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

private function getPeptideHeavy2LightRatioWithTheorInts($light, $heavy, $ints, $heavy2light = false) {

if(! is_array($light) || ! is_array($heavy) || ! is_array($light[0]) || ! is_array($heavy[0])) return "N/A";
if(count($ints) < 2) return "";
// check that first 2 lights are approximately correct ratio
if($light[0][0] == 0 || $ints[0] == 0) return "";
$next = $light[0][1]/$light[0][0];
$next_theor = $ints[1]/$ints[0];
if($this->filter_based_on_theorrelints_lt > 0 && abs($next - $next_theor) > $this->filter_based_on_theorrelints_lt) return "";

if($heavy[0][1] == 0 || $ints[1] == 0) return "";
$next = $heavy[0][2]/$heavy[0][1];
$next_theor = $ints[2]/$ints[1];
if($this->filter_based_on_theorrelints_hv > 0 && abs($next - $next_theor) > $this->filter_based_on_theorrelints_hv) return "";
$final_light = ($light[0][0] + $light[0][1]) / ($ints[0] + $ints[1]);
$final_heavy = (max(0, ($light[0][2] - $ints[2] * $final_light))+ 
				max(0, ($heavy[0][1] - $ints[3] * $final_light)) +
				$heavy[0][2])/($ints[0] + $ints[1] + $ints[2]);


if(! $heavy2light) return $final_heavy > 0 ? number_format($final_light / $final_heavy, 2) : "";
	
return $final_light > 0 ? number_format($final_heavy / $final_light, 2) : "";

}


private function getPeptideHeavy2LightRatio($light, $heavy, $heavy2light = false) {
if(! is_array($light) || ! is_array($heavy) || ! is_array($light[0]) || ! is_array($heavy[0])) return "N/A";
	$next_light3 = max(0, ($light[0][1] > 0 ? ($light[0][2] - $light[0][0] * $heavy[0][1]) / $light[0][1] : floatval($light[0][0])));
	$next_heavy1 = max(0, ($light[0][1] > 0 ? $light[0][0] * $heavy[0][1] / $light[0][1] : floatval($heavy[0][1])));
	$light = $light[0][0] + $light[0][1] + $next_light3;
	$heavy = $next_heavy1 + $heavy[0][1] + $heavy[0][2];
	return $light > 0 ? number_format($heavy / $light, 2, ".", '') : "";
}

// input are the dimers of intensities and sig2noise of the 3 isotopes
private function calculateIntensErrorVersusTheor($obs, $theor, $num) {
$error = 0;
if(count($obs) < $num) return 100;
$tot_obs = 0;
$tot_th = 0;
for($k = 0; $k < $num; $k++) {
	$tot_obs += $obs[$k];
	$tot_th += $theor[$k];
}
if($tot_obs == 0) return 1000;
for($k = 0; $k < $num; $k++) {
	$error += ($obs[$k]/$tot_obs - $theor[$k]/$tot_th) * ($obs[$k]/$tot_obs - $theor[$k]/$tot_th);
}
return number_format(sqrt($error)/$num, 3);

}

private function calculateError($ints, $posdistr, $negdistr, $prior) {
// first get the total intensity
$totint = 0;
for($k = 0; $k < count($ints);  $k++) {
	$totint += $ints[$k];
}
$totcorr = $totint * $prior;
$totincorr = $totint - $totcorr;
$error = 0;
for($k = 0; $k < count($ints);  $k++) {
	$expected_corr = $posdistr[$k] * $totcorr;
	$expected_incorr = $negdistr[$k] * $totincorr;
	$error += ($ints[$k] - $expected_corr - $expected_incorr) * ($ints[$k] - $expected_corr - $expected_incorr);
}
return number_format(sqrt($error) / $totint, 2, ".", "");
}

// returns light intensity, heavy intensity, and error
private function computeRatioFromIsotopeSetsOfPeptide($light_isotope_ints,  $mass_or_peptide, $heavy_light_offset = 2) {
	// must compute the vector based on theoretical intensities
	if(count($light_isotope_ints) < 5) {
		echo "Error: need at least 5 isotope peaks\n";
		exit(1);
	}
	$ints = $this->getIsotopeIntensities($mass_or_peptide, true);
	while(count($ints) < 5) {
		array_push($ints, 0);
	}
	$A = array();
	$nextprod = $ints[0]*$ints[0]+$ints[1]*$ints[1]+$ints[2]*$ints[2];
	$nextxprod = $ints[0]*$ints[2]+$ints[1]*$ints[3]+$ints[2]*$ints[4];
	$next_row1 = array($nextprod + $ints[3]*$ints[3]+$ints[4]*$ints[4], $nextxprod);
	$next_row2 = array($nextxprod, $nextprod);
	array_push($A, $next_row1);
	array_push($A, $next_row2);
	$x = array($ints[0]*$light_isotope_ints[0]+$ints[1]*$light_isotope_ints[1]+$ints[2]*$light_isotope_ints[2]+$ints[3]*$light_isotope_ints[3]+$ints[4]*$light_isotope_ints[4],
		$ints[0]*$light_isotope_ints[2]+$ints[1]*$light_isotope_ints[3]+$ints[2]*$light_isotope_ints[4]);
	$solution = MyStatistics::gaussian_solve($A, $x);
	if($solution[0] < 0) $solution[0] = 0;
	if($solution[1] < 0) $solution[1] = 0;
	// now compute the error
	$error = ($ints[0]*$solution[0] - $light_isotope_ints[0])*($ints[0]*$solution[0] - $light_isotope_ints[0]) +
				($ints[1]*$solution[0] - $light_isotope_ints[1])*($ints[1]*$solution[0] - $light_isotope_ints[1]) +
				($ints[2]*$solution[0]+$ints[0]*$solution[1] - $light_isotope_ints[2])*($ints[2]*$solution[0]+$ints[0]*$solution[1] - $light_isotope_ints[2]) +
				($ints[3]*$solution[0]+$ints[1]*$solution[1] - $light_isotope_ints[3])*($ints[3]*$solution[0]+$ints[1]*$solution[1] - $light_isotope_ints[3]) +
				($ints[4]*$solution[0]+$ints[2]*$solution[1] - $light_isotope_ints[4])*($ints[4]*$solution[0]+$ints[2]*$solution[1] - $light_isotope_ints[4]);
	$error = sqrt($error);
	if($solution[0]>0 || $solution[1] > 0) {
		$error /= ($solution[0]+$solution[1]);
	}
	$output = array(number_format($solution[0], 0, ",", "") . "-" . number_format($solution[1], 0, ",", ""), 
		$solution[0] < 1 && $solution[1] < 1 ? "undef" : $solution[0] == 0 ? "inf" : number_format($solution[1]/$solution[0], 2), number_format($error, 1));
	return $output;
	array_push($solution, $error);
	return $solution;
}


// arrays of light and heavy intensities, the offset (number of isotope peaks that heavy lags behind light), and the sequence or mass of the peptide
private function computeRatioFromIsotopeSets($light_isotope_ints, $heavy_isotope_ints, $heavy_light_offset, $mass_or_peptide) {
/	$ints = $this->getIsotopeIntensities($mass_or_peptide, true);
	while(count($ints) < $heavy_light_offset+1) {
		array_push($ints, 0);
	}
	$next_ratio1 = "N/A";
	if($light_isotope_ints[0] > 0 && $heavy_isotope_ints[0] > 0 && $ints[0] >= 0 && $ints[$heavy_light_offset] >= 0) {
		$next_ratio1 = $light_isotope_ints[0] / ($heavy_isotope_ints[0] - ($ints[$heavy_light_offset] * $light_isotope_ints[0]/$ints[0]));
		if($next_ratio1 < 0) {
			$next_ratio1 = "N/A";
		}
		else $next_ratio1 = number_format($next_ratio1, 2);
	}
	$next_ratio2 = "N/A";
	if($light_isotope_ints[1] > 0 && $heavy_isotope_ints[1] > 0 && $ints[1] >= 0 && $ints[$heavy_light_offset+1] >= 0) {
		$next_ratio2 = $light_isotope_ints[1] / ($heavy_isotope_ints[1] - ($ints[$heavy_light_offset+1] * $light_isotope_ints[1]/$ints[1]));
		if($next_ratio2 < 0) {
			$next_ratio2 = "N/A";
		}
		else $next_ratio2 = number_format($next_ratio2, 2);
	}
	return array($next_ratio1, $next_ratio2, $ints);
}


public function getIsotopeIntensities($seq_or_mass, $normalize = false, $num_isotopes = 5) {
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
	$next = preg_split('/\s+/', $valid[$k]); //explode("\t", $valid[$k]);
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

private function chargeStateDeconvolute($mzXML, $ppm_tol, $h_l_massdiff, $myscan = NULL) {
$hk2_file = $mzXML;
$nextpos = strpos($hk2_file, ".mzXML");
if($nextpos!==false) $hk2_file = substr($hk2_file, 0, $nextpos) . ".hk2";
if(! file_exists($hk2_file)) {
	echo "Error: {$hk2_file} does not exist\n";
	exit(1);
}
$myfile = fopen($hk2_file, "r");
$scan = 0;
$spectrum = array();
$output = array();
$ints = array();
while($line=rtrim(fgets($myfile))){
	$next = explode("\t", $line);
	if($next[0]=="S") {
		if(count($spectrum) > 0) {
			$masses = array();
			usort($spectrum, array($this, "sortBySecondThenThirdIndex"));
			$verbose = false; // $scan == "10310";
			for($k = 0; $k < count($spectrum); $k++) {
				if($verbose) echo "{$scan}: Next mass {$spectrum[$k][1]} with m/z {$spectrum[$k][4]}\n";
				$masses[$spectrum[$k][1]] = $spectrum[$k][2];
				$ints[$spectrum[$k][1]] = array($spectrum[$k][3]);
				$index = $k;
				while($index < count($spectrum)-1 && $spectrum[$index+1][2] != $spectrum[$index][2] &&
					($spectrum[$index+1][1] - $spectrum[$index][1]) / $spectrum[$index+1][1] * 1000000 <= $ppm_tol) {
					$masses[$spectrum[$k][1]] .= "," . $spectrum[$index+1][2];
					array_push($ints[$spectrum[$k][1]], $spectrum[$index+1][3]);
					$index++; // skip it, already have that charge
				}
				$k = $index;
			}
			$mass_vals = array_keys($masses);

			for($k = 0; $k < count($mass_vals) - 3; $k++) {
				// check for partner ahead
				$theor_ints = $this->getIsotopeIntensities($mass_vals[$k], true);
				$exp_rel_int = $theor_ints[1]/$theor_ints[0];
				$obs_rel_int = $ints[$mass_vals[$k]][0] > 0 ? $ints[$mass_vals[$k+1]][0]/$ints[$mass_vals[$k]][0] : 999;
			
				$j = $k + 1;
					if(abs($mass_vals[$j] - $mass_vals[$k] - 1.00335) / $mass_vals[$k] * 1000000 <= $ppm_tol &&
					
					
						abs($mass_vals[$j+1] - $mass_vals[$k] - $h_l_massdiff) / $mass_vals[$k] * 1000000 <= $ppm_tol &&
						abs($mass_vals[$j+2] - $mass_vals[$k] - $h_l_massdiff - 1.00335) / $mass_vals[$k] * 1000000 <= $ppm_tol &&
						abs($obs_rel_int - $exp_rel_int) < 0.2
						// add on expected relative intensities based on mass
						
						
						) {
						// can check intensities here.....
						if(is_null($myscan)) {
							if(! array_key_exists($scan, $output)) $output[$scan] = array();
							$output[$scan][$mass_vals[$k]] = array($masses[$mass_vals[$k]], $mass_vals[$j], $masses[$mass_vals[$j]], $ints[$mass_vals[$k]]);
						}
						else {
							$output[$mass_vals[$k]] = array($masses[$mass_vals[$k]], $mass_vals[$j], $masses[$mass_vals[$j]], $ints[$mass_vals[$k]]);
						}
						$k = $j + 2;
						$j = count($mass_vals);
					}
			}
			if(! is_null($myscan)) {
				fclose($myfile);
				return $output;
			}
		}
		$scan = $next[1];
		$spectrum = array();
	}
	else if($next[0]=="P") {
		if(is_null($myscan) || $myscan == $scan)
			array_push($spectrum, $next);
	}
}
fclose($myfile);
if(count($spectrum) > 0) {
	$masses = array();
	usort($spectrum, array($this, "sortBySecondThenThirdIndex"));
	for($k = 0; $k < count($spectrum); $k++) {
		$masses[$spectrum[$k][1]] = $spectrum[$k][2];
		$ints[$spectrum[$k][1]] = array($spectrum[$k][3]);
		$index = $k;
		while($index < count($spectrum)-1 && $spectrum[$index+1][2] != $spectrum[$index][2] &&
			($spectrum[$index+1][1] - $spectrum[$index][1]) / $spectrum[$index+1][1] * 1000000 <= $ppm_tol) {
			$masses[$spectrum[$k][1]] .= "," . $spectrum[$index+1][2];
			array_push($ints[$spectrum[$k][1]], $spectrum[$index+1][3]);
			$index++; // skip it, already have that charge
		}
		$k = $index;
	}
	ksort($masses);
	// now look for differences
	$mass_vals = array_keys($masses);
	$theor_ints = $this->getIsotopeIntensities($mass_vals[$k], true);
	$exp_rel_int = $theor_ints[1]/$theor_ints[0];
	$obs_rel_int = $ints[$mass_vals[$k]][0] > 0 ? $ints[$mass_vals[$k+1]][0]/$ints[$mass_vals[$k]][0] : 999;
	for($k = 0; $k < count($mass_vals) - 3; $k++) {
		// check for partner ahead
		$j = $k + 1;
		
		if(abs($mass_vals[$j] - $mass_vals[$k] - 1.00335) / $mass_vals[$k] * 1000000 <= $ppm_tol &&
					
					
			abs($mass_vals[$j+1] - $mass_vals[$k] - $h_l_massdiff) / $mass_vals[$k] * 1000000 <= $ppm_tol &&
			abs($mass_vals[$j+2] - $mass_vals[$k] - $h_l_massdiff - 1.00335) / $mass_vals[$k] * 1000000 <= $ppm_tol &&
			abs($obs_rel_int - $exp_rel_int) < 0.2
			// add on expected relative intensities based on mass
						
						
			) {

				if(is_null($myscan)) {
					if(! array_key_exists($scan, $output)) $output[$scan] = array();
					$output[$scan][$mass_vals[$k]] = array($masses[$mass_vals[$k]], $mass_vals[$j], $masses[$mass_vals[$j]], $ints[$mass_vals[$k]]);
				}
				else {
					$output[$mass_vals[$k]] = array($masses[$mass_vals[$k]], $mass_vals[$j], $masses[$mass_vals[$j]], $ints[$mass_vals[$k]]);
				}
				$j = count($mass_vals);
			}
	}
}
return $output;
}

private function isHeavy($modpep) {
	return strpos($modpep, "K[327")!==false;
}

private function convertHeavyToLight($modpep) {
	$output = preg_replace("/K\[327.13\]/", "K[325.13]", $modpep);
	$output = preg_replace("/K\[327.14\]/", "K[325.13]", $output);
	if($this->deadend)
		$output = preg_replace("/K\[327\]/", "K[325]", $output);
	return $output;
}


function getScanratioFilename($is_deadend = false, $light_heavy = false) {
	$file = $is_deadend && array_key_exists("deadendprophetfile", $this->iqpir_params) ?
		$this->iqpir_params['deadendprophetfile'] : $this->iqpir_params['xlinkprophetfile'];
	if(! file_exists($file)) {
		echo "Error: xlinkprophet pepxml file {$file} does not exist\n";
		exit(1);
	}
	$output_dir = $this->iqpir_params['iqpir_output_dir'];
	$file_suffix = "-iqpir.xls"; 
	$outfile = $file; 
	
	$nextpos = strrpos($file, "/");
	if($nextpos!==false) {
		$outfile = $output_dir . substr($file, $nextpos+1);
	}
	$nextpos = strpos($outfile, ".pep.xml");
	if($nextpos!==false) $outfile = substr($outfile, 0, $nextpos);
	$outfile .= $file_suffix;
	if($light_heavy) return substr($outfile, 0, strlen($outfile)-4)."_l-h.txt";
	return $outfile;
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

public function reuseIqpirfiles($files) {
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
			$next = explode("\t", $line);
			$next_scan = $next[0];
			$next_spec = substr($next[1], 0, strlen($next[1])-2);
			$nextpos = strrpos($next_spec, $next_scan);
			if($nextpos === strlen($next_spec) - strlen($next_scan)) { //
			// must check whether scan is in spectrum (Mango) or not (ReACT)
				$next_unique = $next_spec . "\t" . $this->convertHeavyToLight($next[3]) . "\t" .
					$next[4] . "\t" . $next[5] . "\t" . $this->convertHeavyToLight($next[6]) . "\t" . 
					$next[7] . "\t" . $next[8];
			}
			else { // ReACT... add scan to base name
				$nextpos = strpos($next_spec, ".");
				if($nextpos !== false) {
					$next_unique = substr($next_spec, 0, $nextpos) . "." . $next_scan . "\t" . $this->convertHeavyToLight($next[3]) . "\t" .
						$next[4] . "\t" . $next[5] . "\t" . $this->convertHeavyToLight($next[6]) . "\t" . 
						$next[7] . "\t" . $next[8];
				}
				else {
					echo "Error: cannot parse spectrum {$next_spec} in {$line}\n";
					exit(1);
				}
			}
			if(array_key_exists($next_unique, $output)) {
			}
			else {
				$output[$next_unique] = $line;
			}
		}
	}
	fclose($myfile);
}
return $output;

}

public function readPepXML($spec_filter = "", $deconvolute = false) {

	if($this->deadend) return $this->quantifyDeadends($spec_filter, $deconvolute);

function log2ratio($val, $log2bound = "") { //10) {
	if($val === "") return "";
	else if($val === 0) return $log2bound === "" ? $log2bound : -1 * $log2bound;
	else if($val == "inf") return $log2bound;
	
	return number_format(log($val)/log(2), 2, ".", "");;
}

	$file = $this->iqpir_params['xlinkprophetfile'];
	$fdr = $this->iqpir_params['fdr'];
	$level = $this->iqpir_params['filter_crit'];
	$heavy_mod = $this->iqpir_params['heavy_modmass'];
	$output_dir = $this->iqpir_params['iqpir_output_dir'];

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
		$spec_filter = explode(",", $spec_filter);
	}
	else $spec_filter = array();
	
	$output_index = 8;
	$file_suffix = "-iqpir.xls"; //$deconvolute ? ".quantd_test{$output_index}{$filter_suff}.xls" : ".quant_test{$output_index}{$filter_suff}.xls";
	$find_residuals = false; // wheterh to look for peaks with expected h/l spacing
	$pursue_all_pep_charges = true;
	
	$report_log2ratios = true;
	$protein_seqs_and_genes = array();

	$ratio_prefix = $report_log2ratios ? "log2" : "";
	$outfile = $this->scan_ratiofile; //getOutfile($file, $output_dir);
	echo "Ready to write to {$outfile} with fdr {$fdr}, filter level {$level}, and heavy mod {$heavy_mod}\n"; 
	$fout = fopen($outfile, "w");
	if($deconvolute) {
		fwrite($fout, "scan\tspectrum\tpeptide1\tprotein1\tpos1\tpeptide2\tprotein2\tpos2\tnoise\tisoints1\tlight1-1\tlight1-2\tlight1-3\theavy1-1\theavy1-2\theavy1-3\tdIntens1\t{$pepRatio}ratio1\t{$pepRatio}ratio1Error");
		fwrite($fout, "\tisoints2\tlight2-1\tlight2-2\tlight2-3\theavy2-1\theavy2-2\theavy2-3\tdIintens2\t{$pepRatio}ratio2\t{$pepRatio}ratio2Error");
		fwrite($fout, "\tisointsreporter\tlightreporter-1\tlightreporter-2\tlightreporter-3\theavyreporter-1\theavyreporter-2\theavyreporter-3\t{$reporterRatio}ratioreporter-1\t{$reporterRatio}ratioreporter-2");
		fwrite($fout, "\tPepA_Int1\tPepA_Int2\tPepB_Int1\tPepB_Int2\n");
	}
	else {
	
	
		fwrite($fout, "scan\tspectrum\tactual_ratio\tpeptide1\tprotein1\tpos1\tpeptide2\tprotein2\tpos2\tnoise\tlight1-1\tlight1-2\tlight1-3\theavy1-1\theavy1-2\theavy1-3\t{$pepRatio}ratio1-1\t{$pepRatio}ratio1-2\t{$pepRatio}ratio1-3\tdeconv_h-l_ratio");

		if($find_residuals || $pursue_all_pep_charges) {
			for($k = 1; $k <= 3; $k++) {
				if($find_residuals) {
					fwrite($fout, "\tresiduals\tresidual_sum");
				}
				if($pursue_all_pep_charges) {
					fwrite($fout, "\tdeconv_h-l_ratio\tl-h_mz\tints{$k}_1-1\tintratio{$k}_1-1\tlightratio{$k}_1-1\tlight{$k}_1-1\tlight{$k}_1-2\tlight{$k}_1-3\theavy{$k}_1-1\theavy{$k}_1-2\theavy{$k}_1-3\t{$pepRatio}ratio{$k}_1-1\t{$pepRatio}ratio{$k}_1-2\t{$pepRatio}ratio{$k}_1-3");
				}

			}
		}
		fwrite($fout, "\tlight2-1\tlight2-2\tlight2-3\theavy2-1\theavy2-2\theavy2-3\t{$pepRatio}ratio2-1\t{$pepRatio}ratio2-2\t{$pepRatio}ratio2-3\tdeconv_h-l_ratio");
		if($find_residuals || $pursue_all_pep_charges) {
			for($k = 1; $k <= 3; $k++) {
				if($find_residuals) {
					fwrite($fout, "\tresiduals\tresidual_sum");
				}
				if($pursue_all_pep_charges) {
					fwrite($fout, "\tdeconv_h-l_ratio\tl-h_mz\tints{$k}_2-1\tintratio{$k}_2-1\tlightratio{$k}_2-1\tlight{$k}_2-1\tlight{$k}_2-2\tlight{$k}_2-3\theavy{$k}_2-1\theavy{$k}_2-2\theavy{$k}_2-3\t{$pepRatio}ratio{$k}_2-1\t{$pepRatio}ratio{$k}_2-2\t{$pepRatio}ratio{$k}_2-3");
				}
			}
		}
		fwrite($fout, "\tlightreporter-1\tlightreporter-2\tlightreporter-3\theavyreporter-1\theavyreporter-2\theavyreporter-3\t{$reporterRatio}ratioreporter-1\t{$reporterRatio}ratioreporter-2\t{$reporterRatio}ratioreporter-3");
	}
	fwrite($fout, "\tactual_{$ratio_prefix}ratio");
	fwrite($fout, "\treporter_ids\treporter_log2ratios\treporter_ints\tpeptide_ids\tpeptide_log2ratios\tpeptide_ints");
	fwrite($fout, "\tpeptide_mzs\tpeptide_theor_relints\tpetpide_obs_relints");
	fwrite($fout, "\tfrag_ids\tfrag_log2ratios\tfrag_ints");
	fwrite($fout, "\tfrag_mzs\tfrag_theor_relints\tfrag_obs_relints\tfrag_scans\tfrag_noise");
	fwrite($fout, "\tlongarm_reporter_ids\tlongarm_reporter_log2ratios\tlongarm_reporter_ints\tlongarm_scans");
if(! $this->use_optimal_ratios) {			
			fwrite($fout, "\toptimal_peptide_log2ratios\toptimal_frag_log2ratios");
}	
	fwrite($fout, "\tresidue_pair");
	fwrite($fout, "\tprotein_pair");
	fwrite($fout, "\tnoise\n");
	
	$existing_iqpirquant = array();
	if(array_key_exists("existing_xl_iqpirfiles", $this->iqpir_params) && $this->iqpir_params["existing_xl_iqpirfiles"] != "") {
		$existing_iqpirfiles = explode(",", $this->iqpir_params["existing_xl_iqpirfiles"]); //"/net/gs/vol4/shared/brucelab/search/xiaoting/Spectrast_test/interact_spectrast-xl-iqpir.xls", "/net/gs/vol4/shared/brucelab/search/xiaoting/iqPIR_mito_Tac_sham_mango_6pairs_135runs/iprophet-110920-xl-iqpir.xls");
		for($k = 0; $k < count($existing_iqpirfiles); $k++) {
			if(strpos($existing_iqpirfiles[$k], "-xl-iqpir.xls") !== strlen($existing_iqpirfiles[$k]) - strlen("-xl-iqpir.xls")) {
				echo "Error: existing iqpirfile {$existing_iqpirfiles[$k]} is not of correct format ending in -xl-iqpir.xls\n";
				exit(1);
			}
			else if(! file_exists($existing_iqpirfiles[$k])) {
				echo "Error: existing iqpirfile {$existing_iqpirfiles[$k]} is not found\n";
				exit(1);
			}
		}
		$existing_iqpirquant = $this->reuseIqpirfiles($existing_iqpirfiles);
		echo "Read in existing quantitation for ".count($existing_iqpirquant)." spectral cross-link scans\n"; 
	}
	
	$pepxml = new MyPepXML($file, $fdr, $level);
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
	$basename = "";
	while($pepxml->nextLine()) {
		if($pepxml->hasStartTag("msms_run_summary")) { // end of searh result
			$mzxml = $pepxml->getTagValue("base_name") . ".mzXML";
			$noisefile = $pepxml->getTagValue("base_name") . ".noise";
			$basename = $pepxml->getTagValue("base_name");
			$nextpos = strrpos($basename, "/");
			if($nextpos !== false) $basename = substr($basename, $nextpos + 1);
			if(! file_exists($mzxml)) {
				echo "Error: mzxml {$mzxml} does not exist\n";
				exit(1);
			}
			if(! file_exists($noisefile)) {
				echo "Error: noisefile {$noisefile} does not exist\n";
				exit(1);
			}
		}
		else if($pepxml->hasStartTag("search_summary")) { // end of searh result
			$analysis_type = $pepxml->getAnalysisType();
		}
		else if($pepxml->hasStartTag("spectrum_query")) { // end of searh result
			$spectrum = $pepxml->getTagValue("spectrum");
		
			if($analysis_type == "Mango" || $analysis_type == "SpectraST") {
				$scan = $pepxml->getTagValue("start_scan");
			}
		}
		else if($analysis_type == "ReACT" && $pepxml->hasStartTag("search_hit")) { // end of searh result
			$scan = $pepxml->getTagValue("ms2scan");
		}
		else if($pepxml->hasStartTag("linked_peptide")) { // end of searh result
			array_push($masses, $pepxml->getTagValue("calc_neutral_pep_mass"));
			array_push($charges, $pepxml->getTagValue("assumed_charge"));
			if($analysis_type == "ReACT") { // end of searh result
				array_push($ms3scans, $pepxml->getTagValue("ms3_scan"));
			}
		}
		else if(count($protein_seqs_and_genes)==0 && $pepxml->hasStartTag("search_database")) { // end of searh result
			$protein_seqs_and_genes = $this->getDatabaseProteinSequencesAndGenes($pepxml->getTagValue("local_path"));
		}
		else if($pepxml->hasStartTag("modification_info")) { // end of searh result
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
			$next_unique = ($analysis_type == "ReACT" ? $basename . "." .  $scan : $spec_without_ch) . "\t" . $this->convertHeavyToLight($pepxml->getPreservedOrderedModifiedPeptidePair());
			if($pepxml->aboveFilterThreshold() && array_key_exists($next_unique, $seen)) {
				$heavy = strpos($peps[0], $heavy_mod)!==false;
				$next_index = $heavy ? 1 : 0;
				$seen[$next_unique][$next_index]++; 
			}
			else if($pepxml->aboveFilterThreshold() && ! array_key_exists($next_unique, $seen)) {
			
			
				$seen[$next_unique] = array(0, 0);
				$noise = 0;
				$noise_valid = array();
				exec("grep -P '^".$scan."\t' {$noisefile}", $noise_valid);
				if(count($noise_valid)===1) {
					$next_noise = explode("\t", $noise_valid[0]);
					$noise = $next_noise[3];
				}
				$heavy = strpos($peps[0], $heavy_mod)!==false;
				$next_index = $heavy ? 1 : 0;
				$seen[$next_unique][$next_index]++; 

				if(count($existing_iqpirquant) > 0 && array_key_exists($next_unique, $existing_iqpirquant)) {
					fwrite($fout, $existing_iqpirquant[$next_unique] . "\n");
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
								$modpeps2[0] = $this->convertHeavyToLight($modpeps2[0]); //preg_replace("/K\[327.13\]/", "K[325.13]", $modpeps[0]);
								$modpeps2[1] = $this->convertHeavyToLight($modpeps2[1]); //preg_replace("/K\[327.13\]/", "K[325.13]", $modpeps[1]);
							}
							
				$pep1_stumpfrags = $this->getStumpFragments($modpeps2[0], 325.13); 
				$pep2_stumpfrags = $modpeps2[0] == $modpeps2[1] ? 
					array() : $this->getStumpFragments($modpeps2[1], 325.13); 
				$frag_ratios = array();
				$frag_ints = array();
				$seen_fragmasses = array(); // make sure only use once
				$fragment_relints = array();
				$pep1_frag_ratios = array();
				$pep1_frag_ints = array();
				$pep2_frag_ratios = array();
				$pep2_frag_ints = array();
				$frag_fragments = array();
				$frag_scans = array();
				$frag_mzs = array();
				$react_ms3 = true;
				$next_scans = $analysis_type == "ReACT" && $react_ms3 ? $ms3scans : array($scan, $scan);
				$frag_noise = $analysis_type == "ReACT" && $react_ms3 ? 0 : $noise; // either equal to noise for MS2 or 0 for MS3

				$this->ms2_scan = $scan;
				$longarm_scans = array_key_exists('longarm_reporters', $this->iqpir_params) && 
					$this->iqpir_params['longarm_reporters'] == "true" && $analysis_type == "ReACT" ?
					array($ms3scans[0] + 2, $ms3scans[1] + 2) : array();
					
				$this->setScanPeaklist($mzxml, $ms3scans, $longarm_scans);
				$fragment_theor_ints = array();
				$peptide_theor_ints = array();
				$peptide_relints = array();
				
				$peptide_pair = $pepxml->getPreservedOrderedModifiedPeptidePair(false, true);
				
				$samepeps = $peptide_pair[0]==$peptide_pair[3];
				
				$min_peptide_massdiff = 0.5; // daltons?
				$indistinct_pepmasses = abs($masses[0]-$masses[1]) < $min_peptide_massdiff;
				foreach($pep1_stumpfrags as $key => $value) {	
					$next_mass = strval($value[0]);
					$nextpos = strpos($next_mass, ".");
					if($nextpos!==false) $next_mass = substr($next_mass, 0, $nextpos+4); // keep only first 2 digits
					if(array_key_exists($next_mass, $seen_fragmasses)) continue;
					$seen_fragmasses[$next_mass] = 1;
					$light_fragentries[$key] = $this->getIsotopeQuant($value[0], 1, $next_scans[0], $mzxml, $noise, "frag_light{$key}_1 ");
					$heavy_fragentries[$key] = $this->getIsotopeQuant($value[0] + $pep_massdiff, 1, $next_scans[0], $mzxml, $noise, "frag_heavy{$key}_1 ");
					// now compute the ratios
					$next_ints_fragentries = $this->getIsotopeIntensities($value[1], true);
					if($next_ints_fragentries[0]==0) {
						echo "1. No theoreteical intensities obtained with calcisotopes for fragment ion \"{$value[1]}\": ".join(",", $next_ints_fragentries).", possibly due to unexpected peptide modification (contact Jimmy).  Will instead pass mass {$value[0]}\n";
						$next_ints_fragentries = $this->getIsotopeIntensities($value[0], true);
					}
					if($this->use_optimal_ratios) {
						$next_fragdeconv = "N/A";
						$next_fragintens = $light_fragentries[$key][0][0] + $light_fragentries[$key][0][1] + $light_fragentries[$key][0][2] + $heavy_fragentries[$key][0][1] + $heavy_fragentries[$key][0][2];
						if($next_fragintens > 0) {
							$next_fragrelints = array(number_format($light_fragentries[$key][0][0]/$next_fragintens, 3, ".", ""), number_format($light_fragentries[$key][0][1]/$next_fragintens, 3, ".", ""), number_format($light_fragentries[$key][0][2]/$next_fragintens, 3, ".", ""), 
								number_format($heavy_fragentries[$key][0][1]/$next_fragintens, 3, ".", ""), number_format($heavy_fragentries[$key][0][2]/$next_fragintens, 3, ".", ""));
						
							$next_fragdeconv = $this->getOptimalRatio($next_ints_fragentries, $next_fragrelints);
						}
					}
					else {
						$next_fragdeconv = $this->getPeptideHeavy2LightRatioWithTheorInts($light_fragentries[$key], $heavy_fragentries[$key], $next_ints_fragentries);
					}
					if($next_fragdeconv != "N/A" && $next_fragdeconv !== "") {
							$next_fragintens = $light_fragentries[$key][0][0] + $light_fragentries[$key][0][1] + $light_fragentries[$key][0][2] + $heavy_fragentries[$key][0][1] + $heavy_fragentries[$key][0][2];
							if($next_fragintens > 0) {
								array_push($frag_ints, $next_fragintens);
								array_push($frag_ratios, $next_fragdeconv);
								array_push($pep1_frag_ints, $next_fragintens);
								array_push($pep1_frag_ratios, $next_fragdeconv);
								if($samepeps && $next_scans[0] == $next_scans[1]) array_push($frag_fragments, "AB:" . $key . "_1");
								else array_push($frag_fragments, "A:" . $key . "_1");
							
								array_push($fragment_theor_ints, join(",", $next_ints_fragentries));
		
								$next_fragrelints = array(number_format($light_fragentries[$key][0][0]/$next_fragintens, 3, ".", ""), number_format($light_fragentries[$key][0][1]/$next_fragintens, 3, ".", ""), number_format($light_fragentries[$key][0][2]/$next_fragintens, 3, ".", ""), 
									number_format($heavy_fragentries[$key][0][1]/$next_fragintens, 3, ".", ""), number_format($heavy_fragentries[$key][0][2]/$next_fragintens, 3, ".", ""));
								array_push($fragment_relints, join(",", $next_fragrelints));
								array_push($frag_scans, $next_scans[0]);
								array_push($frag_mzs, ($value[0] + $proton_mass) . "-" . ($value[0] + $pep_massdiff + $proton_mass));
							}
					}

				}
				foreach($pep2_stumpfrags as $key => $value) {		
					$next_mass = strval($value[0]);
					$nextpos = strpos($next_mass, ".");
					if($nextpos!==false) $next_mass = substr($next_mass, 0, $nextpos+4);
					if(array_key_exists($next_mass, $seen_fragmasses)) continue;
					$seen_fragmasses[$next_mass] = 1;
					$light_fragentries[$key] = $this->getIsotopeQuant($value[0], 1, $next_scans[1], $mzxml, $noise, "frag_light{$key}_1 ");
					$heavy_fragentries[$key] = $this->getIsotopeQuant($value[0] + $pep_massdiff, 1, $next_scans[1], $mzxml, $noise, "frag_heavy{$key}_1 ");
					// now compute the ratios
					$next_ints_fragentries = $this->getIsotopeIntensities($value[1], true);
					if($next_ints_fragentries[0]==0)  {
						echo "2. No theoreteical intensities obtained with calcisotopes for fragment ion \"{$value[1]}\": ".join(",", $next_ints_fragentries).", possibly due to unexpected peptide modification (contact Jimmy).  Will instead pass mass {$value[0]}\n";
						$next_ints_fragentries = $this->getIsotopeIntensities($value[0], true);
					}
					if($this->use_optimal_ratios) {
						$next_fragdeconv = "N/A";
						$next_fragintens = $light_fragentries[$key][0][0] + $light_fragentries[$key][0][1] + $light_fragentries[$key][0][2] + $heavy_fragentries[$key][0][1] + $heavy_fragentries[$key][0][2];
						if($next_fragintens > 0) {
							$next_fragrelints = array(number_format($light_fragentries[$key][0][0]/$next_fragintens, 3, ".", ""), number_format($light_fragentries[$key][0][1]/$next_fragintens, 3, ".", ""), number_format($light_fragentries[$key][0][2]/$next_fragintens, 3, ".", ""), 
								number_format($heavy_fragentries[$key][0][1]/$next_fragintens, 3, ".", ""), number_format($heavy_fragentries[$key][0][2]/$next_fragintens, 3, ".", ""));
						
							$next_fragdeconv = $this->getOptimalRatio($next_ints_fragentries, $next_fragrelints);
						}
					}
					else {
						$next_fragdeconv = $this->getPeptideHeavy2LightRatioWithTheorInts($light_fragentries[$key], $heavy_fragentries[$key], $next_ints_fragentries);
					}
					if(! $samepeps && $next_fragdeconv != "N/A" && $next_fragdeconv !== "" && $next_fragdeconv != "inf" && $next_fragdeconv > 0) {  // 011020
							$next_fragintens = $light_fragentries[$key][0][0] + $light_fragentries[$key][0][1] + $light_fragentries[$key][0][2] + $heavy_fragentries[$key][0][1] + $heavy_fragentries[$key][0][2];
							if($next_fragintens > 0) {
								array_push($frag_ints, $next_fragintens);
								array_push($frag_ratios, $next_fragdeconv);
								array_push($pep2_frag_ints, $next_fragintens);
								array_push($pep2_frag_ratios, $next_fragdeconv);
								array_push($frag_fragments, "B:" . $key . "_1");
								array_push($fragment_theor_ints, join(",", $next_ints_fragentries));
								$next_fragrelints = array(number_format($light_fragentries[$key][0][0]/$next_fragintens, 3, ".", ""), number_format($light_fragentries[$key][0][1]/$next_fragintens, 3, ".", ""), number_format($light_fragentries[$key][0][2]/$next_fragintens, 3, ".", ""), 
									number_format($heavy_fragentries[$key][0][1]/$next_fragintens, 3, ".", ""), number_format($heavy_fragentries[$key][0][2]/$next_fragintens, 3, ".", ""));
								array_push($fragment_relints, join(",", $next_fragrelints));
								array_push($frag_scans, $next_scans[1]);
								array_push($frag_mzs, ($value[0] + $proton_mass) . "-" . ($value[0] + $pep_massdiff + $proton_mass));
							}
					}



				}
				if(count($frag_ratios) > 0) {
					$nextFragInfo = $this->getWeightAverageWithRatio($frag_ratios, $frag_ints, NULL, NULL, $report_log2ratios);
				}
				
				$light_entries1 = $pursue_all_pep_charges ? array("1" => array(), "2" => array(), "3" => array()) : array(); // charges 1,2,3
				$heavy_entries1 = array("1" => array(), "2" => array(), "3" => array()); // charges 1,2,3
				$dlight1_entries1 = array("1" => array(), "2" => array(), "3" => array()); // charges 1,2,3
				$ratios_entries1 = array("1" => array(), "2" => array(), "3" => array()); // charges 1,2,3
				$ints_entries1 = array("1" => array(), "2" => array(), "3" => array()); // charges 1,2,3
				$deconv_ratio_entries1 = array("1" => array(), "2" => array(), "3" => array()); // charges 1,2,3
				$mzs1 = array();  // from charge to array of light/heavy
				$ppmtolerance = $this->iqpir_ppmtolerance; //25;
				$next_hk_info = $find_residuals ? $this->chargeStateDeconvolute($mzxml, $ppmtolerance, $pep_massdiff, $scan) : array();
				// take only ones that are smaller than greater of both peptide masses
				foreach($next_hk_info as $key => $value) {
					if($key > $masses[0] && $key > $masses[1]) unset($next_hk_info[$key]);
					else if(abs($key - $masses[0]) <= 0.2) unset($next_hk_info[$key]);
					else if(abs($key - $masses[1]) <= 0.2) unset($next_hk_info[$key]);
				}				
				foreach($light_entries1 as $key => $value) {
					$next_charge = intval($key);
					$symbol = $next_charge == $charges[0] ? "*" : "";
					if($heavy) {
						if($masses[0]/$next_charge >= $min_mz && $masses[0]/$next_charge <= $max_mz) {
							$light_entries1[$key] = $this->getIsotopeQuant($masses[0] - $pep_massdiff, $next_charge, $scan, $mzxml, $noise, "{$symbol}light{$key}_1 ");
							$heavy_entries1[$key] = $this->getIsotopeQuant($masses[0], $next_charge, $scan, $mzxml, $noise, "{$symbol}heavy{$key}_1 ");
							if($deconvolute) $dlight1_entries1[$key] = $this->getIsotopeQuant($masses[0] - $pep_massdiff, $next_charge, $scan, $mzxml, $noise, "light{$key}_1 ", $num_isotopes);
							
							$mzs1[$key] = array(($masses[0] - $pep_massdiff + $next_charge * $proton_mass)/$next_charge, ($masses[0] + $next_charge * $proton_mass)/$next_charge);
						}
					}
					else {
						if(($masses[0]+$pep_massdiff)/$next_charge >= $min_mz && ($masses[0]+$pep_massdiff)/$next_charge <= $max_mz) {
							$light_entries1[$key] = $this->getIsotopeQuant($masses[0], $next_charge, $scan, $mzxml, $noise, "light{$key}_1 ");
							$heavy_entries1[$key] = $this->getIsotopeQuant($masses[0] + $pep_massdiff , $next_charge, $scan, $mzxml, $noise, "{$symbol}heavy{$key}_1 ");
							if($deconvolute) $dlight1_entries1[$key] = $this->getIsotopeQuant($masses[0], $next_charge, $scan, $mzxml, $noise, "{$symbol}light{$key}_1 ", $num_isotopes);
							$mzs1[$key] = array(($masses[0] + $next_charge * $proton_mass)/$next_charge, ($masses[0] + $pep_massdiff + $next_charge * $proton_mass)/$next_charge);
						}
					}
					if(count($light_entries1[$key]) > 0 && count($heavy_entries1[$key]) > 0) {
						$ratios_entries1[$key] = $this->printRatios($light_entries1[$key][0], $heavy_entries1[$key][0], "Scan {$scan} Pep1 {$key}_Ratios: ");

						if($deconvolute) {
							$modpeps = $pepxml->getModifiedPeptidePair();
							if($heavy) { // must change sequence to light version
								$modpeps[0] = $this->convertHeavyToLight($modpeps[0]); //preg_replace("/K\[327.13\]/", "K[325.13]", $modpeps[0]);
								$modpeps[1] = $this->convertHeavyToLight($modpeps[1]); //preg_replace("/K\[327.13\]/", "K[325.13]", $modpeps[1]);
							}

							$deconv_ratio_entries1[$key] = $this->computeRatioFromIsotopeSetsOfPeptide($dlight1_entries1[$key][0],  $modpeps[0]);

						}
						$ints1 = array();
						if($deconvolute) {
							$ratios_entries1[$key] = $this->computeRatioFromIsotopeSets($light_entries1[$key][0], $heavy_entries1[$key][0], 2, $modpeps[0]);
							$ints_entries1[$key] = array_pop($ratios_entries1[$key]);
							
						}
						else {
							$modpeps = $pepxml->getModifiedPeptidePair();
							if($heavy) { // must change sequence to light version
								$modpeps[0] = $this->convertHeavyToLight($modpeps[0]); //preg_replace("/K\[327.13\]/", "K[325.13]", $modpeps[0]);
								$modpeps[1] = $this->convertHeavyToLight($modpeps[1]); //preg_replace("/K\[327.13\]/", "K[325.13]", $modpeps[1]);
							}
							$ints_entries1[$key] = $this->getIsotopeIntensities($modpeps[0], true);
							if($ints_entries1[$key][0]==0) {
								echo "3. No theoreteical intensities obtained with calcisotopes for peptide ion \"{$modpeps[0]}\": ".join(",", $ints_entries1[$key]).", possibly due to unexpected peptide modification (contact Jimmy).  Will instead pass mass {$masses[0]}\n";
								echo "here with {$key} {$peps[0]} and {$peps[1]} and heavy ? ".($heavy ? "yes" : "no")." with peps ".join(",", $peps)."\n";
								$ints_entries1[$key] = $this->getIsotopeIntensities($masses[0], true);
							}
						}
					}
				}
				
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
				
					if($deconvolute) {
						$modpeps = $pepxml->getModifiedPeptidePair();
						if($heavy) { // must change sequence to light version
							$modpeps[0] = $this->convertHeavyToLight($modpeps[0]); //preg_replace("/K\[327.13\]/", "K[325.13]", $modpeps[0]);
							$modpeps[1] = $this->convertHeavyToLight($modpeps[1]); //preg_replace("/K\[327.13\]/", "K[325.13]", $modpeps[1]);
						}
					}
				
				$ints1 = array();
				if($deconvolute) {
					$ratios1 = $this->computeRatioFromIsotopeSets($light1ints[0], $heavy1ints[0], 2, $modpeps[0]);
					$ints1 = array_pop($ratios1);
				}

				$dlight2 = NULL;
				
				$light_entries2 = array();
				$heavy_entries2 = array();
				if(! $samepeps) {
					$light_entries2 = $pursue_all_pep_charges ? array("1" => array(), "2" => array(), "3" => array()) : array(); // charges 1,2,3
					$heavy_entries2 = array("1" => array(), "2" => array(), "3" => array()); // charges 1,2,3
				}
				$dlight1_entries2 = array("1" => array(), "2" => array(), "3" => array()); // charges 1,2,3
				$ratios_entries2 = array("1" => array(), "2" => array(), "3" => array()); // charges 1,2,3
				$ints_entries2 = array("1" => array(), "2" => array(), "3" => array()); // charges 1,2,3
				$deconv_ratio_entries2 = array("1" => array(), "2" => array(), "3" => array()); // charges 1,2,3
				$mzs2 = array();  // from charge to array of light/heavy
				foreach($light_entries2 as $key => $value) {
					$next_charge = intval($key);
					$symbol = $next_charge == $charges[1] ? "*" : "";
					if($heavy) {
						if($masses[1]/$next_charge >= $min_mz && $masses[1]/$next_charge <= $max_mz) {
							$light_entries2[$key] = $this->getIsotopeQuant($masses[1] - $pep_massdiff, $next_charge, $scan, $mzxml, $noise, "{$symbol}light{$key}_1 ");
							$heavy_entries2[$key] = $this->getIsotopeQuant($masses[1], $next_charge, $scan, $mzxml, $noise, "{$symbol}heavy{$key}_1 ");
							if($deconvolute) $dlight1_entries2[$key] = $this->getIsotopeQuant($masses[1] - $pep_massdiff, $next_charge, $scan, $mzxml, $noise, "light{$key}_1 ", $num_isotopes);
							$mzs2[$key] = array(($masses[1] - $pep_massdiff + $next_charge * $proton_mass)/$next_charge, ($masses[1] + $next_charge * $proton_mass)/$next_charge);
						}
					}
					else {
						if(($masses[1]+$pep_massdiff)/$next_charge >= $min_mz && ($masses[1]+$pep_massdiff)/$next_charge <= $max_mz) {
							$light_entries2[$key] = $this->getIsotopeQuant($masses[1], $next_charge, $scan, $mzxml, $noise, "{$symbol}light{$key}_1 ");
							$heavy_entries2[$key] = $this->getIsotopeQuant($masses[1] + $pep_massdiff , $next_charge, $scan, $mzxml, $noise, "{$symbol}heavy{$key}_1 ");
							if($deconvolute) $dlight1_entries2[$key] = $this->getIsotopeQuant($masses[1], $next_charge, $scan, $mzxml, $noise, "light{$key}_1 ", $num_isotopes);
							$mzs2[$key] = array(($masses[1] + $next_charge * $proton_mass)/$next_charge, ($masses[1] + $pep_massdiff + $next_charge * $proton_mass)/$next_charge);
						}
					}
					if(count($light_entries2[$key]) > 0 && count($heavy_entries2[$key]) > 0) {
						$ratios_entries2[$key] = $this->printRatios($light_entries2[$key][0], $heavy_entries2[$key][0], "Scan {$scan} Pep2 {$key}_Ratios: ");

						if($deconvolute) {
							$modpeps = $pepxml->getModifiedPeptidePair();
							if($heavy) { // must change sequence to light version
								$modpeps[0] = $this->convertHeavyToLight($modpeps[0]); //preg_replace("/K\[327.13\]/", "K[325.13]", $modpeps[0]);
								$modpeps[1] = $this->convertHeavyToLight($modpeps[1]); //preg_replace("/K\[327.13\]/", "K[325.13]", $modpeps[1]);
							}
							$deconv_ratio_entries2[$key] = $this->computeRatioFromIsotopeSetsOfPeptide($dlight1_entries2[$key][0],  $modpeps[1]);
						}
				
						$ints1 = array();
						if($deconvolute) {
							$ratios_entries2[$key] = $this->computeRatioFromIsotopeSets($light_entries2[$key][0], $heavy_entries2[$key][0], 2, $modpeps[1]);
							$ints_entries2[$key] = array_pop($ratios_entries2[$key]);
						}
						else {
							$modpeps = $pepxml->getModifiedPeptidePair();
							if($heavy) { // must change sequence to light version
								$modpeps[0] = $this->convertHeavyToLight($modpeps[0]); //preg_replace("/K\[327.13\]/", "K[325.13]", $modpeps[0]);
								$modpeps[1] = $this->convertHeavyToLight($modpeps[1]); //preg_replace("/K\[327.13\]/", "K[325.13]", $modpeps[1]);
							}
							$ints_entries2[$key] = $this->getIsotopeIntensities($modpeps[1], true);
							if($ints_entries2[$key][0]==0) {
								echo "4. No theoreteical intensities obtained with calcisotopes for peptide ion \"{$modpeps[1]}\": ".join(",", $ints_entries2[$key]).", possibly due to unexpected peptide modification (contact Jimmy).  Will instead pass mass {$masses[1]}\n";
								echo "here with {$key} {$peps[0]} and {$peps[1]} and heavy ? ".($heavy ? "yes" : "no")." with peps ".join(",", $peps)."\n";
								$ints_entries2[$key] = $this->getIsotopeIntensities($masses[1], true);
							}
						}
					}
				}

				if($heavy) {
					$light2ints = $this->getIsotopeQuant($masses[1] - $pep_massdiff, $charges[1], $scan, $mzxml, $noise, "light1 ");
					$heavy2ints = $this->getIsotopeQuant($masses[1], $charges[1], $scan, $mzxml, $noise, "heavy1 ");
					if($deconvolute) $dlight2 = $this->getIsotopeQuant($masses[1] - $pep_massdiff, $charges[1], $scan, $mzxml, $noise, "light1 ", $num_isotopes);
				}
				else {
					$light2ints = $this->getIsotopeQuant($masses[1], $charges[1], $scan, $mzxml, $noise, "light1 ");
					$heavy2ints = $this->getIsotopeQuant($masses[1] + $pep_massdiff , $charges[1], $scan, $mzxml, $noise, "heavy1 ");
					if($deconvolute) $dlight2 = $this->getIsotopeQuant($masses[1], $charges[1], $scan, $mzxml, $noise, "light1 ", $num_isotopes);
				}
				$ratios2 = $this->printRatios($light2ints[0], $heavy2ints[0], "Scan {$scan} Pep2 Ratios: ");
				$ints2 = array();
				if($deconvolute) {
					$ratios2 = $this->computeRatioFromIsotopeSets($light2ints[0], $heavy2ints[0], 2, $modpeps[1]);
					$ints2 = array_pop($ratios2);
				}
				$peptide_pairratio = "";
				$deconv_ratio1 = "";
				$deconv_ratio2 = "";
				if($deconvolute) {
					$deconv_ratio1 = $this->computeRatioFromIsotopeSetsOfPeptide($dlight1[0],  $modpeps[0]);
					$deconv_ratio2 = $this->computeRatioFromIsotopeSetsOfPeptide($dlight2[0],  $modpeps[1]);
				}

				$lightreporter = $this->getIsotopeQuant(807.442528, 1, $scan, $mzxml, $noise, "lightreporter ");
				$heavyreporter = $this->getIsotopeQuant(811.455947, 1, $scan, $mzxml, $noise, "heavyreporeter ");
				$ratios3 = $this->printRatios($heavyreporter[0], $lightreporter[0], "Scan {$scan} Reporter Ratios: ");
				$next_ints = $this->getIsotopeIntensities(807.442528, true);
				$longarm_ratios = array();
				$longarm_reporter_ints = array();
				$longarm_reporter_ratios = array();
				$longarm_reporter_ids = array();
				$seen_longarm_scans = array();
				for($l = 0; $l < count($longarm_scans); $l++) {
					$next_lightreporter = $this->getIsotopeQuant(807.442528, 1, $longarm_scans[$l], $mzxml, $noise, "lightreporter ");
					$next_heavyreporter = $this->getIsotopeQuant(811.455947, 1, $longarm_scans[$l], $mzxml, $noise, "heavyreporeter ");
					$next_lightreporter[0] = array(array_sum($next_lightreporter[0]));
					$next_heavyreporter[0] = array(array_sum($next_heavyreporter[0]));
					$next_ratios = $this->printRatios($next_heavyreporter[0], $next_lightreporter[0], "Scan {$scan} Reporter Ratios: ");
					if($next_ratios[0] != "N/A") {
						array_push($longarm_reporter_ints, $next_lightreporter[0][0] + $next_heavyreporter[0][0]);
						array_push($longarm_reporter_ratios, $next_ratios[0]);
						array_push($longarm_reporter_ids, ($l == 0 ? "A" : "B") . ":LRep");
						array_push($seen_longarm_scans, $longarm_scans[$l]);
					}				

				}
				$heavy_error = $this->calculateIntensErrorVersusTheor($heavyreporter[0], $next_ints, 3);
				$light_error = $this->calculateIntensErrorVersusTheor($lightreporter[0], $next_ints, 3);
				$max_error = 0.045;
				
				$success = $heavy_error <= $max_error && $light_error <= $max_error;
				if(! $success) $ratios3 = array("N/A", "N/A", "N/A"); // cancel it
				
				if($deconvolute) $ratios3 = $this->computeRatioFromIsotopeSets($lightreporter[0], $heavyreporter[0], 4, 807.442528);
				$ints3 = array();
				if($deconvolute) {
					$ratios3 = $this->computeRatioFromIsotopeSets($lightreporter[0], $heavyreporter[0], 4, 807.442528);
					$ints3 = array_pop($ratios3);
				}
				$next_actual_ratio = 1;
				if(strpos($spectrum, "_2Sto1R_")!==false) {
					$next_actual_ratio = 0.5;
				}
				else if(strpos($spectrum, "_2Rto1S_")!==false) {
					$next_actual_ratio = 2;
				}
				fwrite($fout, $scan . "\t" . $spectrum . "\t" . $next_actual_ratio . "\t" . join("\t", $peptide_pair)."\t" . $noise . "\t");
				if($deconvolute) fwrite($fout, join(",", $ints1)."\t"); 
				fwrite($fout, number_format($light1ints[0][0], 0) . "\t" . number_format($light1ints[0][1], 0) . "\t" . number_format($light1ints[0][2], 0) . "\t");
				fwrite($fout, number_format($heavy1ints[0][0], 0) . "\t" . number_format($heavy1ints[0][1], 0) . "\t" . number_format($heavy1ints[0][2], 0) . "\t");
				if($deconvolute) {
						fwrite($fout, join("\t", $deconv_ratio1) . "\t");
				}
				else {
					fwrite($fout, join("\t", $ratios1) . "\t");
				}
				$next_deconv = $this->getPeptideHeavy2LightRatioWithTheorInts($light1ints, $heavy1ints, $this->getIsotopeIntensities($modpeps[0], true));
				
				
				fwrite($fout, $next_deconv . "\t");
	
				// compute the residuals
				$peptide_ints = array();
				$peptide_ratios = array();
				$peptide_ids = array();
				$reporter_ints = array();
				$reporter_ratios = array();
				$reporter_ids = array();
				$first = true;
				$peptide_mzs = array();
				foreach($light_entries1 as $key => $value) {
					if(true || $first) {
						$next_charge = intval($key);
						$symbol = $next_charge == $charges[0] ? "*" : "";
						fwrite($fout, "--{$symbol}{$key}-pep1-->");
						$first = false;
					}
					if(count($light_entries1[$key]) > 0 && count($heavy_entries1[$key]) > 0) {
						if($find_residuals) {
							$next_resid = array($light_entries1[$key][0][0], $light_entries1[$key][0][1], $light_entries1[$key][0][2], $heavy_entries1[$key][0][1], $heavy_entries1[$key][0][2]);
							$int_sum = array_sum($next_resid);

							if($light_entries1[$key][0][0] != "" && $light_entries1[$key][0][0] != 0 && $int_sum > 0) {
								for($r = 0; $r < 5; $r++) {
									$next_resid[$r] /= $int_sum;
								}
								$next_resid = array(($ints_entries1[$key][0]-$next_resid[0]) * ($ints_entries1[$key][0]-$next_resid[0]),
											($ints_entries1[$key][1]-$next_resid[1]) * ($ints_entries1[$key][1]-$next_resid[1]),
											($ints_entries1[$key][2]-$next_resid[2]) * ($ints_entries1[$key][2]-$next_resid[2]),
											($ints_entries1[$key][3]-$next_resid[3]) * ($ints_entries1[$key][3]-$next_resid[3]),
											($ints_entries1[$key][4]-$next_resid[4]) * ($ints_entries1[$key][4]-$next_resid[4]));		
								fwrite($fout, join(",", $next_resid) . "\t" . number_format(sqrt(array_sum($next_resid))/5, 3) . "\t");	
							}	
							else {
								fwrite($fout, "\t\t");
							}	
						}
						if($this->use_optimal_ratios) {
							$next_deconv = "N/A";
							$next_intens = $light_entries1[$key][0][0] + $light_entries1[$key][0][1] + $light_entries1[$key][0][2] + $heavy_entries1[$key][0][1] + $heavy_entries1[$key][0][2];
							if($next_intens > 0) {
								$next_peprelints = array(number_format($light_entries1[$key][0][0]/$next_intens, 3, ".", ""), number_format($light_entries1[$key][0][1]/$next_intens, 3, ".", ""), number_format($light_entries1[$key][0][2]/$next_intens, 3, ".", ""), 
									number_format($heavy_entries1[$key][0][1]/$next_intens, 3, ".", ""), number_format($heavy_entries1[$key][0][2]/$next_intens, 3, ".", ""));
						
								$next_deconv = $this->getOptimalRatio($ints_entries1[$key], $next_peprelints);
							}
						}
						else {

							$next_deconv = $this->getPeptideHeavy2LightRatioWithTheorInts($light_entries1[$key], $heavy_entries1[$key], $ints_entries1[$key]);
						}
						fwrite($fout, $next_deconv . "\t");
						if($next_deconv != "N/A" && $next_deconv !== "") {
							$next_intens = $light_entries1[$key][0][0] + $light_entries1[$key][0][1] + $light_entries1[$key][0][2] + $heavy_entries1[$key][0][1] + $heavy_entries1[$key][0][2];
							if($next_intens > 0) {
								array_push($peptide_ints, $next_intens);
								array_push($peptide_ratios, $next_deconv);
								if($samepeps || $indistinct_pepmasses) array_push($peptide_ids, "AB:{$key}");
								else array_push($peptide_ids, "A:{$key}");
								array_push($peptide_theor_ints, join(",", $ints_entries1[$key]));

								$next_peprelints = array(number_format($light_entries1[$key][0][0]/$next_intens, 3, ".", ""), number_format($light_entries1[$key][0][1]/$next_intens, 3, ".", ""), number_format($light_entries1[$key][0][2]/$next_intens, 3, ".", ""), 
									number_format($heavy_entries1[$key][0][1]/$next_intens, 3, ".", ""), number_format($heavy_entries1[$key][0][2]/$next_intens, 3, ".", ""));
								array_push($peptide_relints, join(",", $next_peprelints));

								array_push($peptide_mzs, $mzs1[$key]);
							}
						}
						
						
						fwrite($fout, join(",", $mzs1[$key]) . "\t");
						fwrite($fout, join(",", $ints_entries1[$key])."\t");
						fwrite($fout, (count($ints_entries1[$key]) > 0 && $ints_entries1[$key][0] > 0 ? number_format($ints_entries1[$key][1]/$ints_entries1[$key][0], 2) : "")."\t"); // 0 not found
						fwrite($fout, ($value[0][0] > 0 ? number_format($value[0][1]/$value[0][0], 2) : "")."\t");
						fwrite($fout, number_format($value[0][0], 0) . "\t" . number_format($value[0][1], 0) . "\t" . number_format($value[0][2], 0) . "\t");
						fwrite($fout, number_format($heavy_entries1[$key][0][0], 0) . "\t" . number_format($heavy_entries1[$key][0][1], 0) . "\t" . number_format($heavy_entries1[$key][0][2], 0) . "\t");
						if($deconvolute) {
							fwrite($fout, join("\t", $deconv_ratio_entries1[$key]) . "\t");
						}
						else {
							fwrite($fout, join("\t", $ratios_entries1[$key]) . "\t");
						}
														
					}
					else {
						fwrite($fout, "\t\t\t\t\t\t\t\t\t\t\t\t\t\t");
						if($find_residuals) {
							fwrite($fout, "\t\t");
						}
					}
				}

				if($deconvolute) fwrite($fout, join(",", $ints2)."\t"); 
				fwrite($fout, "pep2 => ".number_format($light2ints[0][0], 0) . "\t" . number_format($light2ints[0][1], 0) . "\t" . number_format($light2ints[0][2], 0) . "\t");
				fwrite($fout, number_format($heavy2ints[0][0], 0) . "\t" . number_format($heavy2ints[0][1], 0) . "\t" . number_format($heavy2ints[0][2], 0) . "\t");
				if($deconvolute) {
					fwrite($fout, join("\t", $deconv_ratio2) . "\t");
				}
				else {

					fwrite($fout, join("\t", $ratios2) . "\t");
				}
				$next_deconv = $this->getPeptideHeavy2LightRatioWithTheorInts($light2ints, $heavy2ints, $this->getIsotopeIntensities($modpeps[1], true));
				fwrite($fout, $next_deconv . "\t");

				$first = true;
				foreach($light_entries2 as $key => $value) {
					if(true || $first) {
						$next_charge = intval($key);
						$symbol = $next_charge == $charges[1] ? "*" : "";
						fwrite($fout, "--{$symbol}{$key}-pep2-->");
						$first = false;
					}
					if(count($light_entries2[$key]) > 0 && count($heavy_entries2[$key]) > 0) {
						if($find_residuals) {
							$next_resid = array($light_entries2[$key][0][0], $light_entries2[$key][0][1], $light_entries2[$key][0][2], $heavy_entries2[$key][0][1], $heavy_entries2[$key][0][2]);
							$int_sum = array_sum($next_resid);
							if($light_entries2[$key][0][0] != "" && $light_entries2[$key][0][0] != 0 && $int_sum > 0) {
								for($r = 0; $r < 5; $r++) $next_resid[$r] =  $next_resid[$r] / $int_sum;
									$next_resid = array(($ints_entries2[$key][0]-$next_resid[0]) * ($ints_entries2[$key][0]-$next_resid[0]),
											($ints_entries2[$key][1]-$next_resid[1]) * ($ints_entries2[$key][1]-$next_resid[1]),
											($ints_entries2[$key][2]-$next_resid[2]) * ($ints_entries2[$key][2]-$next_resid[2]),
											($ints_entries2[$key][3]-$next_resid[3]) * ($ints_entries2[$key][3]-$next_resid[3]),
											($ints_entries2[$key][4]-$next_resid[4]) * ($ints_entries2[$key][4]-$next_resid[4]));		
											
								fwrite($fout, join(",", $next_resid) . "\t" . number_format(sqrt(array_sum($next_resid))/5, 3) . "\t");			
							}	
							else {
								fwrite($fout, "\t\t");
							}
						}		
						
						// here make sure relative isotope peak intensities consistent with theoretical
						if($this->use_optimal_ratios) {
							$next_deconv = "N/A";
							$next_intens = $light_entries2[$key][0][0] + $light_entries2[$key][0][1] + $light_entries2[$key][0][2] + $heavy_entries2[$key][0][1] + $heavy_entries2[$key][0][2];
							if($next_intens > 0) {
								$next_peprelints = array(number_format($light_entries2[$key][0][0]/$next_intens, 3, ".", ""), number_format($light_entries2[$key][0][1]/$next_intens, 3, ".", ""), number_format($light_entries2[$key][0][2]/$next_intens, 3, ".", ""), 
									number_format($heavy_entries2[$key][0][1]/$next_intens, 3, ".", ""), number_format($heavy_entries2[$key][0][2]/$next_intens, 3, ".", ""));
						
								$next_deconv = $this->getOptimalRatio($ints_entries2[$key], $next_peprelints);
							}
						}
						else {
							$next_deconv = $this->getPeptideHeavy2LightRatioWithTheorInts($light_entries2[$key], $heavy_entries2[$key], $ints_entries2[$key]);
						}
						fwrite($fout, $next_deconv . "\t");
						if($next_deconv != "N/A" && $next_deconv !== "" && ! $samepeps && ! $indistinct_pepmasses) {
							$next_intens = $light_entries2[$key][0][0] + $light_entries2[$key][0][1] + $light_entries2[$key][0][2] + $heavy_entries2[$key][0][1] + $heavy_entries2[$key][0][2];
							if($next_intens > 0) {
								array_push($peptide_ints, $next_intens);
								array_push($peptide_ratios, $next_deconv);
								array_push($peptide_ids, "B:{$key}");
								array_push($peptide_theor_ints, join(",", $ints_entries2[$key]));

								$next_peprelints = array(number_format($light_entries2[$key][0][0]/$next_intens, 3, ".", ""), number_format($light_entries2[$key][0][1]/$next_intens, 3, ".", ""), number_format($light_entries2[$key][0][2]/$next_intens, 3, ".", ""), 
									number_format($heavy_entries2[$key][0][1]/$next_intens, 3, ".", ""), number_format($heavy_entries2[$key][0][2]/$next_intens, 3, ".", ""));
								array_push($peptide_relints, join(",", $next_peprelints));

								array_push($peptide_mzs, $mzs2[$key]);
							}
						}
						fwrite($fout, join(",", $mzs2[$key]) . "\t");
						fwrite($fout, join(",", $ints_entries2[$key])."\t");
						fwrite($fout, (count($ints_entries2[$key]) > 0 && $ints_entries2[$key][0] > 0 ? number_format($ints_entries2[$key][1]/$ints_entries2[$key][0], 2) : "")."\t"); // NOT FOUND
						fwrite($fout, ($value[0][0] > 0 ? number_format($value[0][1]/$value[0][0], 2) : "")."\t");
						fwrite($fout, number_format($value[0][0], 0) . "\t" . number_format($value[0][1], 0) . "\t" . number_format($value[0][2], 0) . "\t");
						fwrite($fout, number_format($heavy_entries2[$key][0][0], 0) . "\t" . number_format($heavy_entries2[$key][0][1], 0) . "\t" . number_format($heavy_entries2[$key][0][2], 0) . "\t");
						if($deconvolute) {
							fwrite($fout, join("\t", $deconv_ratio_entries2[$key]) . "\t");
						}
						else {
							fwrite($fout, join("\t", $ratios_entries2[$key]) . "\t");
						}
					}
					else {
						fwrite($fout, "\t\t\t\t\t\t\t\t\t\t\t\t\t\t");
						if($find_residuals) {
							fwrite($fout, "\t\t");
						}
					}
				}
if($samepeps) {
				$first = true;
				foreach($light_entries1 as $key => $value) {
					if(true || $first) {
						fwrite($fout, "\t");
						$first = false;
					}
					fwrite($fout, "\t\t\t\t\t\t\t\t\t\t\t\t\t");
					if($find_residuals) {
						fwrite($fout, "\t\t");
					}
				}

}
				fwrite($fout, "--Reporter---");
				if($deconvolute) fwrite($fout, join(",", $ints3)."\t"); 
				fwrite($fout, number_format($lightreporter[0][0], 0) . "\t" . number_format($lightreporter[0][1], 0) . "\t" . number_format($lightreporter[0][2], 0) . "\t");
				fwrite($fout, number_format($heavyreporter[0][0], 0) . "\t" . number_format($heavyreporter[0][1], 0) . "\t" . number_format($heavyreporter[0][2], 0) . "\t");
				fwrite($fout, join("\t", $ratios3));
	
				for($z = 0; $z < 3; $z++) {
					if($ratios3[$z] != "N/A") {
						// can also check for relative intensities of first and second
					
						array_push($reporter_ints, $lightreporter[0][$z] + $heavyreporter[0][$z]);
						array_push($reporter_ratios, $ratios3[$z]);
						array_push($reporter_ids, "R:".($z+1));
					}
				}
				
				
				
				
	
				if($deconvolute) {
					fwrite($fout, "\t" . number_format($light1ints[0][0]/$ints1[0], 0, ",", "") . "\t" . number_format($light1ints[0][1]/$ints1[1], 0, ",", "") . "\t" 
						. number_format($light2ints[0][0]/$ints2[0], 0, ",", "") . "\t" . number_format($light2ints[0][1]/$ints2[1], 0, ",", ""));
				}
				foreach($next_hk_info as $key => $value) {
					// light mass is key
					for($charge = 1; $charge <= 3; $charge++) {
						if(($key+$pep_massdiff)/$charge >= 500 && ($key+$pep_massdiff)/$charge <= 2000) {
							$light_next = $this->getIsotopeQuant($key, $charge, $scan, $mzxml, $noise, "light{$key}_1 ");
							$heavy_next = $this->getIsotopeQuant($key + $pep_massdiff , $charge, $scan, $mzxml, $noise, "heavy{$key}_1 ");
							if(count($light_next) > 0 && count($heavy_next) > 0 && array_sum($light_next[0]) > 0 && array_sum($heavy_next[0]) > 0) {
							$ratios_next = $this->printRatios($light_next[0], $heavy_next[0], "Scan {$scan} Pep1 {$key}_Ratios: ");

							if($deconvolute) {
								$modpeps = $pepxml->getModifiedPeptidePair();
								$deconv_ratio_next = $this->computeRatioFromIsotopeSetsOfPeptide($light_next[0],  $key);

							}
							$ints1 = array();
							if($deconvolute) {
								$ratios_next = $this->computeRatioFromIsotopeSets($light_next[0], $heavy_next[0], 2, $key);
								$ints_next = array_pop($ratios_next);
							
							}
							else {
								$modpeps = $pepxml->getModifiedPeptidePair();
								$ints_next = $this->getIsotopeIntensities($key, true);
							}
							if($deconvolute) $dlight_next = $this->getIsotopeQuant($key, $charge, $scan, $mzxml, $noise, "light{$key}_1 ", $num_isotopes);
							$mzs2_next = array(($key + $next_charge * $proton_mass)/$charge, ($key + $pep_massdiff + $charge * $proton_mass)/$charge);
							$ints_next = $this->getIsotopeIntensities($key, true);

							$light2ints = $this->getIsotopeQuant($key, 1, $scan, $mzxml, $noise, "light1 ");
							$heavy2ints = $this->getIsotopeQuant($key + $pep_massdiff , 1, $scan, $mzxml, $noise, "heavy1 ");
							if($deconvolute) $dlight2 = $this->getIsotopeQuant($key, 1, $scan, $mzxml, $noise, "light1 ", $num_isotopes);
							fwrite($fout, "\t");
							fwrite($fout, $this->getPeptideHeavy2LightRatio($light_next, $heavy_next) . "\t");
							fwrite($fout, join(",", $mzs2_next) . "\t");
							fwrite($fout, join(",", $ints_next)."\t");
							fwrite($fout, ($ints_next[0] > 0 ? number_format($ints_next[1]/$ints_next[0], 2) : "")."\t");
							fwrite($fout, ($light_next[0][0] > 0 ? number_format($light_next[0][1]/$light_next[0][0], 2) : "")."\t");
							fwrite($fout, number_format($light_next[0][0], 0) . "\t" . number_format($light_next[0][1], 0) . "\t" . number_format($light_next[0][2], 0) . "\t");
							fwrite($fout, number_format($heavy_next[0][0], 0) . "\t" . number_format($heavy_next[0][1], 0) . "\t" . number_format($heavy_next[0][2], 0) . "\t");
							if($deconvolute) {
								fwrite($fout, join("\t", $deconv_ratio_next) . "\t");
							}
							else {
								fwrite($fout, join("\t", $ratios_next) . "\t");
								$next_run = $mzxml;
								$nextpos = strrpos($next_run, "/");
								if($nextpos!==false) $next_run = substr($next_run, $nextpos+1);
								$next_run = substr($next_run, 0, strlen($next_run) - 6);
								

								echo $next_run . "\t" . $scan . "\t" . $pepxml->getPreservedOrderedModifiedPeptidePair()."\t" . number_format($masses[0], 2) . "\t" . 
								number_format($masses[1], 2) . "\t" . 
								$key . "\t" . $charge . "\t" . $value[0] . "\t" . $this->getPeptideHeavy2LightRatio($light_next, $heavy_next) .
									"\t" . join("\t", $light_next[0])."\t".$heavy_next[0][1] . "\t" . $heavy_next[0][2] . "\t" . join(",", $ints_next) . "\n";
							}
						} // count greater for light and heavy
					}


				} // next charge
			}
			fwrite($fout, "\t{$next_actual_ratio}");
			$reporter_ratios2 = array();
			$reporter_ints2 = array();
			$reporter_ids2 = array();
			for($k = 0; $k < count($reporter_ratios); $k++) {
				if($reporter_ratios[$k]=="" || $reporter_ratios[$k] <= 0 || $reporter_ratios[$k]=="inf") continue;   // 011020
				array_push($reporter_ratios2, $reporter_ratios[$k]);
				array_push($reporter_ints2, $reporter_ints[$k]);
				array_push($reporter_ids2, $reporter_ids[$k]);
			}
			$longarm_reporter_ratios2 = array();
			$longarm_reporter_ints2 = array();
			$longarm_reporter_ids2 = array();
			$longarm_scans2 = array();
			for($k = 0; $k < count($longarm_reporter_ratios); $k++) {
				if($longarm_reporter_ratios[$k]=="" || $longarm_reporter_ratios[$k] <= 0 || $longarm_reporter_ratios[$k]=="inf") continue;   // 011020
				array_push($longarm_reporter_ratios2, $longarm_reporter_ratios[$k]);
				array_push($longarm_reporter_ints2, $longarm_reporter_ints[$k]);
				array_push($longarm_reporter_ids2, $longarm_reporter_ids[$k]);
				array_push($longarm_scans2, $seen_longarm_scans[$k]);
			}
			$peptide_ratios2 = array();
			$peptide_ints2 = array();
			$peptide_ids2 = array();
			$peptide_mzs2 = array();
			$peptide_theorints2 = array();
			$peptide_relints2 = array();
			for($k = 0; $k < count($peptide_ratios); $k++) {
				if($peptide_ratios[$k]=="" || $peptide_ratios[$k] <= 0 || $peptide_ratios[$k]=="inf") continue;   // 011020
				array_push($peptide_ratios2, $peptide_ratios[$k]);
				array_push($peptide_ints2, $peptide_ints[$k]);
				array_push($peptide_ids2, $peptide_ids[$k]);
				array_push($peptide_mzs2, join("-", $peptide_mzs[$k]));
				array_push($peptide_theorints2, $peptide_theor_ints[$k]);
				array_push($peptide_relints2, $peptide_relints[$k]);
			}

			fwrite($fout, "\t" . join(";", $reporter_ids2));
			fwrite($fout, "\t" . join(";", array_map("log2ratio", $reporter_ratios2)) . "\t" . join(";", $reporter_ints2));

			fwrite($fout, "\t" . join(";", $peptide_ids2));
			fwrite($fout, "\t" . join(";", array_map("log2ratio", $peptide_ratios2)) . "\t" . join(";", $peptide_ints2));
			fwrite($fout, "\t" . join(";", $peptide_mzs2). "\t" . join(";", $peptide_theorints2) . "\t" . join(";", $peptide_relints2));


			fwrite($fout, "\t" . join(";", $frag_fragments));
			fwrite($fout, "\t" . join(";", array_map("log2ratio", $frag_ratios)) . "\t" . join(";", $frag_ints));
			fwrite($fout, "\t" . join(";", $frag_mzs). "\t" . join(";", $fragment_theor_ints) . "\t" . join(";", $fragment_relints));
			fwrite($fout, "\t" . join(";", $frag_scans));
			fwrite($fout, "\t" . $frag_noise);
			fwrite($fout, "\t" . join(";", $longarm_reporter_ids2));
			fwrite($fout, "\t" . join(";", array_map("log2ratio", $longarm_reporter_ratios2)) . "\t" . join(";", $longarm_reporter_ints2));
			fwrite($fout, "\t" . join(";", $longarm_scans2));

if(! $this->use_optimal_ratios) {			
			fwrite($fout, "\t");
			$optimal_pep_ratios = array();
			for($z = 0; $z < count($peptide_relints2); $z++) {
				$next_th = explode(",", $peptide_theorints2[$z]);
				$next_ob = explode(",", $peptide_relints2[$z]);
				array_push($optimal_pep_ratios, $this->getOptimalRatio($next_th, $next_ob));
			}
			fwrite($fout, join(";", array_map("log2ratio", $optimal_pep_ratios)));
			fwrite($fout, "\t");
			if(count($nextFragInfo) > 0) {
				$optimal_frag_ratios = array();
				for($z = 0; $z < count($fragment_relints); $z++) {
					$next_th = explode(",", $fragment_theor_ints[$z]);
					$next_ob = explode(",", $fragment_relints[$z]);
					array_push($optimal_frag_ratios, $this->getOptimalRatio($next_th, $next_ob));
				}
				fwrite($fout, join(";", array_map("log2ratio", $optimal_frag_ratios)));
			}
}	// only if not use otpimal		
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
			fwrite($fout, "\t" . $nextResA . "-" . $nextResB);
			fwrite($fout, "\t" . $protein_seqs_and_genes[$peptide_pair[1]]['gene'] . "-" . $protein_seqs_and_genes[$peptide_pair[4]]['gene']);			
			fwrite($fout, "\t" . $noise . "\n");

			
			}
			$peps = array();
			$masses = array(); 
			$charges = array();
			$next_scans = array();
			$ms3scans = array();
		}

	}
	fclose($fout);
	echo "Results written to {$outfile}\n";
	echo "\n";
	
	$fout_lh = fopen($this->scan_ratiofile_lh, "w");
	fwrite($fout_lh, "spectrum\tpepA\tproA\tkposA\tpepB\tproB\tkposB\tlight_id\theavy_id\n");
	foreach($seen as $xl => $light_heavy) {
		fwrite($fout_lh, $xl . "\t" . join("\t", $light_heavy) . "\n");
	}
	fclose($fout_lh);
	echo "Light-heavy identification info written to {$this->scan_ratiofile_lh}\n";
	echo "\n";
	
	return;
	foreach($bar_graphs_peptide as $key => $value) {
		$nextval = $report_log2ratios ? $key : log(floatval($key))/log(2);
		$bar_graphs_peptide[$key]->printDistr(STDOUT, "Peptides Log2={$nextval}", true);
		$bar_graphs_reporter[$key]->printDistr(STDOUT, "Reporter Log2={$nextval}", true);
		$bar_graphs_combined[$key]->printDistr(STDOUT, "Peptides and Reporter Log2={$nextval}", true);
		echo "\n";
	}
}

private function getWeightedAverageRatio($ratios, $ints = NULL, $outlier_remove = false, $norm_offset = 0, $verbose = false, $pvalue_with_tstat = false, $no_output_wt = false) {
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
		// now get weighted average of the remainder
		$df = count($vals) - 1 ;
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
			$output = array(number_format($mean, 2, ".", ""), count($vals));
			if(! $no_output_wt) array_push($output, $total_wt);
			array_push($output, number_format($std, 3, ".", ""));
			array_push($output, $this->getPvalue($mean, $std, $df, $pvalue_with_tstat));
			array_push($output, join(",", $rejects));
			return $output;
		}
		else if($total_wt > 0 && count($vals) === 1) {
			$output = array(number_format($vals[0][0], 2, ".", ""), 1);
			if(! $no_output_wt) array_push($output, $vals[0][1]);	
			array_push($output, 0);
			array_push($output, $default_pval);
			array_push($output, "");
			return $output;
		}
		$output =  array("", 0);
		if(! $no_output_wt) array_push($output, 0);			
		array_push($output, 0);
		array_push($output, $default_pval);
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
		
		$output = array(number_format($mean, 2, ".", ""), count($ratios), $total_wt, number_format($std, 3, ".", ""), 
			$this->getPvalue($mean, $std, count($ratios)-1, $pvalue_with_tstat));
		return $output;
	}
	if($total_wt > 0 && count($ratios) === 1)  return array($ratios[0], 1, $ints[0], 0, $default_pval);
	return array("", 0, 0, 0, $default_pval);


}

private function getWeightAverageWithRatio($ratios, $ints, $ratios2 = NULL, $ints2 = NULL, $use_log2 = false, $outlier_remove = false) {
	$total = 0;
	$total_wt = 0;
	$total_contributing = 0;
	$vals = array();
	if($outlier_remove) { // seemed to have little effect on outcome
		for($k = 0;  $k < count($ratios); $k++) {
			if($ratios[$k] > 0 && $ratios[$k] < 999) {
				if($use_log2) $ratios[$k] = log($ratios[$k])/log(2);
				array_push($vals, array($ratios[$k], $ints[$k]));
			}
		}
		if(! is_null($ratios2)) {
			for($k = 0;  $k < count($ratios2); $k++) {
				if($ratios2[$k] > 0 && $ratios2[$k] < 999) {
				if($use_log2) $ratios2[$k] = log($ratios2[$k])/log(2);
					array_push($vals, array($ratios2[$k], $ints2[$k]));
				}
			}
		}	
		$rejects = array();
		MyStatistics::removeOutliersByIndex($vals, $rejects, 0, false); // index 1
		// now get weighted average of the remainder
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
		
			$output = array(number_format($mean, 2, ".", ""), count($vals), $total_wt, number_format($std, 3, ".", ""));
			if($outlier_remove) array_push($output, join(",", $rejects));
			return $output;
		}
		$output =  array("", 0, 0, 0);
		if($outlier_remove) array_push($output, "");
		return $output;
	}
	
	
	for($k = 0;  $k < count($ratios); $k++) {
		if($ratios[$k] > 0 && $ratios[$k] < 999) {
			if($use_log2) $ratios[$k] = log($ratios[$k])/log(2);
			$total += $ratios[$k] * $ints[$k];
			$total_wt += $ints[$k];
			$total_contributing++;
		}
	}
	if(! is_null($ratios2)) {
		for($k = 0;  $k < count($ratios2); $k++) {
			if($ratios2[$k] > 0 && $ratios2[$k] < 999) {
				if($use_log2) $ratios2[$k] = log($ratios2[$k])/log(2);
				$total += $ratios2[$k] * $ints2[$k];
				$total_wt += $ints2[$k];
				$total_contributing++;
			}
		}
	}
	
	if($total_wt > 0)  {
		$std = 0;
		$mean = $total / $total_wt;
		for($k = 0; $k < count($ratios); $k++) {
			$std += $ints[$k] * ($ratios[$k] - $mean) * ($ratios[$k] - $mean);
		}
		if(count($ratios) > 1) $std *= count($ratios) / ((count($ratios)-1) * $total_wt);
		$std = sqrt($std);
		$output = array(number_format($mean, 2, ".", ""), $total_contributing, $total_wt, number_format($std, 3, ".", ""));
		return $output;
	}
	return array("", 0, 0, 0);
}


private function getPvalue($mean, $std, $df, $pval_with_tstat = false) {
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
	$next = explode(" ", $modpep);
	if(count($next)!==2) {
		echo "Error: have ".count($next)." crosslink modifictaions in {$modpep} (orig {$copy})\n";
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

public function getLightHeavyInfo() {
$light_heavy_info = array();
if(! file_exists($this->scan_ratiofile_lh)) {
	echo "Warning: light heavy file {$this->scan_ratiofile_lh} does not exist\n";
	return $light_heavy_info;
}
$first = true;
$myfile = fopen($this->scan_ratiofile_lh, "r");
while($line = rtrim(fgets($myfile))) {
	if($first) {
		$first = false;
	}
	else {
		$next = explode("\t", $line);
		if($this->deadend) $next_xlid = join("_", $this->getModifiedPeptideWithoutCrosslinkmodAndKpos($next[1])); 
		else $next_xlid = join("_", $this->getModifiedPeptideWithoutCrosslinkmodAndKpos($next[1])) . "_" .
			join("_", $this->getModifiedPeptideWithoutCrosslinkmodAndKpos($next[4]));
		$light_heavy_info[$next[0] . ":" . $next_xlid] = array($next[7], $next[8]);
		
	}
}
return $light_heavy_info;
}


public function computeCombinedCrosslinkQuant($samples = array()) {

$file = $this->scan_ratiofile;
if(! file_exists($file)) {
	echo "Error: {$file} does not exist.  Please make sure to run RUN step first before LOGRATIOS\n";
	exit(1);
}
if(array_key_exists('use_current_normalization',  $this->iqpir_params) &&  
	$this->iqpir_params['use_current_normalization'] == "true") $this->setUseCurrentNorm();
$normfile = $this->iqpir_params['iqpir_output_dir'] . "sample_biorep_normfactors.txt";
$offset = $this->deadend ? -55 : 0; // account for missing peptide 2 info in deadend files
echo "setting normalization file to {$normfile} with normalization {$this->iqpir_params['normalize']}\n";
// set normalization factor to $this->iqpir_params['sample']['biorep']....
foreach($this->iqpir_params['samples'] as $sample => $bioreps) {
	if(count($samples) > 0 && ! array_key_exists($sample, $samples)) continue;

echo "Ready to analyze {$file} for sample {$sample} with bioreps ".join(",", array_keys($bioreps))."\n";
	if($this->iqpir_params['normalize'] == "true") $sample_biorepnorms[$sample] = array(); // will accumulate all reatios

	$first = true;
	$data = array(); // hashed by crosslink, holds peptide ratios, ints, frag ratios, ints, reporter ratios, ints .... and added summary information of contributing pep and frags
	$reporter_data = array();
	$myfile = fopen($file, "r");
	$xl_prots = array();
	$protpair_unis = array();
	$xl_residuepairs = array();
	$output_genes = true; // rather than prots
	$subtract_noise = false;
	$xlinkdb_ids = array();
	$pepfragratios = array();
	$pepfragoptimalratios = array();
	$samplefile = $this->iqpir_params['iqpir_output_dir'] . "sample." . $sample . ".ratios.txt";
	$chromfile = $this->iqpir_params['iqpir_output_dir'] . "sample." . $sample . ".chroms.txt";
	$respairfile = $this->iqpir_params['iqpir_output_dir'] . "sample." . $sample . ".respairs.txt";
	$protpairfile = $this->iqpir_params['iqpir_output_dir'] . "sample." . $sample . ".protpairs.txt";
	$reshubfile = $this->iqpir_params['iqpir_output_dir'] . "sample." . $sample . ".reshubs.txt";
	$numidsfile = $this->iqpir_params['iqpir_output_dir'] . "sample." . $sample . ".numids.txt";
	$pepfraginfo_keys = array();
	$sample_biorepnorms = array();
	$sample_optimal_biorepnorms = array();
	$num_RH_SH_ids = array();
	$light_heavy_info = $this->getLightHeavyInfo();
	$include_reporters = $this->deadend && array_key_exists("include_reporters", $this->iqpir_params) && 
		$this->iqpir_params['include_reporters'] == "true";
	$include_longarm_reporters = array_key_exists('longarm_reporters', $this->iqpir_params) && 
					$this->iqpir_params['longarm_reporters'] == "true";
	$line_no = 0;
	while($line=fgets($myfile)){
		$line = preg_replace( "/\r|\n/", "", $line);

		$line_no ++;
		$next = explode("\t", $line);
		if($first) {
			if($next[123 + $offset]!=="actual_log2ratio") {
				echo "\n\nYour {$file} was generated by an old version of this program.\nAs a result, you must re-execute the RUN step before LOGRATIOS.  Sorry for the inconvienience.\n\n";
				exit(1);
			}
			$first =  false;
			continue;
		}
		if(count($next) < 128 + $offset) continue;
		
		if($this->intra_only) {
			$prots = explode("-", $next[count($next)-2]);
			if($prots[0]!=$prots[1]) continue;
			$res = explode("-",  $next[count($next)-3]);
			if($res[0]==$res[1]) continue;
		}
		
		$next_xl = $this->deadend ? $this->convertHeavyToLight($next[3]) : 
			$this->convertHeavyToLight($next[3] . "_" . $next[6]);
		$rawfile = $next[1];
		$nextpos = strpos($rawfile, ".");
		if($nextpos!==false) $rawfile = substr($rawfile, 0, $nextpos);
		if(! array_key_exists($rawfile, $this->iqpir_params['rawfiles']) || 
			$this->iqpir_params['rawfiles'][$rawfile]['sample']!= $sample) continue; // not now
		if($this->iqpir_params['normalize'] == "true" && 
			! array_key_exists($this->iqpir_params['rawfiles'][$rawfile]['biorep'], $sample_biorepnorms)) {
			$sample_biorepnorms[$this->iqpir_params['rawfiles'][$rawfile]['biorep']] = array();
			$sample_optimal_biorepnorms[$this->iqpir_params['rawfiles'][$rawfile]['biorep']] = array();
		}
		$next_noise = $next[count($next)-1]; // the last column is noise
			$protpair_unis[$next[count($next)-2]] = $next[4];
			if(! $this->deadend) $protpair_unis[$next[count($next)-2]] .= "-" . $next[7];
			if(! $this->deadend) {
				$nextprots = explode("-", $next[count($next)-2]);
				$protpair_unis[$nextprots[1] . "-" . $nextprots[0]] = $next[7] . "-" . $next[4];
			}
			$xl_prots[$next_xl] = $next[count($next)-2]; //$prot_genes[$next[4]] . "-" . $prot_genes[$next[7]];
			$xl_residuepairs[$next_xl] = $next[count($next)-3];
		}

		// can add on spectrum, pep_ids, pep_mzs, frag_ids, frag_mzs
		$next_pep_ids = $next[127 + $offset] == "" ? array() : explode(";", $next[127 + $offset]);
		$next_pep_mzs = explode(";", $next[130 + $offset]);
		$next_frag_mzs = explode(";", $next[136 + $offset]);
		$next_frag_ids = $next[133 + $offset] == "" ? array() : explode(";", $next[133 + $offset]);
		$next_frag_theorints = explode(";", $next[137 + $offset]);
		$next_pep_theorints = explode(";", $next[131 + $offset]);
		$next_frag_obs_relints = explode(";", $next[138 + $offset]);
		$next_frag_scans = array(); //explode(";", $next[139 + $offset]);
		$next_frag_noise = $next_noise; //$next[140 + $offset];
		if(count($next) < 143 + $offset) {
			for($k = 0; $k < count($next_frag_ids); $k++) {
				array_push($next_frag_scans, $next[0]);
			}
		}
		else {
			$next_frag_scans = explode(";", $next[139 + $offset]);
			$next_frag_noise = $next[140 + $offset];
		}
		$next_larm_rep_ids = array();
		$next_larm_rep_ratios = array();
		$next_larm_rep_ints = array();
		$next_larm_rep_scans = array();
		if($include_longarm_reporters) {
			$next_larm_rep_ids = $next[141 + $offset] == "" ? array() : explode(";", $next[141 + $offset]);
			$next_larm_rep_ratios = explode(";", $next[142 + $offset]);
			$next_larm_rep_ints = explode(";", $next[143 + $offset]);
			$next_larm_rep_scans = explode(";", $next[144 + $offset]);
		}
		
		
		$next_pep_obs_relints = explode(";", $next[132 + $offset]);
		$xl_peps = explode("_", $next_xl);
		$next_reporter_ids = ! $include_reporters ? array() : explode(";", $next[124 + $offset]);
		for($z = 0; $z < count($xl_peps); $z++) {
			$xl_peps[$z] = join("_", $this->getModifiedPeptideWithoutCrosslinkmodAndKpos($xl_peps[$z]));
		}
		$xlinkdb_ids[$next_xl] = join("_", $xl_peps);
		$scan = $next[0];
		$next_spec = substr($next[1], 0, strlen($next[1])-2);
		// here can count up crosslinks heavy/light
		if(count($light_heavy_info) > 0) {
			if(! array_key_exists($next_xl, $num_RH_SH_ids)) $num_RH_SH_ids[$next_xl] = array("RH" => 0, "SH" => 0, "log2ratio" => "");
			if(array_key_exists($next_spec . ":" . $xlinkdb_ids[$next_xl], $light_heavy_info)) {
				if($this->iqpir_params['rawfiles'][$rawfile]['orientation'] == "rev") {
					$num_RH_SH_ids[$next_xl]["RH"] += $light_heavy_info[$next_spec . ":" . $xlinkdb_ids[$next_xl]][1];
					$num_RH_SH_ids[$next_xl]["SH"] += $light_heavy_info[$next_spec . ":" . $xlinkdb_ids[$next_xl]][0];
				}
				else {
					$num_RH_SH_ids[$next_xl]["RH"] += $light_heavy_info[$next_spec . ":" . $xlinkdb_ids[$next_xl]][0];
					$num_RH_SH_ids[$next_xl]["SH"] += $light_heavy_info[$next_spec . ":" . $xlinkdb_ids[$next_xl]][1];
				}
			}
		}
		if(count($next_pep_ids)==0 && count($next_frag_ids) == 0 && count($next_reporter_ids) == 0) continue;

		if(! array_key_exists($next_xl, $data)) $data[$next_xl] = array(); 
		$next_info = array("pep_log2ratios" => explode(";", $next[128 + $offset]), "pep_ints" => explode(";", $next[129 + $offset]), "frag_log2ratios" => explode(";", $next[134 + $offset]), 
			"frag_ints" => explode(";", $next[135 + $offset]), "reporter_log2ratios" => explode(";", $next[125 + $offset]), "reporter_ints" => explode(";", $next[126 + $offset]));
		if(! $this->use_optimal_ratios) {
			$next_info["pep_optimal_ratios"] = explode(";", $next[139 + $offset]);
			$next_info["frag_optimal_ratios"] = explode(";", $next[140 + $offset]);
		}
		$multiplier = $this->iqpir_params['rawfiles'][$rawfile]['orientation'] == "rev" ? -1 : 1;

		if(count($next_reporter_ids) > 0) {
			if(! array_key_exists($next_xl, $reporter_data)) $reporter_data[$next_xl] = array();
			$next_reporter_ints = explode(";", $next[126+$offset]);
			$next_reporter_ratios = explode(";", $next[125+$offset]);
			for($r = 0; $r < count($next_reporter_ids); $r++) {
				if($this->min_sig2noise > 0 && $next_noise > 0 && $next_reporter_ints[$r] < 5 * $this->min_sig2noise * $next_noise) continue;
				array_push($reporter_data[$next_xl], array($next_reporter_ratios[$r] * $multiplier, 
					$next_reporter_ids[$r], $next_reporter_ints[$r], $this->iqpir_params['rawfiles'][$rawfile]['biorep']));
			}
		}

		for($k = 0; $k < count($next_pep_ids); $k++) {
			$excl = "FALSE";
			if($this->min_sig2noise > 0 && $next_noise > 0 && $next_info['pep_ints'][$k] < 5 * $this->min_sig2noise * $next_noise) {
				echo "Ignoring peptide ratio of {$next_pep_ids[$k]} of {$next_xl} due to intensity {$next_info['pep_ints'][$k]} being less than ".(5*$this->min_sig2noise)." times noise {$next_noise} (".($this->min_sig2noise * 5*$next_noise).") in scan {$next[0]}\n";
				$excl = "SIG2NOISE"; // avoid the ints
			}
			$next_info['pep_log2ratios'][$k] *= $multiplier;
		    if(! $this->use_optimal_ratios) {
				$next_info['pep_optimal_ratios'][$k] *= $multiplier;
				array_push($data[$next_xl], array("rawfile" => $rawfile, "scan" => $scan, "id" => $next_pep_ids[$k], "logratio" =>  $next_info['pep_log2ratios'][$k], "mzs" => $next_pep_mzs[$k], 
				"total_intensity" => $subtract_noise && $next_info['pep_ints'][$k] > $next_noise * 5 ? $next_info['pep_ints'][$k] - $next_noise * 5 : $next_info['pep_ints'][$k], 
				"excluded" => $excl, "theor_relints" => $next_pep_theorints[$k], "obs_relints" => $next_pep_obs_relints[$k], "biorep" => $this->iqpir_params['rawfiles'][$rawfile]['biorep'],
				"orientation" => $this->iqpir_params['rawfiles'][$rawfile]['orientation'], "optimal_logratio" => $next_info['pep_optimal_ratios'][$k], "noise" => $next_noise));
			}
			else {
			if(count($pepfraginfo_keys)===0) $pepfraginfo_keys = array_keys($data[$next_xl][0]);
			if($this->iqpir_params['normalize'] == "true" && $excl == "FALSE") {
				array_push($sample_biorepnorms[$this->iqpir_params['rawfiles'][$rawfile]['biorep']], $next_info['pep_log2ratios'][$k]);
				if(! $this->use_optimal_ratios) array_push($sample_optimal_biorepnorms[$this->iqpir_params['rawfiles'][$rawfile]['biorep']], $next_info['pep_optimal_ratios'][$k]);
			}
		}
		for($k = 0; $k < count($next_frag_ids); $k++) {
			$excl = "FALSE";
			if($this->min_sig2noise > 0 && $next_frag_noise > 0 && $next_info['frag_ints'][$k] < 5 * $this->min_sig2noise * $next_frag_noise) {
				echo "Ignoring fragment ratio of {$next_frag_ids[$k]} of {$next_xl} due to intensity {$next_info['frag_ints'][$k]} being less than ".(5*$this->min_sig2noise)." times noise {$next_frag_noise} (".
					($this->min_sig2noise * 5*$next_frag_noise).") in {$rawfile} scan {$next_frag_scans[$k]}\n";
				$excl = "SIG2NOISE"; // avoid the ints
			}
			$next_info['frag_log2ratios'][$k] *= $multiplier;
			if(! $this->use_optimal_ratios) {
				$next_info['frag_optimal_ratios'][$k] *= $multiplier;
				array_push($data[$next_xl], array("rawfile" => $rawfile, "scan" => $next_frag_scans[$k], "id" => $next_frag_ids[$k], "logratio" =>  $next_info['frag_log2ratios'][$k], "mzs" => $next_frag_mzs[$k], 
				"total_intensity" => $subtract_noise && $next_info['frag_ints'][$k] > $next_frag_noise * 5 ? $next_info['frag_ints'][$k] - $next_frag_noise[$k] * 5 : $next_info['frag_ints'][$k], 
				"excluded" => $excl, "theor_relints" => $next_frag_theorints[$k], "obs_relints" => $next_frag_obs_relints[$k], "biorep" => $this->iqpir_params['rawfiles'][$rawfile]['biorep'],
				"orientation" => $this->iqpir_params['rawfiles'][$rawfile]['orientation'], "optimal_logratio" => $next_info['frag_optimal_ratios'][$k], "noise" => $next_frag_noise));
			}
			else {
				array_push($data[$next_xl], array("rawfile" => $rawfile, "scan" => $next_frag_scans[$k], "id" => $next_frag_ids[$k], "logratio" =>  $next_info['frag_log2ratios'][$k], "mzs" => $next_frag_mzs[$k], 
				"total_intensity" => $subtract_noise && $next_info['frag_ints'][$k] > $next_frag_noise * 5 ? $next_info['frag_ints'][$k] - $next_frag_noise[$k] * 5 : $next_info['frag_ints'][$k], 
				"excluded" => $excl, "theor_relints" => $next_frag_theorints[$k], "obs_relints" => $next_frag_obs_relints[$k], "biorep" => $this->iqpir_params['rawfiles'][$rawfile]['biorep'],
				"orientation" => $this->iqpir_params['rawfiles'][$rawfile]['orientation'], "noise" => $next_frag_noise));
			}
			if(count($pepfraginfo_keys)===0) $pepfraginfo_keys = array_keys($data[$next_xl][0]);
			if($this->iqpir_params['normalize'] == "true" && $excl == "FALSE") {
				array_push($sample_biorepnorms[$this->iqpir_params['rawfiles'][$rawfile]['biorep']], $next_info['frag_log2ratios'][$k]);
				if(! $this->use_optimal_ratios) array_push($sample_optimal_biorepnorms[$this->iqpir_params['rawfiles'][$rawfile]['biorep']], $next_info['frag_optimal_ratios'][$k]);
			}
		}
		for($k = 0;  $k < count($next_larm_rep_ids); $k++) {
			$excl = "FALSE";
			if($this->min_sig2noise > 0 && $next_frag_noise > 0 && $next_larm_rep_ints[$k] < 5 * 
				$this->min_sig2noise * $next_frag_noise) {
				echo "Ignoring fragment ratio of {$next_larm_rep_ids[$k]} of {$next_xl} due to intensity {$next_larm_rep_ints[$k]} being less than ".(5*$this->min_sig2noise)." times noise {$next_frag_noise} (".
					($this->min_sig2noise * 5*$next_frag_noise).") in {$rawfile} scan {$next_larm_rep_scans[$k]}\n";
				$excl = "SIG2NOISE"; // avoid the ints
			}
			$next_larm_rep_ratios[$k] *= $multiplier;
			if(! $this->use_optimal_ratios) {
				array_push($data[$next_xl], array("rawfile" => $rawfile, "scan" => $next_larm_rep_scans[$k], 
				"id" => $next_larm_rep_ids[$k], "logratio" =>  $next_larm_rep_ratios[$k], "mzs" => "807.442528,811.455947", 
				"total_intensity" => $subtract_noise && $next_larm_rep_ints[$k] > $next_frag_noise[$k] * 5 ? 
				$next_larm_rep_ints[$k] - $next_frag_noise * 5 :$next_larm_rep_ints[$k], 
				"excluded" => $excl, "theor_relints" => "0.759,0.199,0.037,0.004,0.001", "obs_relints" => $next_larm_rep_ints[$k], "biorep" => $this->iqpir_params['rawfiles'][$rawfile]['biorep'],
				"orientation" => $this->iqpir_params['rawfiles'][$rawfile]['orientation'], 
				"optimal_logratio" => $next_larm_rep_ratios[$k], "noise" => $next_frag_noise));
			}
			else {
				array_push($data[$next_xl], array("rawfile" => $rawfile, "scan" => $next_larm_rep_scans[$k], 
				"id" => $next_larm_rep_ids[$k], "logratio" =>  $next_larm_rep_ratios[$k], "mzs" => "807.442528,811.455947", 
				"total_intensity" => $subtract_noise && $next_larm_rep_ints[$k] > $next_frag_noise * 5 ? 
				$next_larm_rep_ints[$k] - $next_frag_noise * 5 : $next_larm_rep_ints[$k], 
				"excluded" => $excl, "theor_relints" => "0.759,0.199,0.037,0.004,0.001", "obs_relints" =>  $next_larm_rep_ints[$k], "biorep" => $this->iqpir_params['rawfiles'][$rawfile]['biorep'],
				"orientation" => $this->iqpir_params['rawfiles'][$rawfile]['orientation'], "noise" => $next_frag_noise));
			}
			if($this->iqpir_params['normalize'] == "true" && $excl == "FALSE") {
				array_push($sample_biorepnorms[$this->iqpir_params['rawfiles'][$rawfile]['biorep']], $next_larm_rep_ratios[$k]);
				if(! $this->use_optimal_ratios) array_push($sample_optimal_biorepnorms[$this->iqpir_params['rawfiles'][$rawfile]['biorep']], $next_larm_rep_ratios[$k]);
			}
		
		}
	}
	fclose($myfile);
	echo "Read in info for ".count($data)." cross-links from file {$file} for sample {$sample}\n"; 

	if($this->iqpir_params['normalize'] == "true" && ! $this->use_current_normalization) {
        foreach($sample_biorepnorms as $biorep => $ratios) {
			if($this->normalization_method == "mode") $this->iqpir_params['samples'][$sample][$biorep] = array(number_format(-1 * MyStatistics::getDistrMode($ratios), 2, ".", ""));
			
			else $this->iqpir_params['samples'][$sample][$biorep] = array(number_format(-1 * MyStatistics::median($ratios), 2, ".", ""));
			if(! $this->use_optimal_ratios) array_push($this->iqpir_params['samples'][$sample][$biorep], 
				number_format(-1 * MyStatistics::median($sample_optimal_biorepnorms[$biorep]), 2, ".", ""));
			else array_push($this->iqpir_params['samples'][$sample][$biorep], $this->iqpir_params['samples'][$sample][$biorep][0]);
			echo "Have {$this->iqpir_params['samples'][$sample][$biorep][0]} norm for {$sample} sample and biorep {$biorep}\n";
			echo "Have optimal {$this->iqpir_params['samples'][$sample][$biorep][1]} norm for {$sample} sample and biorep {$biorep}\n";
		}
	}
	
	$pepFragMeans = array();
	$reporterMeans = array();
	$resPairMeans = array();
	$resPairPepAs = array();
	$resPairPepBs = array();
	$protPairMeans = array();
	$resHubMeans = array();
	$pepFragOptimalMeans = array();
	$pep_and_frag_ratios = array();
	$pep_and_frag_optimal_ratios = array();
	$next_optimal_norm = array();
	$fsample = fopen($samplefile, "w");
	$prot_suff = $this->deadend ? "" : "s";
	$res_suff = $this->deadend ? "" : "_pair";
	if($this->deadend) fwrite($fsample, "deadend-peptide");
	else fwrite($fsample, "cross-link");
	fwrite($fsample, "\txlinkdb-id\tprotein{$prot_suff}\tresidue{$res_suff}\tactual_ratio\tlogratio_mean\tlogratio_num\tlogratio_std\tlogratio_tstat\tlogratio_df\tlogatio_pval\tlogratio_outliers");
	if(! $this->use_optimal_ratios) fwrite($fsample, "\toptimal_logratio_mean\toptimal_num\toptimal_logratio_std\toptimal_logratio_tstat\toptimal_logratio_df\toptimal_logratio_pval\toptimal_logratio_outliers");
	fwrite($fsample, "\n");
	$fchrom = fopen($chromfile, "w");

	if($this->deadend) fwrite($fchrom, "deadend-peptide");
	else fwrite($fchrom, "cross-link");
	fwrite($fchrom, "\txlinkdb-id\tprotein{$prot_suff}\tresidue{$res_suff}\t".join("\t", $pepfraginfo_keys)."\tratio\tratio_error");
	if(! $this->use_optimal_ratios) fwrite($fchrom, "\toptimal_ratio\toptimal_ratio_error");
	fwrite($fchrom, "\tnorm_logratio\tnorm_ratio");
	if(! $this->use_optimal_ratios) fwrite($fchrom, "\tnorm_optiomal_logratio\tnorm_optimal_ratio");
	fwrite($fchrom, "\n");

	foreach($data as $xl => $info) {
		$ratios = array();
		$optimal_ratios = array();
		for($k = 0; $k < count($info); $k++) {
			if($info[$k]['excluded']!= "FALSE") continue;
			$norm = $this->iqpir_params['samples'][$sample][$info[$k]['biorep']];
			array_push($ratios, $info[$k]['logratio'] + $norm[0]);
			if(! $this->use_optimal_ratios) array_push($optimal_ratios, $info[$k]['optimal_logratio'] + $norm[1]);
		}
		if($include_reporters && array_key_exists($xl, $reporter_data)) {
			for($r = 0; $r < count($reporter_data[$xl]); $r++) {
				$norm = $this->iqpir_params['samples'][$sample][$reporter_data[$xl][$r][3]];
				array_push($ratios, $reporter_data[$xl][$r][0] + $norm[0]);
			}
		}
		$nextPeptideAndFragInfo2 = $this->getWeightedAverageRatio($ratios, NULL, true, 0, false, true, true);
		
		$next_actual = array_key_exists($xl_prots[$xl], $actual_ratios) ? $actual_ratios[$xl_prots[$xl]] : 0;

		$outliers = $nextPeptideAndFragInfo2[count($nextPeptideAndFragInfo2)-1] == "" ? array() : explode(",", $nextPeptideAndFragInfo2[count($nextPeptideAndFragInfo2)-1]);
		$next_norminfo_optimal = $this->getWeightedAverageRatio($optimal_ratios, NULL, true, 0, false, true, true); //$xl == "ADK[325.13]LAEEHGS_AMK[325.13]QVAGTMK");
		$outliers_optimal = $next_norminfo_optimal[count($next_norminfo_optimal)-1] == "" ? array() : explode(",", $next_norminfo_optimal[count($next_norminfo_optimal)-1]);
		for($k = 0; $k < count($outliers); $k++) {
			$found = false;
			$logratios = array();
			for($j = 0; $j < count($info); $j++) {
				if(abs($outliers[$k] - $this->iqpir_params['samples'][$sample][$info[$j]['biorep']][0] 
					- $info[$j]['logratio']) < 0.01) {
					$info[$j]['excluded'] = "OUTLIER";
					$found = true;
				}
				array_push($logratios, $info[$j]['logratio'] + $this->iqpir_params['samples'][$sample][$info[$j]['biorep']][0]);
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
				echo "Error: could not find outlier normalized value {$outliers[$k]} among logratios of {$xl}: ".join(",", $logratios)."\n";
				exit(1);
			}
		}
		$next_actual = array_key_exists($xl_prots[$xl], $actual_ratios) ? $actual_ratios[$xl_prots[$xl]] : "";
		$num_RH_SH_ids[$xl]['log2ratio'] = $nextPeptideAndFragInfo2[0];

		fwrite($fsample, $xl . "\t" . $xlinkdb_ids[$xl] . "\t" . $xl_prots[$xl] . "\t" . $xl_residuepairs[$xl] . "\t" . $next_actual . "\t" . join("\t", $nextPeptideAndFragInfo2));
		if(! $this->use_optimal_ratios) fwrite($fsample, "\t" . join("\t", $next_norminfo_optimal));
		
		fwrite($fsample, "\n");
		if($nextPeptideAndFragInfo2[1] === 0) continue; // nothing more to do

		$nextProts = explode("-", $xl_prots[$xl]);
		$nextres = explode("-", $xl_residuepairs[$xl]); 
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
$next_prots = $xl_prots[$xl];
$next_res = $xl_residuepairs[$xl];
		if(! $this->deadend && $xl_prots[$xl]!=$nextProts && array_key_exists($nextProts . "\t" . $nextres, $resPairMeans)) {
			$next_prots  = $nextProts;
			$next_res = $nextres;
			$indexA = 2;
			$indexB = 0; // peptides are in reverse order
		}
		else if(! array_key_exists($xl_prots[$xl] . "\t". $xl_residuepairs[$xl], $resPairMeans)) {
				$resPairMeans[$xl_prots[$xl] . "\t". $xl_residuepairs[$xl]] = array();
				$resPairPepAs[$xl_prots[$xl] . "\t". $xl_residuepairs[$xl]] = array();
				$resPairPepBs[$xl_prots[$xl] . "\t". $xl_residuepairs[$xl]] = array();
		}
		$next_id = $next_prots . "\t" . $next_res;
		array_push($resPairMeans[$next_id], $nextPeptideAndFragInfo2[0]);
		$nextPeps = explode("_", $xlinkdb_ids[$xl]);
		$resPairPepAs[$next_id][$nextPeps[$indexA]] = 1;
		if(! $this->deadend) $resPairPepBs[$next_id][$nextPeps[$indexB]] = 1;

		if(! $this->deadend && $nextres[0] != 0 && $nextres[1] != 0) { // make sure not res 0 meaning peptide not found in protein sequence
			if(! array_key_exists($nextProts, $protPairMeans)) $protPairMeans[$nextProts] = array();
			array_push($protPairMeans[$nextProts], $nextPeptideAndFragInfo2[0]);
		}
		for($k = 0; $k < count($info); $k++) {
			$norm = $this->iqpir_params['samples'][$sample][$info[$k]['biorep']];
			$next_th = explode(",", $info[$k]['theor_relints']);
			$next_ob = explode(",", $info[$k]['obs_relints']);
			$current_ratio = number_format(pow(2, $info[$k]['logratio']), 2, ".", "");
			$current_error = 0;
			$optimal_ratio = $current_ratio;
			$optimal_error = 0;
			if(strpos( $info[$k]['id'], "LRep") === false) {
				$current_error = $this->computeRatioDeconvError($next_th, $next_ob, $current_ratio);
				$optimal_ratio = $this->getOptimalRatio($next_th, $next_ob);
				$optimal_error = $this->computeRatioDeconvError($next_th, $next_ob, $optimal_ratio);
			}
			$norm_logratio = $info[$k]['logratio']+$norm[0];
			$norm_ratio = number_format(pow(2, $norm_logratio), 2, ".", "");
			if(! $this->use_optimal_ratios) {
				$norm_optimal_logratio = $info[$k]['optimal_logratio']+$norm[1];
				$norm_optimal_ratio = number_format(pow(2, $norm_optimal_logratio), 2, ".", "");
			}
			fwrite($fchrom, $xl . "\t" . $xlinkdb_ids[$xl] . "\t" . $xl_prots[$xl] . "\t". $xl_residuepairs[$xl] . "\t" . join("\t", array_values($info[$k])).
				"\t" . $current_ratio . "\t" . $current_error);
			if(! $this->use_optimal_ratios) fwrite($fchrom, "\t" . $optimal_ratio . "\t" . $optimal_error);
			fwrite($fchrom, "\t" . $norm_logratio . "\t" . $norm_ratio);
			if(! $this->use_optimal_ratios) fwrite($fchrom, "\t" . $norm_optimal_logratio . "\t" . $norm_optimal_ratio);
			fwrite($fchrom, "\n");
		}
	}
	fclose($fchrom);
	fclose($fsample);
	if(count($light_heavy_info) > 0) {
		$fid = fopen($numidsfile, "w");
		fwrite($fid, "cross-link\txlinkdb_id\tnum_RH_ids\tnum_SH_ids\tRH2SH_log2ratio\n");
		$num_rh = array();
		$num_sh = array();
		foreach($num_RH_SH_ids as $xl => $info) {
			fwrite($fid, $xl . "\t" . $xlinkdb_ids[$xl] . "\t" . $info['RH'] . "\t" . $info['SH'] . "\t" . $info['log2ratio'] . "\n");
			array_push($num_rh, $info['RH']);
			array_push($num_sh, $info['SH']);
		}
		fclose($fid);
		$stats_rh = MyStatistics::stats($num_rh, 1, 2);
		$stats_sh = MyStatistics::stats($num_sh, 1, 2);
		$text = $this->deadend ? "deadend peptide" : "cross-link";
		echo "File {$numidsfile} written for sample {$sample} with average {$stats_rh[0]} RH and {$stats_sh[0]} SH ids per {$text}\n";
	}
	$have_respairs = false;
	foreach($resPairMeans as $id => $ratios) {
		if(count($ratios) < 2) continue;
		if(! $have_respairs) {
			$fres = fopen($respairfile, "w");
			fwrite($fres, "protein{$prot_suff}\tresidue{$res_suff}\tmean_logratio\tstdev_logratio\tnumreps\t95conf\tpepA_seqs");
			if(! $this->deadend) fwrite($fres, "\tpepB_seqs");
			fwrite($fres, "\n");
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
	else echo "Files {$samplefile} and {$chromfile} written for sample {$sample}\n";
	$first = true;
	$add_zscore2ratiofile = true;
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
			fwrite($fprot, "\n");
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
				fwrite($fout, $line . "\tprotpair_z\n");
				$first = false;
			}
			else {
				$next = explode("\t", $line);
				$nextProts = explode("-", $next[2]);
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
		fclose($fout);
		rename($samplefile . ".tmp", $samplefile);
	} // if add z score
	$first = true;
	uksort($resHubMeans, array($this, "sortByFirstAndSecondIndex"));
	$hub_written = false;
	foreach($resHubMeans as $id => $info) {
		if($first) {
			$fres = fopen($reshubfile, "w");
			fwrite($fres, "protein\tresidue\tmean_logratio\tstdev_logratio\tnumreps\t95conf\tpartner_protein\tpartner_residue\tmean_logratio\n");
			$first = false;
			$hub_written = true;
			continue;
		}
		usort($info, array($this, "sortByThirdIndex")); 
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
	
} // next sample
	if($this->iqpir_params['normalize']=="true" && ! $this->use_current_normalization) {
		// append to previous file, if nec		
		if(count($samples) > 0 && file_exists($normfile)) {
			$first = true;
			$myfile = fopen($normfile, "r");
			while($line = rtrim(fgets($myfile))) {
				if($first) {
					$first = false;
					continue;
				}
				$next = explode("\t", $line);
				if(! array_key_exists($next[0], $samples)) {
					$this->iqpir_params['samples'][$next[0]][$next[1]] = array($next[2], $next[2]);
					echo "Adding pre-existing normalization additive factor {$next[2]} for sample {$next[0]} biorep {$next[1]}\n";
				}
			}
			fclose($myfile);
		}
		$fnorm = fopen($normfile, "w");
		fwrite($fnorm, "sample\tbiorep\tnorm_add_to_log\n");
		foreach($this->iqpir_params['samples'] as $sample => $bioreps) {
			foreach($bioreps as $biorep => $norm) {
				fwrite($fnorm, $sample . "\t" . $biorep . "\t" . number_format($norm[0], 3, ".", "")."\n");
			}
		}
		fclose($fnorm);	
		echo "Normalization file {$normfile} written for all samples: ".join(",", array_keys($this->iqpir_params['samples']))."\n";
	}
	$this->generateIonRatioPlot();

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

public function generateIonRatioPlot() {
	$chromfiles = glob($this->iqpir_params['iqpir_output_dir'] . "sample.*.chroms.txt");
	$outfile = $this->iqpir_params['iqpir_output_dir'] . "ionratio_distrs.txt";
	if(count($chromfiles) == 0) return;
	$data = array();
	$min_val = 999;
	$max_val = -1;
	$offset = 1; // 8
	$leg_length = 0;
	for($k = 0; $k < count($chromfiles); $k++) {
		$next_id = "";
		$nextpos = strrpos($chromfiles[$k], "/");
		if($nextpos !==false) {
			$next_id = substr($chromfiles[$k], $nextpos+$offset);
		}
		else $next_id = substr($chromfiles[$k], $offset);
		$next_id = substr($next_id, 0, strlen($next_id) - 11);
		if(strlen($next_id) > $leg_length) $leg_length = strlen($next_id);
		$data[$next_id] = $this->getIonRatioDistr($chromfiles[$k]);
		if($data[$next_id][2] < $min_val) $min_val = $data[$next_id][2];
		if($data[$next_id][3] > $max_val) $max_val = $data[$next_id][3];
	}
	$fout = fopen($outfile, "w");
	fwrite($fout,  "log2ratio");
	foreach($data as $id => $ratios) fwrite($fout, "\t" . $id);
	fwrite($fout, "\n");
	
    for($k = $min_val; $k <= $max_val; $k++) {
          fwrite($fout, ($k*0.2 -5));
          foreach($data as $id => $ratios) {
          		fwrite($fout, "\t");
				if(array_key_exists($k, $ratios[0])) fwrite($fout, number_format($ratios[0][$k]/$ratios[1], 3, ".", ""));
				else fwrite($fout, "0");
          	}
          	fwrite($fout, "\n");
	}
	fclose($fout);
	$width = 500 + 9.3 * $leg_length;
	$par_right = 20; //2 * $leg_length;
	$par_right = 5 * $leg_length/22 + 13.2;
	$scriptfile = $this->iqpir_params['iqpir_output_dir'] . "rscript.R";
	$pngfile = $this->iqpir_params['iqpir_output_dir'] . "ionratio_distrs.png";
	$colors = array("blue","red","green","orange","brown","pink","black","grey");
	$color_ind = 0;
	$linetype = 1;
	$col_colors = array();
	$col_linetypes = array();
	$fr = fopen($scriptfile, "w");
	fwrite($fr, "mydf <- read.table('{$outfile}', sep = '\\t', header = T)\n");
	fwrite($fr, "png('{$pngfile}', width = {$width}, height = 450)\n");
	fwrite($fr, "par(mar = c(5,4,4,{$par_right}), xpd=TRUE)\n");
	$samples = array_keys($data);
	for($k = 0; $k < count($samples); $k++) {
		if($color_ind >= count($samples)) {
			$color_ind = 0;
			$linetype++;
		}
		array_push($col_colors, $colors[$color_ind]);
		array_push($col_linetypes, $linetype);
		if($k == 0) fwrite($fr, "plot(mydf[, c(1)], mydf[, c(".($k+2).")], pch = 19, col = '{$colors[$color_ind]}', type='l', main='Ion Log2ratio Distributions', xlab='Log2ratio', ylab='Fraction')\n");
		else fwrite($fr, "points(mydf[, c(1)], mydf[, c(".($k+2).")], col = '{$colors[$color_ind]}', type='l')\n");
		$color_ind++;
	}
	fwrite($fr, "legend('topright', legend = c('".join("','", $samples)."'),  lwd = 3, lty = c(".join(",", $col_linetypes).
		"), col = c('".join("','", $col_colors)."'), inset=c(-0.5,0))\n");
	fwrite($fr, "dev.off()\n");
	fclose($fr);
	
	$errfile = $this->iqpir_params['iqpir_output_dir'] . "rscript.err";
	$valid = array();
	system("Rscript {$scriptfile} >& {$errfile}", $valid);

	if(file_exists($pngfile)) {
		if(file_exists($scriptfile)) unlink($scriptfile);
		if(file_exists($errfile)) unlink($errfile);
	}
	echo "Ion Ratio Distributions PGN written to {$pngfile}\n";
}
public function getIonRatioDistr($chromsfile) {
$tot = 0;
$data = array();
$first = true;
$myfile = fopen($chromsfile, "r");
while($line = rtrim(fgets($myfile))) {
	if($first) {
		$first = false;
	}
	else {
		$next = explode("\t", $line);
		if($next[10] != "FALSE") continue;
		$val = number_format(($next[18] +5)/0.2, 0, ".", "");
		if($val < 0) $val = 0;
		if($val > 50) $val = 50;
		if(! array_key_exists($val, $data)) $data[$val] = 0;
		$data[$val]++;
		$tot++;
	}
}
fclose($myfile);
ksort($data);
$vals = array_keys($data);
return array($data, $tot, $vals[0], $vals[count($vals)-1]);
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
	public function getSampleNonIdents() {
		$xl_ids = array();
		$outfile = $this->iqpir_params['iqpir_output_dir'] . "nonidents.txt";
		$fout = fopen($outfile, "w");
		foreach($this->iqpir_params['samples'] as $key => $value) {
			$samplefile = $this->iqpir_params['iqpir_output_dir'] . "sample." . $key . ".ratios.txt";
			if(! file_exists($samplefile)) continue;
			$myfile = fopen($samplefile, "r");
			$first = true;
			while($line = rtrim(fgets($myfile))) {
				if($first) {
					$first = false;
					continue;
				}
				$next = explode("\t", $line);
				if(! array_key_exists($next[1] . "\t" . $next[2], $xl_ids)) $xl_ids[$next[1] . "\t" . $next[2]] = $this->getSampleArray();
				unset($xl_ids[$next[1] . "\t" . $next[2]][$key]);
			
			}	
			fclose($myfile);
		} // next sample
		ksort($xl_ids);
		if(! $this->deadend) fwrite($fout, "cross-link\tprotein_pair");
		else fwrite($fout, "deadend_peptide\tprotein");
		fwrite($fout, "\tnon-ident_samples_".join("-", array_keys($this->getSampleArray()))."\n");
		foreach($xl_ids as $xl => $deficits) {
			if(count($deficits)===0) continue;
			fwrite($fout,  $xl . "\t" . join(",", array_keys($deficits)) . "\n");
		}
		fclose($fout);
		echo "Non-ident information written to {$outfile}\n";
	}

public function quantifyDeadends($spec_filter = "", $deconvolute = false) {
if(! $this->deadend) {
	echo "Error: quantifyDeadends only posisble in deadend mode.\n";
	exit(1);
}
function log2ratio($val, $log2bound = "") { //10) {
	if($val === "") return "";
	else if($val === 0) return $log2bound === "" ? $log2bound : -1 * $log2bound;
	else if($val == "inf") return $log2bound;
	
	return number_format(log($val)/log(2), 2, ".", "");;
}
// get the modmass info from pepprophet file, if present
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

	$file = $this->iqpir_params['deadendprophetfile']; 
	if(! file_exists($file)) {
		echo "Error: deadend pepXML file {$file} not found\n";
		exit(1);
	}
	echo "reading {$file}\n"; 
	$fdr = $this->iqpir_params['fdr'];
	$level = "probability";
	$heavy_mod = $this->iqpir_params['heavy_modmass'];
	$output_dir = $this->iqpir_params['iqpir_output_dir'];

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
		$spec_filter = explode(",", $spec_filter);
	}
	else $spec_filter = array();
	$outfile = $this->scan_ratiofile;

	$find_residuals = false; // wheterh to look for peaks with expected h/l spacing
	$pursue_all_pep_charges = true;
	
	$report_log2ratios = true;
	$protein_seqs_and_genes = array();

	$ratio_prefix = $report_log2ratios ? "log2" : "";

	echo "Ready to write to {$outfile} with fdr {$fdr}, filter level {$level}, and heavy mod {$heavy_mod}\n"; 
	$fout = fopen($outfile, "w");
	if($deconvolute) {
		fwrite($fout, "scan\tspectrum\tpeptide1\tprotein1\tpos1\tpeptide2\tprotein2\tpos2\tnoise\tisoints1\tlight1-1\tlight1-2\tlight1-3\theavy1-1\theavy1-2\theavy1-3\tdIntens1\t{$pepRatio}ratio1\t{$pepRatio}ratio1Error");
		fwrite($fout, "\tisoints2\tlight2-1\tlight2-2\tlight2-3\theavy2-1\theavy2-2\theavy2-3\tdIintens2\t{$pepRatio}ratio2\t{$pepRatio}ratio2Error");
		fwrite($fout, "\tisointsreporter\tlightreporter-1\tlightreporter-2\tlightreporter-3\theavyreporter-1\theavyreporter-2\theavyreporter-3\t{$reporterRatio}ratioreporter-1\t{$reporterRatio}ratioreporter-2");
		fwrite($fout, "\tPepA_Int1\tPepA_Int2\tPepB_Int1\tPepB_Int2\n");
	}
	else {
	
	
		fwrite($fout, "scan\tspectrum\tactual_ratio\tpeptide\tprotein\tpos\tnoise\tlight-1\tlight-2\tlight-3\theavy-1\theavy-2\theavy-3\t{$pepRatio}ratio1-1\t{$pepRatio}ratio1-2\t{$pepRatio}ratio1-3\tdeconv_h-l_ratio");
		if($find_residuals || $pursue_all_pep_charges) {
			for($k = 1; $k <= 3; $k++) {
				if($find_residuals) {
					fwrite($fout, "\tresiduals\tresidual_sum");
				}
				if($pursue_all_pep_charges) {
					fwrite($fout, "\tdeconv_h-l_ratio\tl-h_mz\tints{$k}_1-1\tintratio{$k}_1-1\tlightratio{$k}_1-1\tlight{$k}_1-1\tlight{$k}_1-2\tlight{$k}_1-3\theavy{$k}_1-1\theavy{$k}_1-2\theavy{$k}_1-3\t{$pepRatio}ratio{$k}_1-1\t{$pepRatio}ratio{$k}_1-2\t{$pepRatio}ratio{$k}_1-3");
				}

			}
		}
		fwrite($fout, "\tlightreporter-1\tlightreporter-2\tlightreporter-3\theavyreporter-1\theavyreporter-2\theavyreporter-3\t{$reporterRatio}ratioreporter-1\t{$reporterRatio}ratioreporter-2\t{$reporterRatio}ratioreporter-3");
	}
	fwrite($fout, "\tactual_{$ratio_prefix}ratio");
	fwrite($fout, "\treporter_ids\treporter_log2ratios\treporter_ints\tpeptide_ids\tpeptide_log2ratios\tpeptide_ints");
	fwrite($fout, "\tpeptide_mzs\tpeptide_theor_relints\tpetpide_obs_relints");
	fwrite($fout, "\tfrag_ids\tfrag_log2ratios\tfrag_ints");
	fwrite($fout, "\tfrag_mzs\tfrag_theor_relints\tfrag_obs_relints");
	if(! $this->use_optimal_ratios) {			
		fwrite($fout, "\toptimal_peptide_log2ratios\toptimal_frag_log2ratios");
	}	
	fwrite($fout, "\tresidue");
	fwrite($fout, "\tprotein");
	fwrite($fout, "\tnoise\n");

	$existing_iqpirquant = array();
	if(array_key_exists("existing_deadend_iqpirfiles", $this->iqpir_params) && $this->iqpir_params["existing_deadend_iqpirfiles"] != "") {
		$existing_iqpirfiles = explode(",", $this->iqpir_params["existing_deadend_iqpirfiles"]); //"/net/gs/vol4/shared/brucelab/search/xiaoting/Spectrast_test/interact_spectrast-xl-iqpir.xls", "/net/gs/vol4/shared/brucelab/search/xiaoting/iqPIR_mito_Tac_sham_mango_6pairs_135runs/iprophet-110920-xl-iqpir.xls");
		for($k = 0; $k < count($existing_iqpirfiles); $k++) {
			if(strpos($existing_iqpirfiles[$k], "-iqpir.xls") !== strlen($existing_iqpirfiles[$k]) - strlen("-iqpir.xls")) {
				echo "Error: existing iqpirfile {$existing_iqpirfiles[$k]} is not of correct format ending in -iqpir.xls\n";
				exit(1);
			}
			else if(! file_exists($existing_iqpirfiles[$k])) {
				echo "Error: existing iqpirfile {$existing_iqpirfiles[$k]} is not found\n";
				exit(1);
			}
		}
		$existing_iqpirquant = $this->reuseIqpirfiles($existing_iqpirfiles);
		echo "Read in existing quantitation for ".count($existing_iqpirquant)." spectral deadend scans\n"; 
	}




	$pepxml = new MyDeadendPepXML($file, $fdr, array("modmass_info" => $modmass_info));
	$scan = 0;
	$pep = "";
	$mass = ""; 
	$charge = "";
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
	$analyze = false;
	$is_decoy = false;
	while($pepxml->nextLine()) {
		if($pepxml->hasStartTag("msms_run_summary")) { // end of searh result
			$mzxml = $pepxml->getTagValue("base_name") . ".mzXML";
			$noisefile = $pepxml->getTagValue("base_name") . ".noise";
			if(! file_exists($mzxml)) {
				echo "Error: mzxml {$mzxml} does not exist\n";
				exit(1);
			}
			if(! file_exists($noisefile)) {
				echo "Error: noisefile {$noisefile} does not exist\n";
				exit(1);
			}
		}
		else if($pepxml->hasStartTag("spectrum_query")) { // end of searh result
			$spectrum = $pepxml->getTagValue("spectrum");
			$charge = $pepxml->getTagValue("assumed_charge");
			$scan = $pepxml->getTagValue("start_scan");
			$analyze = true;
		}
		else if($analyze && $pepxml->hasStartTag("search_hit")) { // end of searh result
			$mass = $pepxml->getTagValue("calc_neutral_pep_mass");
		}
		else if($analyze && $pepxml->hasEndTag("search_hit")) { // end of searh result
			$analyze = false;
		}
		else if(count($protein_seqs_and_genes)==0 && $pepxml->hasStartTag("search_database")) { // end of searh result
			$protein_seqs_and_genes = $this->getDatabaseProteinSequencesAndGenes($pepxml->getTagValue("local_path"));
		}
		else if($analyze && $pepxml->hasStartTag("modification_info")) { // end of searh result
			$pep = $pepxml->getTagValue("modified_peptide");
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
			$next_unique = $spec_without_ch . ":" . $this->convertHeavyToLight($pep);
			if($pepxml->aboveFilterThreshold() && array_key_exists($next_unique, $seen)) {
				$heavy = strpos($pep, $heavy_mod)!==false;
				$next_index = $heavy ? 1 : 0;
				$seen[$next_unique][$next_index]++; 
			}
			else if($pepxml->aboveFilterThreshold() && ! array_key_exists($next_unique, $seen)) {
				$seen[$next_unique] = array(0, 0);
			
				$noise = 0;
				$noise_valid = array();
				exec("grep -P '^".$scan."\t' {$noisefile}", $noise_valid);
				if(count($noise_valid)===1) {
					$next_noise = explode("\t", $noise_valid[0]);
					$noise = $next_noise[3];
				}
				$heavy = strpos($pep, $heavy_mod)!==false;
				$next_index = $heavy ? 1 : 0;
				$seen[$next_unique][$next_index]++; 
				if(count($existing_iqpirquant) > 0 && array_key_exists($next_unique, $existing_iqpirquant)) {
					fwrite($fout, $existing_iqpirquant[$next_unique] . "\n");
					$peps = array();
					$masses = array(); 
					$charges = array();
					$next_scans = array();
					$next_noise = array();
					$ms3scans = array();
					continue;
				}

				$num_isotopes = 5;
				$dlight1 = NULL;
				$proton_mass = 1.00727647;
				$min_mz = 500;
				$max_mz = 2000;
				

							$modpep = $pepxml->getModifiedPeptide(); 
							if($heavy) { // must change sequence to light version
								$modpep = $this->convertHeavyToLight($modpep); 
							}
							
							
				$pep1_stumpfrags = $this->getStumpFragments($modpep, 325.13); 
				$frag_ratios = array();
				$frag_ints = array();
				$seen_fragmasses = array(); // make sure only use once
				$fragment_relints = array();
				$pep1_frag_ratios = array();
				$pep1_frag_ints = array();
				$frag_fragments = array();
				$frag_mzs = array();
				$react_ms3 = true;
				$ms3scans = $react_ms3 ? array($scan + 2) : array();
				$next_scans = $analysis_type == "ReACT" && $react_ms3 ? array($scan + 2) : array($scan); //, $scan);
				$next_noise = array();

				$longarm_scans = array_key_exists('longarm_reporters', $this->iqpir_params) && 
					$this->iqpir_params['longarm_reporters'] == "true" && $analysis_type == "ReACT" ?
					array($scan + 4) : array();

				for($n = 0; $n < count($next_scans); $n++) {
					if($next_scans[$n] == $scan) array_push($next_noise, $noise);
					else {
						$noise_valid = array();
						exec("grep -P '^".$next_scans[$n]."\t' {$noisefile}", $noise_valid);
						if(count($noise_valid)===1) {
							$scan_noise = explode("\t", $noise_valid[0]);
							array_push($next_noise, max(1, $scan_noise[3]));
						}
						else array_push($next_noise, 1); // don't know
					}
				}				
				$this->ms2_scan = $scan;
				$this->setScanPeaklist($mzxml, $ms3scans, $longarm_scans);
				$fragment_theor_ints = array();
				$peptide_theor_ints = array();
				$peptide_relints = array();
				$frag_scans = array();
				$frag_noise = array();
				
				$min_peptide_massdiff = 0.5; // daltons?
				foreach($pep1_stumpfrags as $key => $value) {	
					$next_mass = strval($value[0]);
					$nextpos = strpos($next_mass, ".");
					if($nextpos!==false) $next_mass = substr($next_mass, 0, $nextpos+4); // keep only first 2 digits
					if(array_key_exists($next_mass, $seen_fragmasses)) continue;
					$seen_fragmasses[$next_mass] = 1;
					
					$light_fragentries[$key] = $this->getIsotopeQuant($value[0], 1, $next_scans[0], $mzxml, $noise, "frag_light{$key}_1 ");
					$heavy_fragentries[$key] = $this->getIsotopeQuant($value[0] + $pep_massdiff, 1, $next_scans[0], $mzxml, $noise, "frag_heavy{$key}_1 ");
					// now compute the ratios
					$next_ints_fragentries = $this->getIsotopeIntensities($value[1], true);
					if($next_ints_fragentries[0]==0) {
						echo "1. No theoreteical intensities obtained with calcisotopes for fragment ion \"{$value[1]}\": ".join(",", $next_ints_fragentries).", possibly due to unexpected peptide modification (contact Jimmy).  Will instead pass mass {$valid[0]}\n";
						$next_ints_fragentries = $this->getIsotopeIntensities($value[0], true);
					}
					if($this->use_optimal_ratios) {
						$next_fragdeconv = "N/A";
						$next_fragintens = $light_fragentries[$key][0][0] + $light_fragentries[$key][0][1] + $light_fragentries[$key][0][2] + $heavy_fragentries[$key][0][1] + $heavy_fragentries[$key][0][2];
						if($next_fragintens > 0) {
							$next_fragrelints = array(number_format($light_fragentries[$key][0][0]/$next_fragintens, 3, ".", ""), number_format($light_fragentries[$key][0][1]/$next_fragintens, 3, ".", ""), number_format($light_fragentries[$key][0][2]/$next_fragintens, 3, ".", ""), 
								number_format($heavy_fragentries[$key][0][1]/$next_fragintens, 3, ".", ""), number_format($heavy_fragentries[$key][0][2]/$next_fragintens, 3, ".", ""));
						
							$next_fragdeconv = $this->getOptimalRatio($next_ints_fragentries, $next_fragrelints);
						}
					}
					else {
						$next_fragdeconv = $this->getPeptideHeavy2LightRatioWithTheorInts($light_fragentries[$key], $heavy_fragentries[$key], $next_ints_fragentries);
					}
					if($next_fragdeconv != "N/A" && $next_fragdeconv !== "") {
							$next_fragintens = $light_fragentries[$key][0][0] + $light_fragentries[$key][0][1] + $light_fragentries[$key][0][2] + $heavy_fragentries[$key][0][1] + $heavy_fragentries[$key][0][2];
							if($next_fragintens > 0) {
								array_push($frag_ints, $next_fragintens);
								array_push($frag_ratios, $next_fragdeconv);
								array_push($pep1_frag_ints, $next_fragintens);
								array_push($pep1_frag_ratios, $next_fragdeconv);
								array_push($frag_fragments, "A:" . $key . "_1");
							
								array_push($fragment_theor_ints, join(",", $next_ints_fragentries));
		
								$next_fragrelints = array(number_format($light_fragentries[$key][0][0]/$next_fragintens, 3, ".", ""), number_format($light_fragentries[$key][0][1]/$next_fragintens, 3, ".", ""), number_format($light_fragentries[$key][0][2]/$next_fragintens, 3, ".", ""), 
									number_format($heavy_fragentries[$key][0][1]/$next_fragintens, 3, ".", ""), number_format($heavy_fragentries[$key][0][2]/$next_fragintens, 3, ".", ""));
								array_push($fragment_relints, join(",", $next_fragrelints));
								array_push($frag_mzs, ($value[0] + $proton_mass) . "-" . ($value[0] + $pep_massdiff + $proton_mass));
								array_push($frag_scans, $next_scans[0]);
								array_push($frag_noise, $next_noise[0]);
							}
					}
				}
				if(count($frag_ratios) > 0) {
					$nextFragInfo = $this->getWeightAverageWithRatio($frag_ratios, $frag_ints, NULL, NULL, $report_log2ratios);
				}
				$light_entries1 = $pursue_all_pep_charges ? array("1" => array(), "2" => array(), "3" => array()) : array(); // charges 1,2,3
				$heavy_entries1 = array("1" => array(), "2" => array(), "3" => array()); // charges 1,2,3
				$dlight1_entries1 = array("1" => array(), "2" => array(), "3" => array()); // charges 1,2,3
				$ratios_entries1 = array("1" => array(), "2" => array(), "3" => array()); // charges 1,2,3
				$ints_entries1 = array("1" => array(), "2" => array(), "3" => array()); // charges 1,2,3
				$deconv_ratio_entries1 = array("1" => array(), "2" => array(), "3" => array()); // charges 1,2,3
				$mzs1 = array();  // from charge to array of light/heavy
				$ppmtolerance = $this->iqpir_ppmtolerance; //25;
				$next_hk_info = $find_residuals ? $this->chargeStateDeconvolute($mzxml, $ppmtolerance, $pep_massdiff, $scan) : array();
				foreach($next_hk_info as $key => $value) {
					if($key > $masses[0] && $key > $masses[1]) unset($next_hk_info[$key]);
					else if(abs($key - $mass) <= 0.2) unset($next_hk_info[$key]);
				}				
				foreach($light_entries1 as $key => $value) {
					$next_charge = intval($key);
					$symbol = $next_charge == $charge ? "*" : "";
					if($heavy) {
						if($mass/$next_charge >= $min_mz && $mass/$next_charge <= $max_mz) {
							$light_entries1[$key] = $this->getIsotopeQuant($mass - $pep_massdiff, $next_charge, $scan, $mzxml, $noise, "{$symbol}light{$key}_1 ");
							$heavy_entries1[$key] = $this->getIsotopeQuant($mass, $next_charge, $scan, $mzxml, $noise, "{$symbol}heavy{$key}_1 ");
							if($deconvolute) $dlight1_entries1[$key] = $this->getIsotopeQuant($mass - $pep_massdiff, $next_charge, $scan, $mzxml, $noise, "light{$key}_1 ", $num_isotopes);
							
							$mzs1[$key] = array(($mass - $pep_massdiff + $next_charge * $proton_mass)/$next_charge, ($mass + $next_charge * $proton_mass)/$next_charge);
						}
					}
					else {
						if(($mass+$pep_massdiff)/$next_charge >= $min_mz && ($mass+$pep_massdiff)/$next_charge <= $max_mz) {
							$light_entries1[$key] = $this->getIsotopeQuant($mass, $next_charge, $scan, $mzxml, $noise, "light{$key}_1 ");
							$heavy_entries1[$key] = $this->getIsotopeQuant($mass + $pep_massdiff , $next_charge, $scan, $mzxml, $noise, "{$symbol}heavy{$key}_1 ");
							if($deconvolute) $dlight1_entries1[$key] = $this->getIsotopeQuant($mass, $next_charge, $scan, $mzxml, $noise, "{$symbol}light{$key}_1 ", $num_isotopes);
							$mzs1[$key] = array(($mass + $next_charge * $proton_mass)/$next_charge, ($mass + $pep_massdiff + $next_charge * $proton_mass)/$next_charge);
						}
					}
					if(count($light_entries1[$key]) > 0 && count($heavy_entries1[$key]) > 0) {
						$ratios_entries1[$key] = $this->printRatios($light_entries1[$key][0], $heavy_entries1[$key][0], "Scan {$scan} Pep1 {$key}_Ratios: ");

						if($deconvolute) {
							$modpeps = $pepxml->getModifiedPeptidePair();
							if($heavy) { // must change sequence to light version
								$modpeps[0] = $this->convertHeavyToLight($modpeps[0]); //preg_replace("/K\[327.13\]/", "K[325.13]", $modpeps[0]);
								$modpeps[1] = $this->convertHeavyToLight($modpeps[1]); //preg_replace("/K\[327.13\]/", "K[325.13]", $modpeps[1]);
							}

							$deconv_ratio_entries1[$key] = $this->computeRatioFromIsotopeSetsOfPeptide($dlight1_entries1[$key][0],  $modpeps[0]);

						}
				
						$ints1 = array();
						if($deconvolute) {
							$ratios_entries1[$key] = $this->computeRatioFromIsotopeSets($light_entries1[$key][0], $heavy_entries1[$key][0], 2, $modpeps[0]);
							$ints_entries1[$key] = array_pop($ratios_entries1[$key]);
							
						}
						else {
							$ints_entries1[$key] = $this->getIsotopeIntensities($modpep, true);
							if($ints_entries1[$key][0]==0) {
								echo "3. No theoreteical intensities obtained with calcisotopes for peptide ion \"{$modpep}\": ".
								join(",", $ints_entries1[$key]).", possibly due to unexpected peptide modification (contact Jimmy).  
								Will instead pass mass {$mass}\n";
								$ints_entries1[$key] = $this->getIsotopeIntensities($mass, true);
							}
						}
					}
				}
				
				if($heavy) {
					$light1ints = $this->getIsotopeQuant($mass - $pep_massdiff, $charge, $scan, $mzxml, $noise, "light1 ");
					$heavy1ints = $this->getIsotopeQuant($mass, $charge, $scan, $mzxml, $noise, "heavy1 ");
					if($deconvolute) $dlight1 = $this->getIsotopeQuant($mass - $pep_massdiff, $charge, $scan, $mzxml, $noise, "light1 ", $num_isotopes);
				}
				else {
					$light1ints = $this->getIsotopeQuant($mass, $charge, $scan, $mzxml, $noise, "light1 ");
					$heavy1ints = $this->getIsotopeQuant($mass + $pep_massdiff , $charge, $scan, $mzxml, $noise, "heavy1 ");
					if($deconvolute) $dlight1 = $this->getIsotopeQuant($mass, $charge, $scan, $mzxml, $noise, "light1 ", $num_isotopes);
				}
				$ratios1 = $this->printRatios($light1ints[0], $heavy1ints[0], "Scan {$scan} Pep1 Ratios: ");
					if($deconvolute) {
						$modpeps = $pepxml->getModifiedPeptidePair();
						if($heavy) { // must change sequence to light version
							$modpeps[0] = $this->convertHeavyToLight($modpeps[0]); //preg_replace("/K\[327.13\]/", "K[325.13]", $modpeps[0]);
							$modpeps[1] = $this->convertHeavyToLight($modpeps[1]); //preg_replace("/K\[327.13\]/", "K[325.13]", $modpeps[1]);
						}
					}
				$ints1 = array();
				if($deconvolute) {
					$ratios1 = $this->computeRatioFromIsotopeSets($light1ints[0], $heavy1ints[0], 2, $modpeps[0]);
					$ints1 = array_pop($ratios1);
				}
				$lightreporter = $this->getIsotopeQuant(807.442528, 1, $scan, $mzxml, $noise, "lightreporter ");
				$heavyreporter = $this->getIsotopeQuant(811.455947, 1, $scan, $mzxml, $noise, "heavyreporeter ");
				$ratios3 = $this->printRatios($heavyreporter[0], $lightreporter[0], "Scan {$scan} Reporter Ratios: ");
				$next_ints = $this->getIsotopeIntensities(807.442528, true);

				$longarm_ratios = array();
				$longarm_reporter_ints = array();
				$longarm_reporter_ratios = array();
				$longarm_reporter_ids = array();
				$seen_longarm_scans = array();
				for($l = 0; $l < count($longarm_scans); $l++) {
					$next_lightreporter = $this->getIsotopeQuant(807.442528, 1, $longarm_scans[$l], $mzxml, $noise, "lightreporter ");
					$next_heavyreporter = $this->getIsotopeQuant(811.455947, 1, $longarm_scans[$l], $mzxml, $noise, "heavyreporeter ");
					$next_lightreporter[0] = array(array_sum($next_lightreporter[0]));
					$next_heavyreporter[0] = array(array_sum($next_heavyreporter[0]));
					$next_ratios = $this->printRatios($next_heavyreporter[0], $next_lightreporter[0], "Scan {$scan} Reporter Ratios: ");
					if($next_ratios[0] != "N/A") {
						array_push($longarm_reporter_ints, $next_lightreporter[0][0] + $next_heavyreporter[0][0]);
						array_push($longarm_reporter_ratios, $next_ratios[0]);
						array_push($longarm_reporter_ids, ($l == 0 ? "A" : "B") . ":LRep");
						array_push($seen_longarm_scans, $longarm_scans[$l]);
					}				

				}

				$heavy_error = $this->calculateIntensErrorVersusTheor($heavyreporter[0], $next_ints, 3);
				$light_error = $this->calculateIntensErrorVersusTheor($lightreporter[0], $next_ints, 3);
				$max_error = 0.045;
				
				$success = $heavy_error <= $max_error && $light_error <= $max_error;
				if(! $success) $ratios3 = array("N/A", "N/A", "N/A"); // cancel it
				
				if($deconvolute) $ratios3 = $this->computeRatioFromIsotopeSets($lightreporter[0], $heavyreporter[0], 4, 807.442528);
				$ints3 = array();
				if($deconvolute) {
					$ratios3 = $this->computeRatioFromIsotopeSets($lightreporter[0], $heavyreporter[0], 4, 807.442528);
					$ints3 = array_pop($ratios3);
				}
				$next_actual_ratio = 1;
				if(strpos($spectrum, "_2Sto1R_")!==false) {
					$next_actual_ratio = 0.5;
				}
				else if(strpos($spectrum, "_2Rto1S_")!==false) {
					$next_actual_ratio = 2;
				}
	
				$peptide_pair = $pepxml->getPreservedOrderedModifiedPeptide();
				fwrite($fout, $scan . "\t" . $spectrum . "\t" . $next_actual_ratio . "\t" . join("\t", $peptide_pair) ."\t" . $noise . "\t");
				if($deconvolute) fwrite($fout, join(",", $ints1)."\t"); 
				fwrite($fout, number_format($light1ints[0][0], 0) . "\t" . number_format($light1ints[0][1], 0) . "\t" . number_format($light1ints[0][2], 0) . "\t");
				fwrite($fout, number_format($heavy1ints[0][0], 0) . "\t" . number_format($heavy1ints[0][1], 0) . "\t" . number_format($heavy1ints[0][2], 0) . "\t");
				if($deconvolute) {
						fwrite($fout, join("\t", $deconv_ratio1) . "\t");
				}
				else {
					fwrite($fout, join("\t", $ratios1) . "\t");
				}
				$next_deconv = $this->getPeptideHeavy2LightRatioWithTheorInts($light1ints, $heavy1ints, $this->getIsotopeIntensities($modpep, true));
				
				
				fwrite($fout, $next_deconv . "\t");
	
				// compute the residuals
				$peptide_ints = array();
				$peptide_ratios = array();
				$peptide_ids = array();
				$reporter_ints = array();
				$reporter_ratios = array();
				$reporter_ids = array();
				$first = true;
				$peptide_mzs = array();
				foreach($light_entries1 as $key => $value) {
						$next_charge = intval($key);
						$symbol = $next_charge == $charge ? "*" : "";
						fwrite($fout, "--{$symbol}{$key}-pep-->");
						$first = false;
					if(count($light_entries1[$key]) > 0 && count($heavy_entries1[$key]) > 0) {
						if($find_residuals) {
							$next_resid = array($light_entries1[$key][0][0], $light_entries1[$key][0][1], $light_entries1[$key][0][2], $heavy_entries1[$key][0][1], $heavy_entries1[$key][0][2]);
							$int_sum = array_sum($next_resid);

							if($light_entries1[$key][0][0] != "" && $light_entries1[$key][0][0] != 0 && $int_sum > 0) {
								for($r = 0; $r < 5; $r++) {
									$next_resid[$r] /= $int_sum;
								}
								$next_resid = array(($ints_entries1[$key][0]-$next_resid[0]) * ($ints_entries1[$key][0]-$next_resid[0]),
											($ints_entries1[$key][1]-$next_resid[1]) * ($ints_entries1[$key][1]-$next_resid[1]),
											($ints_entries1[$key][2]-$next_resid[2]) * ($ints_entries1[$key][2]-$next_resid[2]),
											($ints_entries1[$key][3]-$next_resid[3]) * ($ints_entries1[$key][3]-$next_resid[3]),
											($ints_entries1[$key][4]-$next_resid[4]) * ($ints_entries1[$key][4]-$next_resid[4]));		
								fwrite($fout, join(",", $next_resid) . "\t" . number_format(sqrt(array_sum($next_resid))/5, 3) . "\t");	
							}	
							else {
								fwrite($fout, "\t\t");
							}	
						}
						if($this->use_optimal_ratios) {
							$next_deconv = "N/A";
							$next_intens = $light_entries1[$key][0][0] + $light_entries1[$key][0][1] + $light_entries1[$key][0][2] + $heavy_entries1[$key][0][1] + $heavy_entries1[$key][0][2];
							if($next_intens > 0) {
								$next_peprelints = array(number_format($light_entries1[$key][0][0]/$next_intens, 3, ".", ""), number_format($light_entries1[$key][0][1]/$next_intens, 3, ".", ""), number_format($light_entries1[$key][0][2]/$next_intens, 3, ".", ""), 
									number_format($heavy_entries1[$key][0][1]/$next_intens, 3, ".", ""), number_format($heavy_entries1[$key][0][2]/$next_intens, 3, ".", ""));
						
								$next_deconv = $this->getOptimalRatio($ints_entries1[$key], $next_peprelints);
							}
						}
						else {

							$next_deconv = $this->getPeptideHeavy2LightRatioWithTheorInts($light_entries1[$key], $heavy_entries1[$key], $ints_entries1[$key]);
						}
						fwrite($fout, $next_deconv . "\t");
						if($next_deconv != "N/A" && $next_deconv !== "") {
							$next_intens = $light_entries1[$key][0][0] + $light_entries1[$key][0][1] + $light_entries1[$key][0][2] + $heavy_entries1[$key][0][1] + $heavy_entries1[$key][0][2];
							if($next_intens > 0) {
								array_push($peptide_ints, $next_intens);
								array_push($peptide_ratios, $next_deconv);
								array_push($peptide_ids, "A:{$key}");
								array_push($peptide_theor_ints, join(",", $ints_entries1[$key]));

								$next_peprelints = array(number_format($light_entries1[$key][0][0]/$next_intens, 3, ".", ""), number_format($light_entries1[$key][0][1]/$next_intens, 3, ".", ""), number_format($light_entries1[$key][0][2]/$next_intens, 3, ".", ""), 
									number_format($heavy_entries1[$key][0][1]/$next_intens, 3, ".", ""), number_format($heavy_entries1[$key][0][2]/$next_intens, 3, ".", ""));
								array_push($peptide_relints, join(",", $next_peprelints));

								array_push($peptide_mzs, $mzs1[$key]);
							}
						}
						
						
						fwrite($fout, join(",", $mzs1[$key]) . "\t");
						fwrite($fout, join(",", $ints_entries1[$key])."\t");
						fwrite($fout, (count($ints_entries1[$key]) > 0 && $ints_entries1[$key][0] > 0 ? number_format($ints_entries1[$key][1]/$ints_entries1[$key][0], 2) : "")."\t"); // 0 not found
						fwrite($fout, ($value[0][0] > 0 ? number_format($value[0][1]/$value[0][0], 2) : "")."\t");
						fwrite($fout, number_format($value[0][0], 0) . "\t" . number_format($value[0][1], 0) . "\t" . number_format($value[0][2], 0) . "\t");
						fwrite($fout, number_format($heavy_entries1[$key][0][0], 0) . "\t" . number_format($heavy_entries1[$key][0][1], 0) . "\t" . number_format($heavy_entries1[$key][0][2], 0) . "\t");
						if($deconvolute) {
							fwrite($fout, join("\t", $deconv_ratio_entries1[$key]) . "\t");
						}
						else {
							fwrite($fout, join("\t", $ratios_entries1[$key]) . "\t");
						}
														
					}
					else {
						fwrite($fout, "\t\t\t\t\t\t\t\t\t\t\t\t\t\t");
						if($find_residuals) {
							fwrite($fout, "\t\t");
						}
					}
				}
				fwrite($fout, "--Reporter---");
				if($deconvolute) fwrite($fout, join(",", $ints3)."\t"); 
				fwrite($fout, number_format($lightreporter[0][0], 0) . "\t" . number_format($lightreporter[0][1], 0) . "\t" . number_format($lightreporter[0][2], 0) . "\t");
				fwrite($fout, number_format($heavyreporter[0][0], 0) . "\t" . number_format($heavyreporter[0][1], 0) . "\t" . number_format($heavyreporter[0][2], 0) . "\t");
				fwrite($fout, join("\t", $ratios3));
				for($z = 0; $z < 3; $z++) {
					if($ratios3[$z] != "N/A") {
						array_push($reporter_ints, $lightreporter[0][$z] + $heavyreporter[0][$z]);
						array_push($reporter_ratios, $ratios3[$z]);
						array_push($reporter_ids, "R:".($z+1));
					}
				
				}
	
				if($deconvolute) {
					fwrite($fout, "\t" . number_format($light1ints[0][0]/$ints1[0], 0, ",", "") . "\t" . number_format($light1ints[0][1]/$ints1[1], 0, ",", "") . "\t" 
						. number_format($light2ints[0][0]/$ints2[0], 0, ",", "") . "\t" . number_format($light2ints[0][1]/$ints2[1], 0, ",", ""));
				}

				foreach($next_hk_info as $key => $value) {
					// light mass is key
					for($charge = 1; $charge <= 3; $charge++) {
						if(($key+$pep_massdiff)/$charge >= 500 && ($key+$pep_massdiff)/$charge <= 2000) {
							$light_next = $this->getIsotopeQuant($key, $charge, $scan, $mzxml, $noise, "light{$key}_1 ");
							$heavy_next = $this->getIsotopeQuant($key + $pep_massdiff , $charge, $scan, $mzxml, $noise, "heavy{$key}_1 ");
							if(count($light_next) > 0 && count($heavy_next) > 0 && array_sum($light_next[0]) > 0 && array_sum($heavy_next[0]) > 0) {
							$ratios_next = $this->printRatios($light_next[0], $heavy_next[0], "Scan {$scan} Pep1 {$key}_Ratios: ");

							if($deconvolute) {
								$modpeps = $pepxml->getModifiedPeptidePair();
								$deconv_ratio_next = $this->computeRatioFromIsotopeSetsOfPeptide($light_next[0],  $key);

							}
							$ints1 = array();
							if($deconvolute) {
								$ratios_next = $this->computeRatioFromIsotopeSets($light_next[0], $heavy_next[0], 2, $key);
								$ints_next = array_pop($ratios_next);
							
							}
							else {
								$modpeps = $pepxml->getModifiedPeptidePair();
								$ints_next = $this->getIsotopeIntensities($key, true);
							}
							if($deconvolute) $dlight_next = $this->getIsotopeQuant($key, $charge, $scan, $mzxml, $noise, "light{$key}_1 ", $num_isotopes);
							$mzs2_next = array(($key + $next_charge * $proton_mass)/$charge, ($key + $pep_massdiff + $charge * $proton_mass)/$charge);
							$ints_next = $this->getIsotopeIntensities($key, true);

							$light2ints = $this->getIsotopeQuant($key, 1, $scan, $mzxml, $noise, "light1 ");
							$heavy2ints = $this->getIsotopeQuant($key + $pep_massdiff , 1, $scan, $mzxml, $noise, "heavy1 ");
							if($deconvolute) $dlight2 = $this->getIsotopeQuant($key, 1, $scan, $mzxml, $noise, "light1 ", $num_isotopes);
							fwrite($fout, "\t");
							fwrite($fout, $this->getPeptideHeavy2LightRatio($light_next, $heavy_next) . "\t");
							fwrite($fout, join(",", $mzs2_next) . "\t");
							fwrite($fout, join(",", $ints_next)."\t");
							fwrite($fout, ($ints_next[0] > 0 ? number_format($ints_next[1]/$ints_next[0], 2) : "")."\t");
							fwrite($fout, ($light_next[0][0] > 0 ? number_format($light_next[0][1]/$light_next[0][0], 2) : "")."\t");
							fwrite($fout, number_format($light_next[0][0], 0) . "\t" . number_format($light_next[0][1], 0) . "\t" . number_format($light_next[0][2], 0) . "\t");
							fwrite($fout, number_format($heavy_next[0][0], 0) . "\t" . number_format($heavy_next[0][1], 0) . "\t" . number_format($heavy_next[0][2], 0) . "\t");
							if($deconvolute) {
								fwrite($fout, join("\t", $deconv_ratio_next) . "\t");
							}
							else {
								fwrite($fout, join("\t", $ratios_next) . "\t");
								$next_run = $mzxml;
								$nextpos = strrpos($next_run, "/");
								if($nextpos!==false) $next_run = substr($next_run, $nextpos+1);
								$next_run = substr($next_run, 0, strlen($next_run) - 6);

								echo $next_run . "\t" . $scan . "\t" . $pepxml->getPreservedOrderedModifiedPeptidePair()."\t" . number_format($masses[0], 2) . "\t" . 
								number_format($masses[1], 2) . "\t" . 
								$key . "\t" . $charge . "\t" . $value[0] . "\t" . $this->getPeptideHeavy2LightRatio($light_next, $heavy_next) .
									"\t" . join("\t", $light_next[0])."\t".$heavy_next[0][1] . "\t" . $heavy_next[0][2] . "\t" . join(",", $ints_next) . "\n";
							}
						} // count greater for light and heavy
					}


				} // next charge
			}
			fwrite($fout, "\t{$next_actual_ratio}");
			$reporter_ratios2 = array();
			$reporter_ints2 = array();
			$reporter_ids2 = array();
			for($k = 0; $k < count($reporter_ratios); $k++) {
				if($reporter_ratios[$k]=="" || $reporter_ratios[$k] <= 0 || $reporter_ratios[$k]=="inf") continue;   // 011020
				array_push($reporter_ratios2, $reporter_ratios[$k]);
				array_push($reporter_ints2, $reporter_ints[$k]);
				array_push($reporter_ids2, $reporter_ids[$k]);
			}
			$longarm_reporter_ratios2 = array();
			$longarm_reporter_ints2 = array();
			$longarm_reporter_ids2 = array();
			$longarm_scans2 = array();
			for($k = 0; $k < count($longarm_reporter_ratios); $k++) {
				if($longarm_reporter_ratios[$k]=="" || $longarm_reporter_ratios[$k] <= 0 || $longarm_reporter_ratios[$k]=="inf") continue;   // 011020
				array_push($longarm_reporter_ratios2, $longarm_reporter_ratios[$k]);
				array_push($longarm_reporter_ints2, $longarm_reporter_ints[$k]);
				array_push($longarm_reporter_ids2, $longarm_reporter_ids[$k]);
				array_push($longarm_scans2, $seen_longarm_scans[$k]);
			}
			$peptide_ratios2 = array();
			$peptide_ints2 = array();
			$peptide_ids2 = array();
			$peptide_mzs2 = array();
			$peptide_theorints2 = array();
			$peptide_relints2 = array();
			for($k = 0; $k < count($peptide_ratios); $k++) {
				if($peptide_ratios[$k]=="" || $peptide_ratios[$k] <= 0 || $peptide_ratios[$k]=="inf") continue;   // 011020
				array_push($peptide_ratios2, $peptide_ratios[$k]);
				array_push($peptide_ints2, $peptide_ints[$k]);
				array_push($peptide_ids2, $peptide_ids[$k]);
				array_push($peptide_mzs2, join("-", $peptide_mzs[$k]));
				array_push($peptide_theorints2, $peptide_theor_ints[$k]);
				array_push($peptide_relints2, $peptide_relints[$k]);
			}

			fwrite($fout, "\t" . join(";", $reporter_ids2));
			fwrite($fout, "\t" . join(";", array_map("log2ratio", $reporter_ratios2)) . "\t" . join(";", $reporter_ints2));

			fwrite($fout, "\t" . join(";", $peptide_ids2));
			fwrite($fout, "\t" . join(";", array_map("log2ratio", $peptide_ratios2)) . "\t" . join(";", $peptide_ints2));
			fwrite($fout, "\t" . join(";", $peptide_mzs2). "\t" . join(";", $peptide_theorints2) . "\t" . join(";", $peptide_relints2));


			fwrite($fout, "\t" . join(";", $frag_fragments));
			fwrite($fout, "\t" . join(";", array_map("log2ratio", $frag_ratios)) . "\t" . join(";", $frag_ints));
			fwrite($fout, "\t" . join(";", $frag_mzs). "\t" . join(";", $fragment_theor_ints) . "\t" . join(";", $fragment_relints));
			fwrite($fout, "\t" . join(";", $frag_scans));
			fwrite($fout, "\t" . join(";", $frag_noise));
			fwrite($fout, "\t" . join(";", $longarm_reporter_ids2));
			fwrite($fout, "\t" . join(";", array_map("log2ratio", $longarm_reporter_ratios2)) . "\t" . join(";", $longarm_reporter_ints2));
			fwrite($fout, "\t" . join(";", $longarm_scans2));

if(! $this->use_optimal_ratios) {			
			fwrite($fout, "\t");
			$optimal_pep_ratios = array();
			for($z = 0; $z < count($peptide_relints2); $z++) {
				$next_th = explode(",", $peptide_theorints2[$z]);
				$next_ob = explode(",", $peptide_relints2[$z]);
				array_push($optimal_pep_ratios, $this->getOptimalRatio($next_th, $next_ob));
			}
			fwrite($fout, join(";", array_map("log2ratio", $optimal_pep_ratios)));
			fwrite($fout, "\t");
			if(count($nextFragInfo) > 0) {
				$optimal_frag_ratios = array();
				for($z = 0; $z < count($fragment_relints); $z++) {
					$next_th = explode(",", $fragment_theor_ints[$z]);
					$next_ob = explode(",", $fragment_relints[$z]);
					array_push($optimal_frag_ratios, $this->getOptimalRatio($next_th, $next_ob));
				}
				fwrite($fout, join(";", array_map("log2ratio", $optimal_frag_ratios)));
			}
}	// only if not use otpimal		
			$next_pepA = preg_replace("/L/", "I", $pepxml->stripMods($pep));
			if(count($protein_seqs_and_genes) > 0 && ! array_key_exists($peptide_pair[1], $protein_seqs_and_genes)) {
				echo "Error: protein {$peptide_pair[1]} not found in protein seqs: \n"; //.join(",", array_keys($protein_seqs_and_genes))."\n";
				exit(1);
			}
			$nextposA = strpos($protein_seqs_and_genes[$peptide_pair[1]]['sequence'], $next_pepA);
			if($nextposA===false) {
				echo "ErrorA: {$next_pepA} not found in {$peptide_pair[1]} sequence {$protein_seqs_and_genes[$peptide_pair[1]]['sequence']}\n";
				exit(1);
			}
			$nextResA = $nextposA + $peptide_pair[2] + 1; //getCrosslinkProteinRes($next_pepA, $peptide_pair[2], $protein_seqs_and_genes[$peptide_pair[1]]);
			fwrite($fout, "\t" . $nextResA);
			fwrite($fout, "\t" . $protein_seqs_and_genes[$peptide_pair[1]]['gene']);			
			fwrite($fout, "\t" . $noise . "\n");


			
			}
			$peps = array();
			$masses = array(); 
			$charges = array();
			$next_scans = array();
			$next_noise = array();
			$ms3scans = array();
		}

	}
	fclose($fout);
	echo "Results written to {$outfile}\n";
	echo "\n";

	$fout_lh = fopen($this->scan_ratiofile_lh, "w");
	fwrite($fout_lh, "spectrum\tpepA\tproA\tkposA\tpepB\tproB\tkposB\tlight_id\theavy_id\n");
	foreach($seen as $xl => $light_heavy) {
		$next = explode(":", $xl); // spectrum modpep
		$nextpos = strpos($next[1], "K[");
		if($nextpos!==false) {
			fwrite($fout_lh, $next[0] . "\t" . $next[1] . "\t--\t{$nextpos}\t{$next[1]}\t--\t{$nextpos}\t" . join("\t", $light_heavy) . "\n");
		}
	}
	fclose($fout_lh);
	echo "Light-heavy identification info written to {$this->scan_ratiofile_lh}\n";
	echo "\n";



}

	public function generateDeadendUploadFile($sampledef) {
		if(! $this->deadend && ! $this->intra_only) {
			echo "Error: cannot generate upload file unless DE or INTRA is included in the option to specify deadends or intra_only\n";
			exit(1);
		}
		$min_num_reps = 3;
		$max_conf = 0.3;
		
		$myfile = fopen($sampledef, "r");
	    $first = true;
	    $exp_ref_conds = array(); // hashed sample to 
		while($line=fgets($myfile)){
			$line = rtrim($line);
			if($first) {
				$first = false;
			}
			else {
				$next = explode("\t", $line);
				$exp_ref_conds[$next[0]] = $next[1] . "_,_" . $next[2] . "_,_iqpir_,_";
			}
		}
		fclose($myfile);
		ksort($exp_ref_conds);
		$protratios = array(); // hashed by protein to ratio, std, numreps, pvalue
		echo "Read in information for ".count($exp_ref_conds). " samples in sample def file {$sampledef}\n"; 
		$outfile = $this->iqpir_params['iqpir_output_dir'] . "proteinratios_".join("_", array_keys($exp_ref_conds)).".txt";
		$prots = array();
		foreach($exp_ref_conds as $sample => $cond) {
			// get the protpair file
			$protpairfile = $this->iqpir_params['iqpir_output_dir'] . "sample." . $sample . ".protpairs.txt";
			if(! file_exists($protpairfile)) {
				echo "Error: {$protpairfile} does not exist for sample {$sample} in {$protpairfile}\n";
				exit(1);
			}
			$protratios[$sample] = array(); // hashed protein to logratio, std, num , pval
			$myfile = fopen($protpairfile, "r");
	    	$first = true;
			while($line=fgets($myfile)){
				$line = rtrim($line);
				if($first) {
					$first = false;
				}
				else {
					$next = explode("\t", $line);
					if($min_num_reps > 0 && $next[4] < $min_num_reps) continue;
					if($max_conf > 0 && $next[5] > $max_conf) continue;
					
					$next_prot = explode("-", $next[1]);
					$protratios[$sample][$next_prot[0]] = array($next[2], $next[3], $next[4], $next[8]);
					$prots[$next_prot[0]] = 1;
				}
			}
			fclose($myfile);
		}
		$fout = fopen($outfile, "w");
		fwrite($fout, "Protein");
		foreach($exp_ref_conds as $sample => $cond) fwrite($fout, "\t{$cond}log2ratio\t{$cond}log2stdev\t{$cond}numreps\t{$cond}pvalue");
		fwrite($fout, "\n");
		ksort($prots);
		foreach($prots as $pro => $value) {
			fwrite($fout, $pro);
			foreach($exp_ref_conds as $sample => $cond) {
				if(array_key_exists($pro, $protratios[$sample])) fwrite($fout, "\t".join("\t", $protratios[$sample][$pro]));
				else fwrite($fout, "\t\t\t\t");
			}
			fwrite($fout, "\n");
		}
		fclose($fout);
		echo "\nProtein ratio information for upload written to {$outfile}\nUse this file to upload proteome quantitation to a dataset on XLinkDB\n\n";
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
				$next = explode("\t", $line);
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
				$next = explode("\t", $line);
				if($next[4] < 3) continue; // num reps
				$nextprots = explode("-", $next[0]);
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
				$next = explode("\t", $line);
				if($next[6] < 1) continue; // num reps
				$nextprots = explode("-", $next[2]);
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
			echo "Composite information written to {$outfile}\n"; 
		} // next sample
	
	}
	
	public function combinedDeadendIntralinkProteinRatios($samples = array()) {
		foreach($this->iqpir_params['samples'] as $sample => $value) {
			if(count($samples) > 0 && ! array_key_exists($sample, $samples)) continue;
			$de_ratiofile = $this->iqpir_params['iqpir_output_dir'] . $this->deadend_subdir . "sample." . $sample . ".ratios.txt";
			$intra_ratiofile = $this->iqpir_params['iqpir_output_dir'] . $this->intra_only_subdir . "sample." . $sample . ".ratios.txt";
			$de_protpairfile = $this->iqpir_params['deadendprophetfile'];
			$output_dir = $this->iqpir_params['iqpir_output_dir'] . $this->deadend_subdir;
			$file_suffix = "-iqpir.xls"; //$deconvolute ? ".quantd_test{$output_index}{$filter_suff}.xls" : ".quant_test{$output_index}{$filter_suff}.xls";

			$nextpos = strrpos($de_protpairfile, "/");
			if($nextpos!==false) {
				$de_protpairfile = $output_dir . substr($de_protpairfile, $nextpos+1); // . $file_suffix;
			}
			$nextpos = strpos($de_protpairfile, ".pep.xml");
			if($nextpos!==false) $de_protpairfile = substr($de_protpairfile, 0, $nextpos);
			$de_protpairfile .= $file_suffix;
			
			$intra_protpairfile = $this->scan_ratiofile; //$this->iqpir_params['iqpir_output_dir'] . $this->intra_only_subdir . "sample." . $sample . ".protpairs.txt";
			if(! file_exists($de_ratiofile)) {
				echo "Error: ratiofile {$de_ratiofile} not found\n";
				exit(1);
			}
			if(! file_exists($intra_ratiofile)) {
				echo "Error: ratiofile {$intra_ratiofile} not found\n";
				exit(1);
			}
			$proteinratios = array(); // biorep to array
			$protpair_unis = array();
			$myfile = fopen($de_ratiofile, "r");
			$first = true;
			while($line = rtrim(fgets($myfile))) {
				if($first) {
					$first = false;
					continue;
				}
				$next = explode("\t", $line);
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
				$next = explode("\t", $line);
				if($next[6] == 0) continue; // num reps
				$nextprots = explode("-", $next[2]);
				if(! array_key_exists($nextprots[0], $proteinratios)) $proteinratios[$nextprots[0]] = array();
				array_push($proteinratios[$nextprots[0]], $next[5]);
			}
			fclose($myfile);
			$myfile = fopen($de_protpairfile, "r");
			$first = true;
			while($line = rtrim(fgets($myfile))) {
				if($first) {
					$first = false;
					continue;
				}
				$next = explode("\t", $line);
				$protpair_unis[$next[count($next)-2]] = $next[4];
			}
			fclose($myfile);
			$myfile = fopen($intra_protpairfile, "r");
			$first = true;
			$is_looplink = false;
			while($line = rtrim(fgets($myfile))) {
				if($first) {
					$first = false;
					continue;
				}
				$next = explode("\t", $line);
				
				$nextprots = explode("-", $next[count($next)-2]);

				$protpair_unis[$nextprots[0]] = $next[4];
				$protpair_unis[$nextprots[1]] = $next[7];
			}
			fclose($myfile);
			echo "Read in ".count($protpair_unis)." unis for ".count($proteinratios)." protein ratios for sample {$sample}\n"; 
			
			$add_zscore2ratiofile = false;
			$prot_written = false;
			$remove_protpair_outliers = true;
			$first = true;
			foreach($proteinratios as $id => $ratios) {
				if(count($ratios) < 2) {
					if($add_zscore2ratiofile) unset($protPairMeans[$id]);
					continue;
				}
				$protpairfile = $this->iqpir_params['iqpir_output_dir'] . $this->deadend_intra_subdir . "sample." . $sample . ".protpairs.txt";
				if(! file_exists($this->iqpir_params['iqpir_output_dir'] . $this->deadend_intra_subdir)) system("mkdir ".$this->iqpir_params['iqpir_output_dir'] . $this->deadend_intra_subdir);
				if($first) {
					$fprot = fopen($protpairfile, "w");
					fwrite($fprot, "protein\tuniprot\tmean_logratio\tstdev_logratio\tnumreps\t95conf\ttstat\tdf\tpvalue");
					if($remove_protpair_outliers) fwrite($fprot, "\toutliers");
					fwrite($fprot, "\n");
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
		} // next sample
	
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
				$next = explode("\t", $line);
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
				$next = explode("\t", $line);
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
		$tally = array();
		foreach($this->iqpir_params['samples'] as $sample => $value) {
			$de_ratiofile = $this->iqpir_params['iqpir_output_dir'] . $this->deadend_subdir . "sample." . $sample . ".ratios.txt";
			$xl_ratiofile = $this->iqpir_params['iqpir_output_dir'] . "sample." . $sample . ".ratios.txt";
			if(! file_exists($de_ratiofile)) {
				echo "Error: ratiofile {$de_ratiofile} not found\n";
				exit(1);
			}
			if(! file_exists($xl_ratiofile)) {
				echo "Error: ratiofile {$xl_ratiofile} not found\n";
				exit(1);
			}
			$deadendsites = array(); // biorep to array
			$myfile = fopen($de_ratiofile, "r");
			$first = true;
			while($line = rtrim(fgets($myfile))) {
				if($first) {
					$first = false;
					continue;
				}
				$next = explode("\t", $line);
				if($next[6] == 0) continue; // num reps
				$deadendsites[$next[2] . "-" . $next[3]] = 1; // protein-res
			}
			fclose($myfile);
			$myfile = fopen($xl_ratiofile, "r");
			$first = true;
			$tot = 0;
			$tally[$sample] = array("tot_des" => count($deadendsites), "tot_xls" => 0, "xl_0" => 0, "xl_1" => 0, "xl_2" => 0, "xl_2hd" => 0); 
			while($line = rtrim(fgets($myfile))) {
				if($first) {
					$first = false;
					continue;
				}
				$next = explode("\t", $line);
				if($next[6] == 0) continue; // num reps
				$prots = explode("-", $next[2]);
				$sites = explode("-", $next[3]);
				$num = 0;
				if(array_key_exists($prots[0] . "-" . $sites[0], $deadendsites)) $num++;
				if(array_key_exists($prots[1] . "-" . $sites[1], $deadendsites)) $num++;
				$is_homodimer = $num == 2 && $prots[0]==$prots[1] && $sites[0] == $sites[1] ? "hd" : "";
				$tally[$sample]["xl_" . $num . $is_homodimer]++;
				$tally[$sample]["tot_xls"]++;
			}
			fclose($myfile);
			if($tally[$sample]["tot_xls"] > 0) {
				$tally[$sample]["xl_0"] = number_format($tally[$sample]["xl_0"] / $tally[$sample]["tot_xls"], 3, ".", "");
				$tally[$sample]["xl_1"] = number_format($tally[$sample]["xl_1"] / $tally[$sample]["tot_xls"], 3, ".", "");
				$tally[$sample]["xl_2"] = number_format($tally[$sample]["xl_2"] / $tally[$sample]["tot_xls"], 3, ".", "");
				$tally[$sample]["xl_2hd"] = number_format($tally[$sample]["xl_2hd"] / $tally[$sample]["tot_xls"], 3, ".", "");
			}
		} // next sample
		$first = true;
		ksort($tally);
		foreach($tally as $sample => $stats) {
			if($first) {
				echo "sample\t".join("\t", array_keys($stats))."\n";
				$first = false;
			}
			echo $sample . "\t" . join("\t", array_values($stats))."\n";
		}
	}

}


?>
