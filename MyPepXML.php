<?php


class MyPepXML
{

	private $pepxml = NULL;
	private $line = NULL; // next line read from pepxml
	private $fhandle = NULL;
	private $analysis_type = ""; // ReACT or Mango
	private $interprophet = false;
	private $header = true;
	private $index = NULL;
	private $modmasses = NULL;
	private $analysis_types = NULL;
	private $base_name = "";
	private $probability = -1;
	private $composite_prob = -1;
	private $expect = 9999;
	private $fdr = NULL; //0.01; // what value to filter by
	private $filter_criteria = NULL;  //"composite_probability"; // or probability or expect
	private $thresh = 0; // the threshold found for filtering data by criterion
	private $composite_on = false;
	private $minprob = 0.5;
	private $maxexpect = 0.2;
	private $decoy_expect_fdr = 100;
	private $xl_peptides = array(); // at start of spectrum query, reset all of these if have members
	private $xl_proteins = array();
	private $xl_modpeps = array(); 
	private $xl_poss = array(); 
	private $xl_masses = array(); 
	private $current_proteins = array();
	private $current_annotations = array();
	private $annotations = array();
	private $proteins = array();
	private $modpep = false;
	private $xl_highres_modpeps = array(); // hashed from pos to mass
	private $use_highres_modpeps = true; // adds two decimal places to all variable modifications
	private $legal_criteria = array("probability" => 1, "composite_probability" => 1, "expect" => 1);
	private $current_peptide_info = NULL;
	private $spectrum_required_text = NULL; //"041918_HeLa_H20nMTXL_LD_F8_4hr_2"; //NULL; //"032218_HeLa_LH_F10_4hr1"; //NULL; // require this in all spectra in order to pass filtering criteria
	private $spectrum = NULL;
	private $common_modmasses = array("K:325.127385", "K:333.141584", "K:329.23223", "K:329.154689", "K:329.152492", "K:327.134095");  // these indicate the crosslink site in peptide, as written in the pepxml mod_aminoacid_mass mass tag
	private $included_runs = NULL;
	private $ordered_pairs = NULL; // preferentially
	private $ordered_peptide_pairs = NULL; 
	private $remove_silac_and_cys_mods = false; // set this to false in future
	private $reporter_mass = 751.40508;
	private $email_address = NULL;
	private $displaySpectrum = false;
	private $modmass_verified = array();
	private $max_verified_moddiff = 0.0001; // if the specified mod mass is more than this different from observed, won't subst observed
	private $remove_silac_lys_arg_mods = true;
	
	public function __construct ($pepxml, $fdr, $filter_crit, $params = array()) {
		echo "Opening pepxml file {$pepxml}";
		if($fdr >= 0) {
	 		echo " with fdr {$fdr} \n";
	 	}
	 	else {
	 		echo " to read header information\n";
		 }
		$this->pepxml = $pepxml;
		$this->fdr = $fdr;
		$this->filter_criteria = $filter_crit;
		echo "Employing filter criterion {$this->filter_criteria}\n";
		if(! array_key_exists($this->filter_criteria, $this->legal_criteria)) {
			echo "Error: {$this->filter_criteria} is illegal filtering criterion.  Please select one of ".join(",", array_keys($this->legal_criteria))."\n";
			exit(1);
		}
		foreach($params as $key => $value) {
			if($key=="use_highres_modpeps") {
				if($value != $this->use_highres_modpeps) {
					$this->use_highres_modpeps = $value;
					echo "Setting use_highres_modpeps to ".($value ? "true" : "false")."\n";
				}
			}
			else if($key=="spectrum_required_text") {
				if($value != $this->spectrum_required_text) {
					$this->spectrum_required_text = $value;
					echo "Setting spectrum_required_text to {$value}\n";
				}
			}
			else if($key=="email_address" && $value != "") {
				$this->email_address = $value;
				echo "Setting email_address to {$value}\n";
			}
			else if($key=="display_spectrum" && $value != "") {
				$this->displaySpectrum = true;
				echo "Setting displaySpectrum to true\n";
			}
			else {
				echo "Error: don't know parameter {$key} for MyPepXML constructor\n";
				exit(1);
			}
		}
		$this->fhandle = fopen($this->pepxml, "r") or exit("The pepxml file {$this->pepxml} is not found");
		$this->modmasses = array();
		$this->analysis_types = array();
		
		if(! $this->hasModAminoAcidVariableTags()) {
			$remove_silac_and_cys_mods = true;
			echo "Setting remove silac and static cys mods to true\n";
		}
	}

	private function fixSuryaProtNames($prot) {
# must get to 10 chars
$altprot = $prot;
if(strpos($prot, "MCOT")===0) {
	$altprot = "M" . substr($prot, 4);
	$nextpos = strrpos($altprot, ".");
	$altprot = substr($altprot, 0, $nextpos) . substr($altprot, $nextpos+1);
}
else if(strpos($prot, "XP_")!==false) {
	$altprot = substr($prot, strpos($prot, "XP_"));
	$altprot = "X" . substr($altprot, 3);
	$altprot = substr($altprot, 0, strlen($altprot)-3);
}
else if(strpos($prot, "WP_")===0) {
	$altprot = "W" . substr($prot, 3);
	$altprot = substr($altprot, 0, strlen($altprot)-2);
}
if(strlen($altprot)>10) {
	echo "Error with {$altprot} from {$prot}\n";
}
return $altprot;
}
	public function getReporterMass() {
		return $this->reporter_mass;
	}

	public function displayModmasses() {
		$output = array();
		for($k = 0; $k < count($this->modmasses); $k++) {
			array_push($output, join(":", $this->modmasses[$k]));
		}
		return join(",", $output);
	}
	public static function getStrippedPeptidePosOrderWithModPos($pep1, $kpos1, $pep2, $kpos2) {
		if($pep1!=$pep2) return MyPepXML::getStrippedPeptideOrder($pep1, $pep2, $kpos1, $kpos2);
		// if pep1 has lowest kpos return true
		if($kpos1 <= $kpos2) return true;
		return false;
	}
	// this one handles cases where the two peptides are the same, by sorting them according to kpos
	public function getStrippedPeptideOrderWithModPos($pep1, $pep2) {
		if($pep1!=$pep2) return MyPepXML::getStrippedPeptideOrder($pep1, $pep2);
		// if pep1 has lowest kpos return true
		if($this->xl_poss[0] <= $this->xl_poss[1]) return true;
		return false;
	}

	// returns whether current order is ok, or if false, whether must be swapped with pep2 firest
	public static function getStrippedPeptideOrder($pep1, $pep2, $kpos1 = 0, $kpos2 = 0) {
		if($pep1 == $pep2) return true;


    	$aamasses = array("A"=>71.03711, "R" => 156.10111, "N"=>114.04293, "D"=>115.02694,"C"=>103.00919,"E"=>129.04259,"Q"=>128.05858,
                                                "G"=>57.02146,"H"=>137.05891,"I"=>113.08406,"L"=>113.0840,"M"=>131.04049,"F"=>147.06841,"P"=>97.05276,
                                                "S"=>87.03203,"T"=>101.04768,"W"=>186.07931, "Y"=>163.06333, "V"=>99.06841, "K" => 128.09496, "n" => 1.00727);
                                                
   		$pep1 = MyPepXML::stripMods($pep1);
		$mass1 = 0;
		for($k = 0; $k < strlen($pep1); $k++) {
			if(! array_key_exists($pep1[$k], $aamasses)) {
				echo "Warning: Don't have mass for {$pep1[$k]} in {$pep1} of crosslink\n";
				continue;
				//exit(1);
			}
			$mass1 += $aamasses[$pep1[$k]];
		}
		$pep2 = MyPepXML::stripMods($pep2);
		$mass2 = 0;
		for($k = 0; $k < strlen($pep2); $k++) {
			if(! array_key_exists($pep2[$k], $aamasses)) {
				echo "Warning: Don't have mass for {$pep2[$k]} in {$pep2} of crosslink\n";
				continue;
			}
			$mass2 += $aamasses[$pep2[$k]];
		}
		if($mass1 < $mass2) return true;
		if($mass1 > $mass2) return false;
		if(strcmp($pep1, $pep2)===0) {
			if($kpos1 <= $kpos2) return true;
			return false;
		}
		return strcmp($pep1, $pep2) < 0;

	}

	public static function peptidePairStrippedForwardOrder($pep1, $pep2) {
		$pep1 = MyPepXML::stripMods($pep1);
		$pep2 = MyPepXML::stripMods($pep2);
		return strcmp($pep1, $pep2) <= 0; // alphabetical order of stripped peptides
	}

	private function hasModAminoAcidVariableTags() {
		$valid = array();
		exec("grep '<mod_aminoacid_mass' {$this->pepxml} | grep -c variable", $valid);
		if(count($valid) === 1) {
			if(rtrim($valid[0]) == 0) {
				return false;
			}
			return true;
		}
		else {
			echo "Error with grep '<mod_aminoacid_mass' {$this->pepxml} | grep -c variable\n";
			exit(1);
		}
	}


	private function getSpectrumPrefix() {
		$nextpos = strpos($this->spectrum, ".");
		if($nextpos===false) die("Error with format of spectrum {$this->spectrum}\n");
		return substr($this->spectrum, 0, $nextpos);
	}

	private function readSampleMapRuns($samplemap) {
		$this->included_runs = array();
		$myfile = fopen($samplemap, "r");
		while($line=fgets($myfile)){
			$line = rtrim($line);
			$next = split("\t", $line);
			$this->included_runs[$next[1]] = 1;
		}
		fclose($myfile);
		echo "Read in ".count($this->included_runs)." runs to filter pepXML from {$samplemap}\n"; //exit(1);
	}

	function grabFromLine($line, $start, $end) {
		$nextpos = strpos($line, $start);
		if($nextpos!== false) {
			$next_index = substr($line, $nextpos+strlen($start));
			// find terminus
			$nextpos = strpos($next_index, $end);
			if($nextpos!==false) {
				return substr($next_index, 0, $nextpos);
			}
		}
		return ""; // error
	}

	public function haveCrosslink() {
		return count($this->xl_peptides)===2 && count($this->xl_proteins)===2 && count($this->xl_poss)===2;
	}

	public function isIntralink() {
		if(strpos($this->xl_modpeps[0], "K[128")!==false) return true; // looplink
		return $this->xl_peptides[0]!==$this->xl_peptides[1] && $this->xl_proteins[0]=== $this->xl_proteins[1];
	}

	public function getUniprot($prot) {
		if(strpos($prot, "sp|")!==false) {
			return 	$this->grabFromLine($prot, "sp|", "|");
		}
		if(strpos($prot, "tr|")!==false) {
			return 	$this->grabFromLine($prot, "tr|", "|");
		}
		return '';
	}
	
	public function removeSilacLysArgMods($entry) {
		$entry = preg_replace( "/K\[178.12\]/", "K[170.11]", $entry); # acetyl Lys
		$entry = preg_replace( "/K\[178.16\]/", "K[170.14]", $entry); # trimethyl Lys
		$entry = preg_replace( "/K\[164.14\]/", "K[156.13]", $entry); # dimethyl Lys
		$entry = preg_replace( "/K\[150.12\]/", "K[142.11]", $entry); # methyl Lys
		$entry = preg_replace( "/R\[176.14\]/", "R[170.12]", $entry); # methyl Arg
		$entry = preg_replace( "/R\[190.15\]/", "R[184.13]", $entry); # dimethyl Arg
		return $entry;
	}
	
	public function addCurrentProteins() {
	
		if(count($this->xl_modpeps)!=2) {
			echo "Error: don't have two modpeps, just ".join(",", $this->xl_modpeps)."\n";
			exit(1);
		}
    		$next_peps = array($this->removeCrosslinkermods($this->xl_modpeps[0]), $this->removeCrosslinkermods($this->xl_modpeps[1]));
    		
    		
    		$next_id = join("_", $next_peps);
		foreach($this->current_proteins as $key => $value) {
			if($this->remove_silac_and_cys_mods) {
				MyPepXML::removeSilacHeavyMods($key);
			}
			if($this->remove_silac_lys_arg_mods) {
				$key = $this->removeSilacLysArgMods($key); //preg_replace( "/K\[178.12\]/", "K[170.11]", $key); 
			}
			if(! array_key_exists($next_id, $this->proteins)) {
				$this->proteins[$next_id] = array(array(), array()); // one for each protein
			}
			for($p = 0; $p < count($next_peps); $p++) {
				if($next_peps[$p] == $key) {
					$nextprots = array_keys($value);
					for($k = 0; $k < count($nextprots); $k++) {
						$this->proteins[$next_id][$p][$nextprots[$k]] = $value[$nextprots[$k]];
					}
				}
			}
			continue;

			if(! array_key_exists($key, $this->proteins)) {
				$this->proteins[$key] = array();
				$nextprots = array_keys($value);
				for($k = 0; $k < count($nextprots); $k++) {
				if(array_key_exists($nextprots[$k], $this->proteins[$key])) echo "Comparing new wt {$value[$nextprots[$k]][2]} with prev {$this->proteins[$key][$nextprots[$k]][2]} for {$key} and {$nextprots[$k]}\n";
					$this->proteins[$key][$nextprots[$k]] = $value[$nextprots[$k]];
				}
			}
		}
	}

	public function addCurrentAnnotation() {
		foreach($this->current_annotations as $key => $value) {
			if(! array_key_exists($key, $this->annotations)) {
				$this->annotations[$key] = $value;
			}
		}
	}

	public function getSpectrumRequiredText() {
		return $this->spectrum_required_text;
	}

	public function getSpectrum() {
		return $this->spectrum;
	}

	public static function removeSilacHeavyMods(&$modpep) {
		$silac_mods2remove = array(array("K", 136.11), array("R", 162.12), array("C", 160.03)); 
		for($k = 0; $k < count($silac_mods2remove); $k++) {
			$modpep = preg_replace( "/{$silac_mods2remove[$k][0]}\[{$silac_mods2remove[$k][1]}\]/", $silac_mods2remove[$k][0], $modpep); 	
		}
	}

	public function getPreservedOrderedModifiedLooplink($verbose = false, $return_array = false) {

		$final_pepA = "";
		for($f = 0; $f < strlen($this->xl_modpeps[0]); $f++) {
			if($f < strlen($this->xl_modpeps[0]) - 8 && substr($this->xl_modpeps[0], $f, 9) == "K[128.09]") {
				$final_pepA .= substr($this->xl_modpeps[1], $f, 9);
				$f += 8;
			}
			else if($f < strlen($this->xl_modpeps[0]) - 8 && substr($this->xl_modpeps[1], $f, 9) == "K[128.09]") {
				$final_pepA .= substr($this->xl_modpeps[0], $f, 9);
				$f += 8;
			}
			else $final_pepA .= $this->xl_modpeps[0][$f];
		}						
		$output =  $final_pepA . "\t" . $this->xl_proteins[0] . "\t" . min($this->xl_poss[0], $this->xl_poss[1]) .
			"\tLOOPLINK\t" .  $this->xl_proteins[0] . "\t" . max($this->xl_poss[0], $this->xl_poss[1]);
		if($return_array) return split("\t", $output);
		return $output;

	}

	public function getPreservedOrderedModifiedPeptidePair($verbose = false, $return_array = false, $subst_looplink = false) {
		if($subst_looplink && strpos($this->xl_modpeps[0], "K[128.09]")!==false) return 
			$this->getPreservedOrderedModifiedLooplink($verbose, $return_array);
		if(is_null($this->ordered_pairs)) {
			$this->ordered_pairs = array(); // hashed by stripped pep1 and pep2, going to prots 1 and 2
			$this->ordered_peptide_pairs = array();
		}
		$next_modpeps = join("_", $this->xl_modpeps);
		if(array_key_exists($next_modpeps, $this->ordered_peptide_pairs)) {
			if($return_array) return split("\t", $this->ordered_peptide_pairs[$next_modpeps]);
			return $this->ordered_peptide_pairs[$next_modpeps];
		}
		$next = $this->xl_peptides[0]."_".$this->xl_peptides[1];
		$prots = array($this->xl_proteins[0], $this->xl_proteins[1]);
		if(array_key_exists($next, $this->ordered_pairs)) {
		
			// first check if there's a conflict....
			$current = split("\t", $this->getOrderedModifiedPeptidePair());
			if($current[0]==$this->xl_modpeps[0] && $current[3]==$this->xl_modpeps[1]) {
				if($current[1]!=$prots[0] || $current[4]!=$prots[1]) {
					echo "Warning: have converting entry with {$this->xl_modpeps[0]} {$current[1]} {$this->xl_modpeps[1]} {$current[4]} to proteins {$prots[0]} and {$prots[1]}\n"; //{$this->ordered_pairs[$next]}\n";
				}
			}
			else if($current[0]==$this->xl_modpeps[1] && $current[3]==$this->xl_modpeps[0]) {
				if($current[1]!=$prots[1] || $current[4]!=$prots[0]) {
					echo "Warning: have converting entry with {$this->xl_modpeps[0]} {$current[1]} {$this->xl_modpeps[1]} {$current[4]} to proteins {$prots[1]} and {$prots[0]}\n"; //{$this->ordered_pairs[$next]}\n";
				}
			}
		}
		else {
			$this->ordered_pairs[$next] = $prots;
		}
		
		$next_pair = "";
		$ref = "";
		if($verbose) echo "Here with {$next} and ".join(",", $this->xl_masses)." xlmasses, ".
			join(",", $this->xl_modpeps)." xlmodpeps, ".join(",", $this->ordered_pairs[$next])." ordered pairs, ".
			join(",", $this->xl_poss)." poss\n";
		if($this->xl_masses[0] < $this->xl_masses[1] || $this->xl_modpeps[0]===$this->xl_modpeps[1]) {
			$next_pair = $this->xl_modpeps[0]."\t".$this->ordered_pairs[$next][0]."\t".$this->xl_poss[0]."\t".$this->xl_modpeps[1]."\t".$this->ordered_pairs[$next][1]."\t".$this->xl_poss[1];
			$ref = 1;
		}
		else if($this->xl_masses[0] > $this->xl_masses[1]) {
			$next_pair = $this->xl_modpeps[1]."\t".$this->ordered_pairs[$next][1]."\t".$this->xl_poss[1]."\t".$this->xl_modpeps[0]."\t".$this->ordered_pairs[$next][0]."\t".$this->xl_poss[0];
			$ref = 2;
		}
		else if(strcmp($this->xl_modpeps[0], $this->xl_modpeps[1]) < 0) {
			$next_pair = $this->xl_modpeps[0]."\t".$this->ordered_pairs[$next][0]."\t".$this->xl_poss[0]."\t".$this->xl_modpeps[1]."\t".$this->ordered_pairs[$next][1]."\t".$this->xl_poss[1];
			$ref = 3;
		}
		else {
			$next_pair = $this->xl_modpeps[1]."\t".$this->ordered_pairs[$next][1]."\t".$this->xl_poss[1]."\t".$this->xl_modpeps[0]."\t".$this->ordered_pairs[$next][0]."\t".$this->xl_poss[0];
			$ref = 4;
		}
		if($this->remove_silac_and_cys_mods) {
			MyPepXML::removeSilacHeavyMods($next_pair);
		}
		if($this->remove_silac_lys_arg_mods) {
			$next_pair = $this->removeSilacLysArgMods($next_pair); //preg_replace( "/K\[178.12\]/", "K[170.11]", $next_pair); 
		}
		$this->ordered_peptide_pairs[$next_modpeps]	= $next_pair;	
		// for heavy silac acetyl mods
		if($return_array) {
			$results = split("\t", $next_pair);
			return $results;
		}
		return $next_pair;
	}
	
	// peptide1\tprotein1\tposs1\tpeptide2\tprotein2\tpsos2
	public function getOrderedModifiedPeptidePair() {
		if($this->xl_masses[0] < $this->xl_masses[1] || $this->xl_modpeps[0]===$this->xl_modpeps[1]) {
			return $this->xl_modpeps[0]."\t".$this->xl_proteins[0]."\t".$this->xl_poss[0]."\t".$this->xl_modpeps[1]."\t".$this->xl_proteins[1]."\t".$this->xl_poss[1];
		}
		if($this->xl_masses[0] > $this->xl_masses[1]) {
			return $this->xl_modpeps[1]."\t".$this->xl_proteins[1]."\t".$this->xl_poss[1]."\t".$this->xl_modpeps[0]."\t".$this->xl_proteins[0]."\t".$this->xl_poss[0];
		}
		if(strcmp($this->xl_modpeps[0], $this->xl_modpeps[1]) < 0) {
			return $this->xl_modpeps[0]."\t".$this->xl_proteins[0]."\t".$this->xl_poss[0]."\t".$this->xl_modpeps[1]."\t".$this->xl_proteins[1]."\t".$this->xl_poss[1];
		}
		return $this->xl_modpeps[1]."\t".$this->xl_proteins[1]."\t".$this->xl_poss[1]."\t".$this->xl_modpeps[0]."\t".$this->xl_proteins[0]."\t".$this->xl_poss[0];
	}
	
	public function getSeenLine() {
		return $this->xl_peptides[0]."_".$this->xl_poss[0]."_".$this->xl_peptides[1]."_".$this->xl_poss[1];
	}
	// peptide1\tprotein1\tposs1\tpeptide2\tprotein2\tpsos2
	public function getOrderedPeptidePair() {
		if($this->xl_masses[0] < $this->xl_masses[1] || $this->xl_peptides[0]===$this->xl_peptides[1]) {
			return $this->xl_peptides[0]."\t".$this->xl_proteins[0]."\t".$this->xl_poss[0]."\t".$this->xl_peptides[1]."\t".$this->xl_proteins[1]."\t".$this->xl_poss[1];
		}
		if($this->xl_masses[0] > $this->xl_masses[1]) {
			return $this->xl_peptides[1]."\t".$this->xl_proteins[1]."\t".$this->xl_poss[1]."\t".$this->xl_peptides[0]."\t".$this->xl_proteins[0]."\t".$this->xl_poss[0];
		}
		if(strcmp($this->xl_peptides[0], $this->xl_peptides[1]) < 0) {
			return $this->xl_peptides[0]."\t".$this->xl_proteins[0]."\t".$this->xl_poss[0]."\t".$this->xl_peptides[1]."\t".$this->xl_proteins[1]."\t".$this->xl_poss[1];
		}
		return $this->xl_peptides[1]."\t".$this->xl_proteins[1]."\t".$this->xl_poss[1]."\t".$this->xl_peptides[0]."\t".$this->xl_proteins[0]."\t".$this->xl_poss[0];
	}
	
	public function hasMods($mod) {
		$next_info = split(":", $mod);
		$valid = array();
		exec("grep -m 1 'mass=\"{$next_info[1]}\"' {$this->pepxml}", $valid);
		return count($valid) > 0;
	}
	
	public function nextLine() {
		if($this->line = fgets($this->fhandle)) {
			if($this->hasStartTag("search_summary")) {
				$comet_version = $this->getTagValue("search_engine_version");
				if($this->getTagValue("search_engine")=="SpectraST") {
					$this->analysis_type = "SpectraST";
				}
				else if($this->getTagValue("search_engine")=="LoopComet") {
					$this->analysis_type = "Looplink";
				}
				else if(strpos($comet_version, "Mango")===false) {
				 	$this->analysis_type = "ReACT";
				 }
				 else {
				 	$this->analysis_type = "Mango";
				}
				if(! array_key_exists($this->analysis_type, $this->analysis_types)) {
					$this->analysis_types[$this->analysis_type] = 1;
					echo "Adding {$this->analysis_type} to analysis list....\n";
				}
				$this->base_name = $this->getTagValue("base_name");
			}
			else if($this->header && ! $this->interprophet && $this->hasStartTag("analysis_summary") && $this->hasTagValue("analysis", "interprophet")) {
				$this->interprophet = true;
				echo "Setting interprophet status to true\n"; //exit(1);
			}
			else if($this->hasStartTag("spectrum_query")) {
				if($this->fdr < 0) {
					fclose($this->fhandle);
					return false; // done already
				}
				$this->header = false;
				$this->index = $this->getTagValue("index");
				// cancel all other arrays.....
				if(count($this->xl_peptides)>0) {
					$this->xl_peptides = array();
		 			$this->xl_proteins = array();
	 				$this->xl_modpeps = array(); 
	 				$this->xl_poss = array(); 
	 				$this->xl_masses = array(); 
	 				$this->xl_modpeps = array();
	 				$this->current_proteins = array();
	 				$this->current_annotations = array();
					$this->modpep = false;
					$this->current_peptide_info = array();
				}
				if($this->displaySpectrum) {
					$this->spectrum = $this->getTagValue("spectrum") . ":" . $this->getTagValue("precursor_neutral_mass") . ":" . $this->getTagValue("assumed_charge");
				}
				else if(! is_null($this->spectrum_required_text) || ! is_null($this->included_runs)) {
					$this->spectrum = $this->getTagValue("spectrum");
				}
			}
			else if($this->header && count($this->modmasses)===0 && $this->hasStartTag("peptideprophet_summary")) {
				// first check whether have reporter mass option
				$options = $this->getTagValue("options");
				$nextpos = strpos($options, "REPORTERMASS=");
				if($nextpos!==false) {
					$this->reporter_mass = floatVal(substr($options, $nextpos + strlen("REPORTERMASS=")));
					echo "Setting reporter mass to {$this->reporter_mass}\n";
				}
				$modmassinfo = split(",", $this->getTagValue("xlinker_stump_modmasses"));
				for($z = 0; $z < count($modmassinfo); $z++) {
					if(strlen($modmassinfo[$z]) < strlen($this->common_modmasses[0])) { // check if one of them
						for($i = 0; $i < count($this->common_modmasses); $i++) {
							if(strpos($this->common_modmasses[$i], $modmassinfo[$z])===0 && $this->hasMods($this->common_modmasses[$i])) {
								echo "Substituting modmass {$this->common_modmasses[$i]} for {$modmassinfo[$z]}\n";
								$modmassinfo[$z] = $this->common_modmasses[$i];
								$i = count($this->common_modmasses);
							}
							else {
								$next_mod = substr($this->common_modmasses[$i], 0, 2).sprintf("%0.2f", substr($this->common_modmasses[$i], 2));
								if($next_mod===$modmassinfo[$z]) {
									echo "Substituting modmass {$this->common_modmasses[$i]} for {$modmassinfo[$z]}\n";
									$modmassinfo[$z] = $this->common_modmasses[$i];
									$i = count($this->common_modmasses);
								}
							}
						}
					}
					else if(strlen($modmassinfo[$z]) > strlen($this->common_modmasses[0])) { // check if one of them
						for($i = 0; $i < count($this->common_modmasses); $i++) {
							if(strpos($modmassinfo[$z], $this->common_modmasses[$i])===0 && $this->hasMods($this->common_modmasses[$i])) {
								echo "Substituting modmass {$this->common_modmasses[$i]} for {$modmassinfo[$z]}\n";
								$modmassinfo[$z] = $this->common_modmasses[$i];
								$i = count($this->common_modmasses);
							}
						}
					}
					$next_info = split(":", $modmassinfo[$z]);
				
					// check how many decimal places to make sure as is written in pepXML
					$next_dec = strpos($next_info[1], ".");
					$num_decs = $next_dec === false ? 0 : strlen(substr($next_info[1], $next_dec+1));
					if($num_decs < 5) { // error, didn't find the full value
						echo "Error: have truncated specified modification mass {$modmassinfo[$z]} without a full length reference among: ".join(",", $this->common_modmasses)."\n";
						echo "Please add the full length modification mass (as in the pepXML mod_aminoacid_mass mass value) to your xlinkprophet.params file and the header xlinker_stump_modmasses of {$this->pepxml}\n";
						exit(1);
					}
								
					array_push($this->modmasses, array($next_info[0], strval($next_info[1]), number_format($next_info[1], 2)));
					array_push($this->modmass_verified, false);
					echo "Put next modmass {$next_info[1]} onto queue...\n";
				}
				if(count($this->modmasses)===0) {
					echo "Error: could not parse modmass info from line {$this->line}\n";
					exit(1);
				}
			}
			else if($this->hasStartTag("search_score") && $this->hasTagValue("name", "composite_probability")) {
				$this->composite_prob = $this->getTagValue("value");
			}
			else if(($this->interprophet && $this->hasStartTag("interprophet_result")) ||
					(! $this->interprophet && $this->hasStartTag("peptideprophet_result"))) {
				$this->probability = $this->getTagValue("probability");
				if($this->probability < 0) {
					$this->probability = 0;
				}
			}
			else if($this->hasStartTag("search_score") && $this->hasTagValue("name", "expect")) {
				$this->expect = $this->getTagValue("value");
			}
			else if($this->thresh===0 && $this->filter_criteria==="probability" && $this->hasStartTag("error_point")) {
				$next_error = $this->getTagValue("error");
				if($next_error==$this->fdr) {
					$this->thresh = $this->getTagValue("min_prob");
					if($this->thresh < $this->minprob) {
						$this->thresh = $this->minprob;
					}
					$this->thresh = number_format($this->thresh, 2);
					echo "Setting minimum probability to ".$this->thresh." to achieve desired FDR of ".$this->fdr."<br/>\n";
				}
			}
			else if($this->filter_criteria==="composite_probability" && ! $this->composite_on && $this->hasStartTag("roc_error_data") && $this->hasTagValue("charge", "composite")) {
				$this->composite_on = true;
			}
			else if($this->filter_criteria==="composite_probability" && $this->composite_on && $this->hasStartTag("error_point")) {
				$next_error = $this->getTagValue("error");
				if($next_error==$this->fdr) {
					$this->thresh = $this->getTagValue("min_prob");
					if($this->thresh < $this->minprob) {
						$this->thresh = $this->minprob;
					}
					$this->thresh = number_format($this->thresh, 2);
					echo "Setting minimum composite probability to ".$this->thresh." to achieve desired composite FDR of ".$this->fdr."<br/>\n";
					$composite_on = false;
					$use_composite_fdr = false; // done
				}
			}
			else if($this->hasStartTag("linked_peptide")) {
				array_push($this->xl_peptides, $this->getTagValue("peptide"));
				
				if($this->use_highres_modpeps) {
					$this->xl_highres_modpeps = array();
				}
				$next_prot = $this->getTagValue("protein");
				$next_prot2 = $this->getUniprot($next_prot); //grabFromLine($next_prot, "sp|", "|");
		
				if($next_prot2 !== "") {
					array_push($this->xl_proteins, $next_prot2);
				}
				else {
					array_push($this->xl_proteins, $next_prot);
				}
				$this->current_peptide_info[$next_prot2=="" ? $next_prot : $next_prot2] = array('wt' => $this->getTagValue("weight"), 'nsp' => $this->getTagValue("nsp"),
																					'peptide_prev_aa' => $this->getTagValue("peptide_prev_aa"), 
																					'peptide_next_aa' => $this->getTagValue("peptide_next_aa"),
																					'gene' => $this->getTagValue("peptide_next_aa"), 
																					'protein_descr' => $this->getTagValue("peptide_next_aa"),
																					'protein' => $next_prot);
				array_push($this->xl_masses, $this->getTagValue("calc_neutral_pep_mass"));
			}
			else if($this->hasEndTag("linked_peptide")) {
				$this->setCurrentProteinInformation(); // submit the current protein and annotation information
			}
			else if(count($this->xl_peptides)>0 && $this->hasStartTag("alternative_protein")) {
				$next_prot = $this->getTagValue("protein");
				$this->current_peptide_info[$next_prot2=="" ? $next_prot : $next_prot2] = array('wt' => $this->getTagValue("weight"), 'nsp' => $this->getTagValue("nsp"),
																					'peptide_prev_aa' => $this->getTagValue("peptide_prev_aa"), 
																					'peptide_next_aa' => $this->getTagValue("peptide_next_aa"),
																					'gene' => $this->getTagValue("gene"), 
																					'protein_descr' => $this->getTagValue("protein_descr"),
																					'protein' => $next_prot);
			}
			else if($this->hasStartTag("mod_aminoacid_mass")) {
				if($this->isCrosslinkModMass()) {
					array_push($this->xl_poss, intval($this->getTagValue("position"))-1);
				}
				else {
					if($this->hasStartTag("variable")) $this->xl_masses[count($this->xl_masses)-1] -= floatval($this->getTagValue("variable"));
					else if($this->hasStartTag("static")) $this->xl_masses[count($this->xl_masses)-1] -= floatval($this->getTagValue("static"));
					else $this->xl_masses[count($this->xl_masses)-1] -= floatval($this->getTagValue("mass")); // old format of Comet
					$this->modpep = true;
				}
				if($this->use_highres_modpeps) { 
					if($this->remove_silac_and_cys_mods || $this->hasTag("variable")) { // removed this 070918 since some search rusults don't have this tag to distinguish variable from static mods
						$this->xl_highres_modpeps[$this->getTagValue("position")] = $this->getTagValue("mass");
					}
				}
			}
			else if($this->hasStartTag("modification_info")) {
				if($this->hasTag("mod_nterm_mass")!==false) {
					if($this->isCrosslinkModMass(true)) {
						array_push($this->xl_poss, 0); # nterm
					}
					else {
						$this->xl_masses[count($this->xl_masses)-1] -= floatval($this->getTagValue("mod_nterm_mass"));
						$this->modpep = true;
					}
					if($this->use_highres_modpeps) {
						$this->xl_highres_modpeps['n'] = $this->getTagValue("mod_nterm_mass");
					}
				}
				if(! $this->use_highres_modpeps) {
					array_push($this->xl_modpeps, $this->getTagValue("modified_peptide"));
				}
			}
			else if($this->use_highres_modpeps && $this->hasEndTag("modification_info")) {
				$nextmod = "";
				if(array_key_exists('n', $this->xl_highres_modpeps)) {
					$nextmod .= 'n[' . $this->xl_highres_modpeps['n']."]";
				}
				for($k = 0; $k < strlen($this->xl_peptides[count($this->xl_peptides)-1]); $k++) {
					$nextmod .= $this->xl_peptides[count($this->xl_peptides)-1][$k];
					if(array_key_exists($k+1, $this->xl_highres_modpeps)) {
						$nextmod .= '[' . sprintf("%0.2f", $this->xl_highres_modpeps[$k+1])."]";
					}
				}
				## compute it here
				array_push($this->xl_modpeps, $nextmod);
			}
			else if($this->filter_criteria==="expect" && $this->hasStartTag("search_score") && $this->hasTagValue("name", "decoy_expect_fdr")) {
				if($this->expect <= $this->maxexpect) {
					$this->decoy_expect_fdr = $this->getTagValue("value");
				}
				else {
					$this->decoy_expect_fdr = 100;
				}
			}
			
			return true;
		}
		fclose($this->fhandle);
		return false;
	}
	
	public function aboveFilterThreshold($minprob_with_mincomposite = false) {
		if(! is_null($this->spectrum_required_text) && strpos($this->spectrum, $this->spectrum_required_text)===false) return false;
		if(! is_null($this->included_runs) && ! array_key_exists($this->getSpectrumPrefix(), $this->included_runs)) return false;
		if(
			($this->filter_criteria==="probability" && $this->probability >= $this->thresh) ||
			($this->filter_criteria==="composite_probability" && $this->composite_prob >= $this->thresh &&
				(! $minprob_with_mincomposite || $this->probability >= $this->minprob)) ||
			($this->filter_criteria==="expect" && $this->decoy_expect_fdr <= $this->fdr)) {
				$this->addCurrentProteins();
				$this->addCurrentAnnotation();
				return true;
		}
		return false;
	}
	public function getPeptidePair() {
		return $this->xl_peptides;
	}
	
	public function getModifiedPeptidePair() {
		return $this->xl_modpeps;
	}
	
	public function getProteinPair() {
		$output = array();
		for($k = 0; $k < count($this->xl_proteins); $k++) {
			if(array_key_exists($this->xl_proteins[$k], $this->current_annotations)) {
				array_push($output, $this->xl_proteins[$k]."|".$this->current_annotations[$this->xl_proteins[$k]][1]);
			}
			else {
				array_push($output, $this->xl_proteins[$k]);
			}
		}
		return $output;
		return $this->xl_proteins;
	}
	
	public function getPositionPair() {
		return $this->xl_poss;
	}
	
	public function getExpect() {
		return $this->expect;
	}
	
	public function getProbability() {
		return $this->probability;
	}
	
	public function getCompositeProbability() {
		return $this->composite_prob;
	}
	
	public function getBase_Name() {
		return $this->base_name;
	}
	
	public function getAnalysisType() {
		return $this->analysis_type;
	}
	
	public function getAnalysisTypes() {
		return $this->analysis_types;
	}
	public function isInterProphet() {
		return $this->interprophet;
	}

	public function hasStartTag($tag) {
		return ! is_null($this->line) && strpos($this->line, "<{$tag} ")!==false;
	}

	public function hasEndTag($tag) {
		return ! is_null($this->line) && strpos($this->line, "</{$tag}")!==false;
	}

	public function hasTagValue($tag, $value) {
		return ! is_null($this->line) && strpos($this->line, " {$tag}=\"{$value}\"")!==false;
	}
	public function hasTag($tag) {
		return ! is_null($this->line) && (strpos($this->line, " {$tag}=\"")!==false || strpos($this->line, "<{$tag}=\"")!==false);
	}
	public function getTagValue($tag) {
		return $this->grabFromLine($this->line, "{$tag}=\"", "\"");
	}

	public function getIndex() {
		return $this->index;
	}
	
	public function isCrosslinkModMass($is_term = false) {
		for($k = 0; $k < count($this->modmasses); $k++) {
			if(! $is_term && $this->xl_peptides[count($this->xl_peptides)-1][$this->getTagValue("position")-1]==$this->modmasses[$k][0] &&
				$this->hasStartTag("mod_aminoacid_mass") && $this->hasTagValue("mass", $this->modmasses[$k][1])) {
				if(! $this->modmass_verified[$k]) $this->modmass_verified[$k] = true;
				return true;
			}
			if($is_term && $this->modmasses[$k][0]==="n" && $this->hasStartTag("modification_info") && $this->hasTagValue("mod_nterm_mass", $this->modmasses[$k][1])) {
				if(! $this->modmass_verified[$k]) $this->modmass_verified[$k] = true;
				return true;
			}
		}
		// still here, get one chance to set the observed mass rather than the one is is listed, if very close
		
		for($k = 0; $k < count($this->modmasses); $k++) {
			if($this->modmass_verified[$k] ) continue;
			if(! $is_term && $this->xl_peptides[count($this->xl_peptides)-1][$this->getTagValue("position")-1]==$this->modmasses[$k][0] &&
				$this->hasStartTag("mod_aminoacid_mass") && abs($this->getTagValue("mass")- $this->modmasses[$k][1]) <= $this->max_verified_moddiff) {
				echo "Setting modmass to ".$this->getTagValue("mass")." based on similarity to ".$this->modmasses[$k][1]."\n"; //exit(1);
				$this->modmasses[$k][1] = $this->getTagValue("mass");
				$this->modmass_verified[$k] = true;
				return true;
			}
			if($is_term && $this->modmasses[$k][0]==="n" && $this->hasStartTag("modification_info") && 
				abs($this->getTagValue("mod_nterm_mass") - $this->modmasses[$k][1]) <= $this->max_verified_moddiff) {
				echo "Setting modmass to ".$this->getTagValue("mass")." based on similarity to ".$this->modmasses[$k][1]."\n";
				$this->modmasses[$k][1] = $this->getTagValue("mass");
				$this->modmass_verified[$k] = true;
				return true;
			}
		}
		
		return false;
	}
	
	public function getModifiedPeptideWithoutCrosslinkmodAndKpos($modpep) {
		$copy = $modpep;
		for($k = 0; $k < count($this->modmasses); $k++) {
			$modpep = preg_replace( "/{$this->modmasses[$k][0]}\[{$this->modmasses[$k][2]}\]/", $this->modmasses[$k][0] . " ", $modpep); 	
		}
		$next = split(" ", $modpep);
		if(count($next)!==2) {
			echo "Error: have ".count($next)." crosslink modifictaions in {$modpep} (orig {$copy})\n";
			exit(1);
		}
		return array(join("", $next), strlen($this->stripMods($next[0]))-1);
	}
	
	public function replaceCrosslinkermods($modpep, $replacement) {
		for($k = 0; $k < count($this->modmasses); $k++) {
			$modpep = preg_replace( "/{$this->modmasses[$k][0]}\[{$this->modmasses[$k][2]}\]/", $replacement, $modpep); 	
		}
		return $modpep;
	}
	
	public function removeCrosslinkermods($modpep) {
		for($k = 0; $k < count($this->modmasses); $k++) {
			$modpep = preg_replace( "/{$this->modmasses[$k][0]}\[{$this->modmasses[$k][2]}\]/", $this->modmasses[$k][0], $modpep); 	
		}
		return $modpep;
	}
	
	public function getProteins() {
		return $this->proteins;
	}
	
	public function getAnnotations() {
		return $this->annotations;
	}
	
	public function isModPep() {
		return $this->modpep;
	}
	
	public function getThreshold() {
		return $this->thresh;
	}
	
	public static function stripMods($modpep) {
		return preg_replace( "/\[.*?\]/", "", $modpep);

	}
	
	private function setCurrentProteinInformation() {
		if(count($this->current_peptide_info)===0) return; // nothing to do
		
		$next_pep = $this->xl_modpeps[count($this->xl_modpeps)-1]; // must store each modpep separately
		
		$next_pep = $this->removeCrosslinkermods($next_pep);
		
		if($next_pep==="" || array_key_exists($next_pep, $this->xl_proteins)) {
			return; // cancel
		}
		if(count($this->xl_peptides) == 1) {
			$this->current_proteins = array($next_pep => array()); // hashed from protein to rest.....
		}
		else {
			$this->current_proteins[$next_pep] = array();
		}
		foreach($this->current_peptide_info as $key => $value) {
			if(! array_key_exists($key, $this->current_proteins[$next_pep])) {
				$this->current_proteins[$next_pep][$key] = array($next_pep, $key, $value['wt'], $value['nsp'], $value['peptide_prev_aa'], $value['peptide_next_aa']);
			}
			if(! array_key_exists($key, $this->annotations)) {
				$next_gn = $value['gene']; //$this->getTagValue("gene");
				$next_descr = $value['protein_descr']; //$this->getTagValue("protein_descr");
				$next_prot = $value['protein'];
				$os = strpos($next_descr, "OS=");
				if($os!==false) {
					if($next_gn==="") {
						$next_gn = $this->grabFromLine($next_descr, "GN=", " ");
					}
					$next_descr = substr($next_descr, 0, $os-1);
				}
				$next_gene_pos = strrpos($next_prot, "|");
				$next_gene = "";
				if($next_gene_pos !== false) {
					$next_gene = substr($next_prot, $next_gene_pos+1);
				}
				$this->current_annotations[$key] = array($next_descr, $next_gene, $next_gn);
			}
		} // next protein to process
		$this->current_peptide_info = array(); // cancel
	}

	
}

?>
