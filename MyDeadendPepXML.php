<?php


class MyDeadendPepXML
{

	private $pepxml = NULL;
	private $line = NULL; // next line read from pepxml
	private $fhandle = NULL;
	private $analysis_type = ""; // ReACT or Mango
	private $interprophet = false;
	private $header = true;
	private $index = NULL;
	private $modmasses = array();
	private $analysis_types = NULL;
	private $base_name = "";
	private $probability = -1;
	private $composite_prob = -1;
	private $expect = 9999;
	private $fdr = NULL; //0.01; // what value to filter by
	private $filter_criteria = "probability";  //"composite_probability"; // or probability or expect
	private $thresh = 0; // the threshold found for filtering data by criterion
	private $composite_on = false;
	private $minprob = 0.5;
	private $maxexpect = 0.2;
	private $decoy_expect_fdr = 100;
	private $deadend_peptide = ""; // at start of spectrum query, reset all of these if have members
	private $deadend_protein = "";
	private $deadend_modpep = ""; 
	private $deadend_pos = ""; 
	private $deadend_mass = ""; 
	private $current_protein = "";
	private $current_proteins = array();
	private $current_annotation = "";
	private $annotations = array();
	private $proteins = array();
	private $modpep = false;
	private $deadend_highres_modpep = array(); // hashed from pos to mass
	private $use_highres_modpeps = true; // adds two decimal places to all variable modifications
	private $legal_criteria = array("probability" => 1, "expect" => 1);
	private $current_peptide_info = NULL;
	private $spectrum_required_text = NULL; //"041918_HeLa_H20nMTXL_LD_F8_4hr_2"; //NULL; //"032218_HeLa_LH_F10_4hr1"; //NULL; // require this in all spectra in order to pass filtering criteria
	private $spectrum = NULL;
	private $common_modmasses = array("K:325.127385", "K:333.141584", "K:329.23223", "K:329.154689", "K:329.152492", "K:327.134095");  // these indicate the crosslink site in peptide, as written in the pepxml mod_aminoacid_mass mass tag
	private $included_runs = NULL;
	private $ordered_pairs = array(); // preferentially 
	private $remove_silac_and_cys_mods = false; // set this to false in future
	private $reporter_mass = 751.40508;
	private $email_address = NULL;
	private $displaySpectrum = false;
	private $modmass_verified = array();
	private $max_verified_moddiff = 0.0001; // if the specified mod mass is more than this different from observed, won't subst observed
	private $remove_silac_lys_arg_mods = true;
	private $current_annotations = array();
	private $analyze = false;
	private $is_decoy = false;
	private $decoy_pref = "rev_";
	private $libra_ints = NULL;
	
	public function __construct ($pepxml, $fdr, $params = array()) {
		echo "Opening pepxml file {$pepxml}";
		if($fdr >= 0) {
	 		echo " with fdr {$fdr} \n";
	 	}
	 	else {
	 		echo " to read header information\n";
		 }
		$this->pepxml = $pepxml;
		$this->fdr = $fdr;
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
			else if($key=="filter_crit" && $value == "expect") {
				$this->filter_criteria = $value;
				echo "Setting filter criterion to {$value}\n";
			}
			else if($key=="display_spectrum" && $value != "") {
				$this->displaySpectrum = true;
				echo "Setting displaySpectrum to true\n";
			}
			else if($key=="libra_quant" && $value != "") {
				$this->libra_ints = array();
				echo "Setting libra quantitation to true\n";
			}
			else if($key=="modmass_info" && $value != "") {
				$modmassinfo = split(",", $value);
				for($z = 0; $z < count($modmassinfo); $z++) {
					if(! array_key_exists("libra_quant", $params)) echo "Here with {$modmassinfo[$z]}....\n";

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
					if(! array_key_exists("libra_quant", $params))  echo "Put next modmass {$next_info[1]} onto queue...\n";
				}
				if(count($this->modmasses)===0) {
					echo "Error: could not parse modmass info from line {$value}\n";
					exit(1);
				}
			}
			else {
				echo "Error: don't know parameter {$key} for MyPepXML constructor\n";
				exit(1);
			}
		}
		$this->fhandle = fopen($this->pepxml, "r") or exit("The pepxml file {$this->pepxml} is not found");
		$this->analysis_types = array();
		
		if(! $this->hasModAminoAcidVariableTags()) {
			$remove_silac_and_cys_mods = true;
			echo "Setting remove silac and static cys mods to true\n";
		}
	}

	public function getReporterMass() {
		return $this->reporter_mass;
	}

	private function hasModAminoAcidVariableTags() {
		$valid = array();
		exec("grep '<mod_aminoacid_mass' {$this->pepxml} | grep -c variable", $valid);
		if(count($valid) === 1) {
			//echo "Have {$valid[0]} counts for {$this->pepxml}\n";
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
		return count($this->deadend_peptide)===2 && count($this->deadend_protein)===2 && count($this->deadend_pos)===2;
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
		foreach($this->current_proteins as $key => $value) {
			if($this->remove_silac_and_cys_mods) {
				MyPepXML::removeSilacHeavyMods($key);
			}
			if($this->remove_silac_lys_arg_mods) {
				$key = $this->removeSilacLysArgMods($key); //preg_replace( "/K\[178.12\]/", "K[170.11]", $key); 
			}
			if(! array_key_exists($key, $this->proteins)) {
				$this->proteins[$key] = array();
				$nextprots = array_keys($value);
				for($k = 0; $k < count($nextprots); $k++) {
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

	public function getPreservedOrderedModifiedPeptidePair($verbose = false, $return_array = false) {
		$next = $this->getPreservedOrderedModifiedPeptide($verbose, false);
		$output = $next;
		for($k = 0; $k < count($next); $k++) array_push($output, $next[$k]);
		if(! $return_array) return join("\t", $output);
		return $output;
	}
	public function getPreservedOrderedModifiedPeptide($verbose = false, $return_array = false) {
		$next = $this->deadend_modpep."\t".$this->deadend_protein . "\t" . $this->deadend_pos;
			$ref = 4;
		if($this->remove_silac_and_cys_mods) {
			MyPepXML::removeSilacHeavyMods($next);
		}
		if($this->remove_silac_lys_arg_mods) {
			$next = $this->removeSilacLysArgMods($next); //preg_replace( "/K\[178.12\]/", "K[170.11]", $next_pair); 
		}
		$results = split("\t", $next);
		return $results;
	}
	

	public function getSeenLine() {
		return $this->deadend_peptide."_".$this->deadend_pos;
	}

	public function hasMods($mod) {

		$next_info = split(":", $mod);
		$valid = array();
		exec("grep -m 1 'mass=\"{$next_info[1]}\"' {$this->pepxml}", $valid);
		return count($valid) > 0;
	}
	
	public function getAnalysisType() { return "Comet"; }


	
	public function nextLine() {
		if($this->line = fgets($this->fhandle)) {
			if(! $this->header && ! $this->analyze && ! $this->hasStartTag("spectrum_query")) return true;
			if($this->hasStartTag("search_summary")) {
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
				$this->analyze = true;
				$this->header = false;
				$this->index = $this->getTagValue("index");
				// cancel all other arrays.....
				if($this->deadend_peptide!="") {
					$this->deadend_peptide = "";
		 			$this->deadend_protein = "";
	 				$this->deadend_modpep = ""; 
	 				$this->deadend_pos = ""; 
	 				$this->deadend_mass = ""; 
	 				$this->deadend_modpep = "";
	 				$this->current_protein = "";
	 				$this->current_proteins = array();
	 				$this->current_annotation = "";
					$this->modpep = false;
					$this->current_peptide_info = array();
					if(! is_null($this->libra_ints)) $this->libra_ints = array();
				}
				if($this->displaySpectrum) {
					$this->spectrum = $this->getTagValue("spectrum") . ":" . $this->getTagValue("precursor_neutral_mass") . ":" . $this->getTagValue("assumed_charge");
				}
				else if(! is_null($this->spectrum_required_text) || ! is_null($this->included_runs)) {
					$this->spectrum = $this->getTagValue("spectrum");
				}
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
			else if($this->hasStartTag("search_hit")) {
				$this->deadend_peptide = $this->getTagValue("peptide");
				
				if($this->use_highres_modpeps) {
					$this->deadend_highres_modpep = "";
				}
				$next_prot = $this->getTagValue("protein");
				$next_prot2 = $this->getUniprot($next_prot); //grabFromLine($next_prot, "sp|", "|");
				$this->is_decoy = $this->decoy_pref != "" && strpos($next_prot, $this->decoy_pref)!==false;
				if($next_prot2 !== "") {
					$this->deadend_protein = $next_prot2;
				}
				else {
					$this->deadend_protein = $next_prot;
				}
				$this->current_peptide_info[$next_prot2=="" ? $next_prot : $next_prot2] = array('wt' => $this->getTagValue("weight"), 'nsp' => $this->getTagValue("nsp"),
																					'peptide_prev_aa' => $this->getTagValue("peptide_prev_aa"), 
																					'peptide_next_aa' => $this->getTagValue("peptide_next_aa"),
																					'gene' => $this->getTagValue("peptide_next_aa"), 
																					'protein_descr' => $this->getTagValue("protein_descr"),
																					'protein' => $next_prot);
				$this->deadend_mass = $this->getTagValue("calc_neutral_pep_mass");
			}
			else if($this->hasEndTag("search_hit")) {
				$this->setCurrentProteinInformation(); // submit the current protein and annotation information
				$this->analyze = false;
			}
			else if($this->deadend_peptide!=="" && $this->hasStartTag("alternative_protein")) {
				$next_prot = $this->getTagValue("protein");
				$next_prot2 = $this->getUniprot($next_prot); //grabFromLine($next_prot, "sp|", "|");
				
				$this->current_peptide_info[$next_prot2=="" ? $next_prot : $next_prot2] = array('wt' => $this->getTagValue("weight"), 'nsp' => $this->getTagValue("nsp"),
																					'peptide_prev_aa' => $this->getTagValue("peptide_prev_aa"), 
																					'peptide_next_aa' => $this->getTagValue("peptide_next_aa"),
																					'gene' => $this->getTagValue("gene"), 
																					'protein_descr' => $this->getTagValue("protein_descr"),
																					'protein' => $next_prot);
			}
			else if($this->hasStartTag("mod_aminoacid_mass")) {
				if($this->isCrosslinkModMass()) {
					$this->deadend_pos = intval($this->getTagValue("position"))-1;
				}
				else {
					$this->deadend_mass -= floatval($this->getTagValue("mass"));
					$this->modpep = true;
				}
				if($this->use_highres_modpeps) { 
					if($this->remove_silac_and_cys_mods || $this->hasTag("variable")) { // removed this 070918 since some search rusults don't have this tag to distinguish variable from static mods
						$this->deadend_highres_modpep[$this->getTagValue("position")] = $this->getTagValue("mass");
					}
				}
			}
			else if($this->hasStartTag("modification_info")) {
				if($this->hasTag("mod_nterm_mass")!==false) {
					if($this->isCrosslinkModMass(true)) {
						$this->deadend_pos =  0; # nterm
					}
					else {
						$this->deadend_mass -= floatval($this->getTagValue("mod_nterm_mass"));
						$this->modpep = true;
					}
					if($this->use_highres_modpeps) {
						$this->deadend_highres_modpep['n'] = $this->getTagValue("mod_nterm_mass");
					}
				}
				if(! $this->use_highres_modpeps) {
					$this->deadend_modpep = $this->getTagValue("modified_peptide");
				}
			}
			else if($this->use_highres_modpeps && $this->hasEndTag("modification_info")) {
				$nextmod = "";
				if(array_key_exists('n', $this->deadend_highres_modpep)) {
					$nextmod .= 'n[' . $this->deadend_highres_modpep['n']."]";
				}
				for($k = 0; $k < strlen($this->deadend_peptide); $k++) {
					$nextmod .= $this->deadend_peptide[$k];
					if(array_key_exists($k+1, $this->deadend_highres_modpep)) {
						$nextmod .= '[' . sprintf("%0.2f", $this->deadend_highres_modpep[$k+1])."]";
					}
				}
				## compute it here
				$this->deadend_modpep = $nextmod;
			}
			else if($this->filter_criteria==="expect" && $this->hasStartTag("search_score") && $this->hasTagValue("name", "decoy_expect_fdr")) {
				if($this->expect <= $this->maxexpect) {
					$this->decoy_expect_fdr = $this->getTagValue("value");
				}
				else {
					$this->decoy_expect_fdr = 100;
				}
			}
			else if(! is_null($this->libra_ints) && $this->hasStartTag("intensity")) {
				array_push($this->libra_ints, $this->getTagValue("absolute"));
			}
			
			return true;
		}
		fclose($this->fhandle);
		return false;
	}
	
	public function aboveFilterThreshold() {
		if($this->is_decoy) return false;
	
		if(! is_null($this->spectrum_required_text) && strpos($this->spectrum, $this->spectrum_required_text)===false) return false;
		if(! is_null($this->included_runs) && ! array_key_exists($this->getSpectrumPrefix(), $this->included_runs)) return false;
		if(
			($this->filter_criteria==="probability" && $this->probability >= $this->thresh) ||
			($this->filter_criteria==="composite_probability" && $this->composite_prob >= $this->thresh && $this->probability >= $this->minprob) ||
			($this->filter_criteria==="expect" && $this->decoy_expect_fdr <= $this->fdr)) {
				$this->addCurrentProteins();
				$this->addCurrentAnnotation();
				return true;
		}
		return false;
	}
	
	public function getLibraQuant() {
		if(is_null($this->libra_ints)) return array();
		$tot = array_sum($this->libra_ints);
		if($tot == 0) return array();
		return array(number_format($tot, 0, ".", ""), join(",", array_map(function($item) use ($tot) 
			{ return number_format($item / $tot, 5, ".", ""); },
    					$this->libra_ints)));
	}
	
	public function getPeptide() {
		return $this->deadend_peptide;
	}
	public function getModifiedPeptidePair() {
		return array($this->deadend_modpep);
	}

	public function getModifiedPeptide() {
		return $this->deadend_modpep;
	}
	
	public function getPosition() {
		return $this->deadend_pos;
	}
	
	public function getExpect() {
		return $this->expect;
	}
	
	public function getProbability() {
		return $this->probability;
	}
	
	public function getBase_Name() {
		return $this->base_name;
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
	//echo "Here with {$tag} and {$value}.....\n";
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
			if(! $is_term && $this->deadend_peptide[$this->getTagValue("position")-1]==$this->modmasses[$k][0] &&
				$this->hasStartTag("mod_aminoacid_mass") && $this->hasTagValue("mass", $this->modmasses[$k][1])) {
				if(! $this->modmass_verified[$k]) $this->modmass_verified[$k] = true;
				return true;
			}
			if($is_term && $this->modmasses[$k][0]==="n" && $this->hasStartTag("modification_info") && $this->hasTagValue("mod_nterm_mass", $this->modmasses[$k][1])) {
				if(! $this->modmass_verified[$k]) $this->modmass_verified[$k] = true;
				return true;
			}
		}
		for($k = 0; $k < count($this->modmasses); $k++) {
			if($this->modmass_verified[$k] ) continue;
			if(! $is_term && $this->deadend_peptide[$this->getTagValue("position")-1]==$this->modmasses[$k][0] &&
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

	public function isModPep() {
		return $this->modpep;
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
			if(strpos($modpep, ".")===false) {
				$next_abbrev = number_format($this->modmasses[$k][2], 0, ".", "");
				$modpep = preg_replace( "/{$this->modmasses[$k][0]}\[{$next_abbrev}\]/", $this->modmasses[$k][0], $modpep); 	
			
			}
			else $modpep = preg_replace( "/{$this->modmasses[$k][0]}\[{$this->modmasses[$k][2]}\]/", $this->modmasses[$k][0], $modpep); 	

		}
		return $modpep;
	}
	
	public function getProtein() {
		return $this->protein;
	}
	
	public function getAnnotations() {
		return $this->annotations;
	}
	
	public function getThreshold() {
		return $this->thresh;
	}
	
	public static function stripMods($modpep) {
		return preg_replace( "/\[.*?\]/", "", $modpep);

	}
	
	private function setCurrentProteinInformation() {
		if(count($this->current_peptide_info)===0) return; // nothing to do
		
		$next_pep = $this->deadend_modpep; // must store each modpep separately
		
		$next_pep = $this->removeCrosslinkermods($next_pep);
		
		if($next_pep==="") { 
			return; // cancel
		}
		if($this->deadend_peptide!=="") {
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
				$next_gn = $value['gene']; 
				$next_descr = $value['protein_descr']; 
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
