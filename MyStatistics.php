<?php


class MyStatistics
{


	public static $index = 0;
	public static $sort_asc = true;

	public static function normalize(&$arr, $num_dec = 5) {
		$tot = array_sum($arr);
		if($tot == 0) return;
		$arr = array_map(function($item) use ($tot, $num_dec) 
				{ return number_format($item / $tot, $num_dec, ".", ""); },
    			$arr);	
	}
	
	public static function removeOutliersByIndex(&$arr, &$rejects, $index, $verbose = false) {
		if(count($arr) < 3) return false;
		if(count($arr) <= 4) {
			$quar = 0;
			$ter = count($arr)-1;
		}
		else {
			$quar = number_format(count($arr)/4, 0, ".", "");
			$ter = number_format(3*count($arr)/4, 0, ".", "");
			if($ter >= count($arr)) {
				$ter--;
			}
		}
		// set this from false if want to do second round
		$first = false; // count($rejects)===0;
		if($ter < 0) $ter = 0;
		if($quar > count($arr)-1) $quar = count($arr)-1;
		$max = max($quar, count($arr) - $ter - 1);
		if($quar < $max) $quar = $max;
		else if(count($arr) - $ter - 1 < $max) $ter = count($arr) - $max - 1;
		
		if(count($arr) > 100) {
			$quar++;
			$ter--;
		}
		// sort by index
		MyStatistics::sortMultidimensionalArrayByIndex($arr, $index);

		$iqr = $arr[$ter][$index]-$arr[$quar][$index];
		$l_index = 0;
		$r_index = count($arr)-1;
		$left_bound = $arr[$quar][$index]-1.5*$iqr;
		$right_bound = $arr[$ter][$index]+1.5*$iqr;
		while($l_index < count($arr) && $arr[$l_index][$index] < $left_bound) {
			array_push($rejects, $arr[$l_index][$index]);
			$l_index++;
		}
		while($r_index > 0 && $arr[$r_index][$index] > $right_bound) {
			if($verbose)
				echo "2. Adding to rejects: ".$arr[$r_index][$index]."\n";
			array_push($rejects, $arr[$r_index][$index]);
			$r_index--;
		}
		if($l_index > 0 || $r_index < count($arr)-1) {
			$arr = array_slice($arr, $l_index, $r_index-$l_index+1);
			if($first) {
				return MyStatistics::removeOutliersByIndex($arr, $rejects, $index, $verbose);
			}
			return true;
		}
		return false;

	}
	
	public static function sortMultidimensionalArrayByIndex(&$arr, $index, $desc = false) {
		MyStatistics::$index = $index;
		if($desc) {
			MyStatistics::$sort_asc = false;
		}
		usort($arr, array("MyStatistics", "sortByIndex"));
		MyStatistics::$index = 0; // set back to default
		if($desc) {
			MyStatistics::$sort_asc = true;
		}
	}

	public static function get_tstatistic($statarr) {
		if($statarr[2] < 2) {
			return "x"; // not available
		}
		if($statarr[1]==0) {
			return 999; // maximum value
		}
		return number_format(sqrt($statarr[2]-1)*$statarr[0]/$statarr[1], 1);
	}

	public static function get_tpvalue($tstat, $df) {
		$tstat_pvalfile = "t_test_pvals.txt";
		if(! file_exists($tstat_pvalfile)) {
			echo "Error: {$tstat_pvalfile} does not exist\n";
			exit(1);
		}
		$tstat = abs($tstat); // don't care which direction
		$valid = array();
		$dfs = array(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
			40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200);
		$dfs = array_flip($dfs);
		if(! array_key_exists($df, $dfs)) {
			$df = sprintf("%0.0f", intval($df)/10)*10;
			if(! array_key_exists($df, $dfs)) {
				$df = 200; // last one
			}
		}
		exec("grep ^{$df} {$tstat_pvalfile}", $valid);
		for($k = 0; $k < count($valid); $k++) {
			$valid[$k] = rtrim($valid[$k]);
			$pvalues = split("\t", $valid[$k]);
			if($pvalues[0] != $df) continue;
			// now need to get header so know which column for tstat
			$valid2 = array();
			exec("head -1 {$tstat_pvalfile}", $valid2);
			if(count($valid2)!==1) {
				echo "Error no head in {$tstat_pvalfile}\n";
				exit(1);
			}
			$valid2[0] = rtrim($valid2[0]);
			$table = array_map(function($item) 
				{ return substr($item, 6); },
    			split("\t", $valid2[0]));
 			// see which columns flank tstat
			if($tstat <= $table[1]) return $pvalues[1];
			for($j = 1; $j < count($table) - 1; $j++) {
				if($tstat >= $table[$j] && $tstat < $table[$j+1]) {
					return sprintf("%0.2e", ($tstat - $table[$j])*(($pvalues[$j+1]-$pvalues[$j])/($table[$j+1]-$table[$j])) + $pvalues[$j]);
				}
			}
			return $pvalues[count($pvalues)-1];
		
		}
		return "N/A";
	}
	
	public static function stats($arr, $mean_numdec = -1, $std_numdec = -1) {
		if(count($arr)===0) {
    		return array(-1, -1, 0);
		}
		if(count($arr)===1) {
			if($mean_numdec > -1 ) return array(number_format($arr[0], $mean_numdec, ".", ""), 0, 1);
			return array($arr[0], 0, 1);
		}
		$mean = 0;
		$meansq = 0;
		for($k = 0; $k < count($arr); $k++) {
        	$mean += $arr[$k];
       		$meansq += $arr[$k] * $arr[$k];
		}
		$mean /= (count($arr));
		$meansq /= (count($arr) - 1);
		$meansq -= $mean * $mean * (count($arr))/(count($arr)-1);

		if($meansq < 0) {
			$meansq = 0;
		}
		$meansq = sqrt($meansq);
		if($mean_numdec > -1) $mean = number_format($mean, $mean_numdec, ".", "");
		if($std_numdec > -1) $meansq = number_format($meansq, $std_numdec, ".", "");
		return array($mean, $meansq, count($arr));
	}
	
	public static function median($arr) {
		if(count($arr)===0) {
			echo "Error: have not data in arr\n";
			exit(1);
		}
		sort($arr);
		return $arr[number_format(count($arr)/2, 0, '.', '')];
	}
	
	public static function oneWayAnovaAlt($stats1, $stats2, $ref = "") {
		if($stats1[2] + $stats2[2] < 3) return array("", "", "");
		$n = $stats1[2] + $stats2[2];
		$mean = (($stats1[2] * $stats1[0]) + ($stats2[2] * $stats2[0]))/$n;
		$df1 = 1;
		$df2 = $n - 2;
		$ms_error = ((($stats1[2]-1) * $stats1[1] * $stats1[1]) + (($stats2[2]-1) * $stats2[1] * $stats2[1]))/$df2;
		if($ms_error == 0) return array("", "", "");
		$ms_groups = ($stats1[2] * ($stats1[0]-$mean) * ($stats1[0] - $mean)) + ($stats2[2] * ($stats2[0]-$mean) * ($stats2[0] - $mean));
		$fratio = $ms_groups/$ms_error;
		$fcrit = MyStatistics::getFcriticalVal5pct($df2);
		$pvalue = MyStatistics::getFPvalueFromRef($fratio, $df2, $ref);
		return array($fcrit, number_format($fratio, 1, ".", ""), sprintf("%0.2e",$pvalue));
	}

	public static function getFcriticalVal5pct($df2) {
		$f = array(0, 0, 0, 10.13, 7.71, 6.61, 5.59, 5.32, 5.12, 4.96, 4.84, 4.75, 4.67, 4.6,
		4.54, 4.49, 4.45, 4.41, 4.38, 4.35, 4.32, 4.3, 4.28, 4.26, 4.25, 4.23, 4.22, 4.2, 4.19, 4.17);
		if($df2 < count($f)) return $f[$df2];
		if($df2 <= 35) return 4.12;
		if($df2 <= 40) return 4.08;
		if($df2 <= 45) return 4.06;
		if($df2 <= 50) return 4.03;
		if($df2 <= 60) return 4;
		if($df2 <= 70) return 3.98;
		if($df2 <= 80) return 3.96;
		if($df2 <= 100) return 3.94;
		if($df2 <= 200) return 3.89;
		return 3.5;
	}
	
	public static function getFPvalueFromRef($fstat, $df2, $ref = "") {
		$df2 = min(30, $df2);
		$fstat = number_format($fstat, 1, ".", "");
		if($ref == "") $ref = "anova_pvals.txt";
		$valid = array();
		exec("grep ^'{$fstat}' {$ref}", $valid);
		for($k = 0; $k < count($valid); $k++) {
			$valid[$k] = rtrim($valid[$k]);
			$next = split("\t", $valid[$k]);
			if($next[0] == $fstat && $next[1] == $df2) return $next[2];
		}
		return "";
	}
	



}

?>