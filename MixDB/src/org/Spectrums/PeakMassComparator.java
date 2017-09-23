package org.Spectrums;

import java.util.Comparator;

public class PeakMassComparator implements Comparator<Peak>{
	public static PeakMassComparator comparator = new PeakMassComparator();
	public int compare(Peak p0, Peak p1) {
		if(p0.getMass()> p1.getMass()){
			return 1;
		}else if(p0.getMass() == p1.getMass()){
			return 0;
		}else{
			return -1;
		}
	}
}
