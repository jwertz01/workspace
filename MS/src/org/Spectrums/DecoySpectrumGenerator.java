package org.Spectrums;

import IO.GenericSpectrumReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import org.apache.commons.collections.MultiMap;
import org.apache.commons.collections.map.MultiValueMap;







public class DecoySpectrumGenerator
{
  public static char[] StaticResidues = { 'K', 'R', 'P' };
  private static String DECOY = "DECOY";
  private Spectrum s;
  public double tolerance = 0.05D;












































  
  public Spectrum generateDecoy(Spectrum s, Random rand) {
    this.s = s;
    Spectrum d = new Spectrum();
    d.spectrumName = s.spectrumName;
    d.scanNumber = s.scanNumber;
    d.charge = this.s.charge;
    d.parentMass = this.s.parentMass;
    d.rt = s.rt;
    Peptide pep = new Peptide(s.peptide, s.charge);
    d.peptide = shufflePeptide(pep, s.getPeptide(), rand);
    Peptide decoyPep = new Peptide(d.peptide, s.charge);
    TheoreticalSpectrum.prefixIons = new String[] { "b", "b-H20", "b-NH3" };
    TheoreticalSpectrum.suffixIons = new String[] { "y", "y-H20", "y-NH3" };
    TheoreticalSpectrum t = new TheoreticalSpectrum(pep);
    TheoreticalSpectrum.addIsotopicPeaks(t, 1);
    SimpleMatchingGraph g = t.getMatchGraph(s, this.tolerance);
    MultiMap map = generateAnnotationMap(g);
    TheoreticalSpectrum t2 = new TheoreticalSpectrum(decoyPep);
    TheoreticalSpectrum.addIsotopicPeaks(t2, 1);
    List<Peak> pList = new ArrayList<Peak>();
    Map<String, LabelledPeak> map2 = new HashMap<String, LabelledPeak>();
    for (Iterator<Peak> it = t2.getPeak().iterator(); it.hasNext(); ) {
      
      LabelledPeak lp = (LabelledPeak)it.next();
      String key = lp.getPeakType();
      map2.put(key, lp);
    } 
    for (Iterator<LabelledPeak> it = g.vertexSet(SimpleMatchingGraph.Observed).iterator(); it.hasNext(); ) {
      
      Peak p = (Peak)it.next();
      List neighs = g.getNeighbors(p);
      LabelledPeak closestP = null;
      double smallestErr = 10000.0D;
      for (int i = 0; i < neighs.size(); i++) {
        
        LabelledPeak currPeak = (LabelledPeak)neighs.get(i);
        double currDiff = Math.abs(currPeak.getMass() - p.getMass());
        closestP = (currDiff < smallestErr) ? currPeak : closestP;
        smallestErr = (currDiff < smallestErr) ? currDiff : smallestErr;
      } 
      if (closestP != null) {
        
        LabelledPeak lp = closestP;
        String key = lp.getPeakType();
        if (map2.containsKey(key)) {
          
          LabelledPeak newlp = (LabelledPeak)map2.get(key);
          Peak newPeak = new Peak(newlp.getMass(), p.getIntensity());
          pList.add(newPeak);
        } 
      } 
    } 
    //addUnannotatedPeaks(g, pList);
    
    Collections.sort(pList, PeakMassComparator.comparator);
    d.setPeaks(pList);
    d.spectrumName = "DECOY_" + d.spectrumName;
    d.protein = "DECOY_" + s.protein;
    return d;
  }


  
  public Spectrum generateTarget(Spectrum s) { return generateTarget(s, new String[] { "b", "b-H20", "b-NH3" }, new String[] { "y", "y-H20", "y-NH3" }); }


  
  public Spectrum generateTarget(Spectrum s, String[] prefixIons, String[] suffixIons) {
    this.s = s;
    Spectrum t = new Spectrum();
    t.spectrumName = s.spectrumName;
    t.scanNumber = s.scanNumber;
    t.charge = this.s.charge;
    t.parentMass = this.s.parentMass;
    t.peptide = s.peptide;
    t.protein = s.protein;
    t.rt = s.rt;
    Peptide pep = new Peptide(s.peptide, s.charge);
    
    TheoreticalSpectrum th = new TheoreticalSpectrum(pep, prefixIons, suffixIons);
    TheoreticalSpectrum.addIsotopicPeaks(th, 1);
    SimpleMatchingGraph g = th.getMatchGraph(s, this.tolerance);
    List<Peak> pList = new ArrayList<Peak>();
    for (Iterator<LabelledPeak> it = g.vertexSet(SimpleMatchingGraph.Observed).iterator(); it.hasNext(); ) {
      
      Peak p = (Peak)it.next();
      List neighs = g.getNeighbors(p);
      LabelledPeak closestP = null;
      double smallestErr = 10000.0D;
      for (int i = 0; i < neighs.size(); i++) {
        
        LabelledPeak currPeak = (LabelledPeak)neighs.get(i);
        double currDiff = Math.abs(currPeak.getMass() - p.getMass());
        closestP = (currDiff < smallestErr) ? currPeak : closestP;
        smallestErr = (currDiff < smallestErr) ? currDiff : smallestErr;
      } 
      neighs.size();
      
      neighs.size();
      if (closestP != null) {
        
        Peak theo = new Peak(closestP.getMass(), p.getIntensity());
        
        pList.add(theo);
      } 
    } 
    addUnannotatedPeaks(g, pList);
    
    Collections.sort(pList, PeakMassComparator.comparator);
    t.setPeaks(pList);
    return t;
  }

  
  private void addUnannotatedPeaks(SimpleMatchingGraph g, List<Peak> pList) {
    int annotated = pList.size();
    List<Peak> unAnnotated = new ArrayList<Peak>();
    for (Iterator<LabelledPeak> it = g.vertexSet(SimpleMatchingGraph.Observed).iterator(); it.hasNext(); ) {
      
      Peak p = (Peak)it.next();
      List neighs = g.getNeighbors(p);
      if (neighs.size() == 0) {
        pList.add(p);
      }
    } 
  }
















  
  private MultiMap generateAnnotationMap(TheoreticalSpectrum t, Spectrum s) { return generateAnnotationMap(t.getMatchGraph(s, 0.5D)); }


  
  private MultiMap generateAnnotationMap(SimpleMatchingGraph mg) {
	    MultiValueMap multiValueMap = new MultiValueMap();
	    for (Iterator<LabelledPeak> it = mg.vertexSet(SimpleMatchingGraph.Theoretical).iterator(); it.hasNext(); ) {
	      LabelledPeak lp = (LabelledPeak)it.next();
	      List neighs = mg.getNeighbors(lp);
	      for (Iterator<Peak> it2 = neighs.iterator(); it2.hasNext();) {
	        multiValueMap.put(lp.getPeakType(), it2.next());
	      }
	    } 
	    return multiValueMap;
   }

  
  public String shufflePeptide(Peptide pep, String modified_pep, Random rand) { return shuffle(pep, modified_pep, rand); }



  
  public String shuffle(Peptide pep_obj, String modified_pep, Random rand) {
	String pep = pep_obj.getPeptide();
	
    StringBuffer shuffle = new StringBuffer(pep);
    List<Integer> index = new ArrayList<Integer>();
    for (int i = 0; i < pep.length(); i++) {
      if (!isStaticResidue(pep.charAt(i))) {
        index.add(Integer.valueOf(i));
      }
    } 
    
    double nterm_mass = 0.0;
    
    // Shuffle corresponding mods along with amino acids
    int[] pos = new int[pep_obj.getPos().length];
    double[] masses = new double[pep_obj.getPtmmasses().length];
    for (int i = 0; i < pos.length; i++) {
    	// Is N-term
    	if (
    			(pep_obj.getPos()[i] == 1) &&
    			(modified_pep.startsWith("+") || modified_pep.startsWith("-") || modified_pep.startsWith("["))
    	) {
    		nterm_mass = pep_obj.getPtmmasses()[i];
    	}
    	else {
    		pos[i] = pep_obj.getPos()[i];
    		masses[i] = pep_obj.getPtmmasses()[i];
    	}
    }

    if (index.size() != 0) {
      for (int i = 0; i < 50; i++) {
        int j = rand.nextInt(index.size());
        int k = rand.nextInt(index.size());
        j = ((Integer)index.get(j)).intValue();
        k = ((Integer)index.get(k)).intValue();
        char c = shuffle.charAt(j);
     	for (int x = 0; x < pos.length; x++) {
 		  if (pos[x] - 1 == k) {
 			pos[x] = j + 1;
 		  }
 		  else if (pos[x] - 1 == j) {
 			pos[x] = k + 1;	
 	      }
   	    }
        shuffle.setCharAt(j, shuffle.charAt(k));
        shuffle.setCharAt(k, c);
        } 
    }
    
    // add mods to peptide string
    String outpepstr = "";
    if (nterm_mass > 0) {
    	outpepstr = "+" + Double.toString(nterm_mass);
    }
    else if (nterm_mass < 0) {
    	outpepstr = Double.toString(nterm_mass);
    }
    for (int i = 0; i < shuffle.length(); i++) {
    	int pos_index = -1;
    	for (int j = 0; j < pos.length; j++) {
    		if (i == pos[j] - 1) {
    			pos_index = j;
    		}
    	}
    	
    	if (pos_index == -1) {
    		outpepstr += shuffle.charAt(i);
    	}
    	else if (masses[pos_index] >= 0) {
    		outpepstr += shuffle.charAt(i);
    		outpepstr += "+" + masses[pos_index];
    	}
    	else {
    		outpepstr += shuffle.charAt(i);
    		outpepstr += masses[pos_index];
    	}
    }
    return outpepstr;
  }

  
  public boolean isStaticResidue(char c) {
    for (int i = 0; i < StaticResidues.length; i++) {
      if (c == StaticResidues[i]) {
        return true;
      }
    } 
    return false;
  }



























  
  public static void testShuffleSpectra(String spectrumFile, String decoyFileName, double tolerance) {
    Random rand = new Random(0L);
    GenericSpectrumReader genericSpectrumReader = new GenericSpectrumReader(spectrumFile);
    int index = spectrumFile.lastIndexOf('.');
    if (decoyFileName == null || decoyFileName.length() < 1) {
      decoyFileName = String.valueOf(spectrumFile.substring(0, index)) + "_plusDecoy.mgf";
    }
    
    try {
      BufferedWriter bo = new BufferedWriter(new FileWriter(decoyFileName));
      int fail = 0;
      while (genericSpectrumReader.hasNext()) {
        
        Spectrum s = (Spectrum)genericSpectrumReader.next();
        
        s.peptide.contains("+");
        
        DecoySpectrumGenerator generator = new DecoySpectrumGenerator();
        generator.tolerance = tolerance;
        
        Spectrum target = generator.generateTarget(s);
        if (target.getPeak().size() >= 10) {

          Spectrum decoy = generator.generateDecoy(s, rand);
          
          double sim = s.cosine(decoy, tolerance);
          int count = 0;
          String best_decoy = decoy.toString();
          double sim_best_decoy = sim;

          while (sim > 0.65D) {
            
            decoy = generator.generateDecoy(s, rand);
            sim = s.cosine(decoy, tolerance);
            count++;
            
            if (sim < sim_best_decoy) {
              best_decoy = decoy.toString();
          	  sim_best_decoy = sim;
            }
            
            if (count > 50) {
              /* changed Nov. 2019: Instead of failing if didn't generate a decoy with low 
              / enough cosine, use the best decoy that didn't pass .65 threshold*/
              fail++;          
              //decoy = null;
              break;
            } 
          }
          
          bo.write(target.toString());
            
          bo.write("\n");
            
          bo.write(best_decoy);
          bo.write("\n");
        } 
      } 
      //System.out.println("cannot generated decoy for " + fail + " spectra in the library");
      System.out.println(
    		  "Forced to create " + fail +
    		  " decoy library spectra that were very similar to target (cosine > .65)"
      );
      bo.flush();
      bo.close();
    }
    catch (IOException ioe) {
      
      System.out.println(ioe.getMessage());
    } 
  }





















































  
  public static void main(String[] args) { testShuffleSpectra(args[0], args[1], Double.parseDouble(args[2])); }
}
