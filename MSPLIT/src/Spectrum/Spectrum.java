package Spectrum;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Vector;

public class Spectrum
  implements Comparable<Spectrum>, Serializable
{
  public static final long serialVersionUID = 1L;
  public static double BINWIDTH = 1.0D;
  public static double MINMASS = 0.5D;
  private static double MAXMASS = 2000.0D;
  private static String GENERICNAME = "No-Spectrum-name";
  private static String GENERICPROTEIN = "No-Protein-name";
  private static String GENERICPEP = "DUMMYSEQ";
  private List<Peak> peaks;
  public String spectrumName = GENERICNAME;
  public String peptide = GENERICPEP;
  public String protein = GENERICPROTEIN;
  public double parentMass;
  public int charge;
  public double modMass;
  public int modPos;
  public double rt;
  double score = 0.0D;
  double upperBound = 0.0D;
  public int scanNumber = 0;
  public int specIndex = 0;
  
  public Spectrum()
  {
    this.peaks = new Vector();
    this.modMass = 0.0D;
    this.modPos = 0;
    this.charge = 1;
    this.parentMass = 10000.0D;
    
    this.scanNumber = 0;
  }
  
  public Spectrum(Vector peaks, String peptide, double pm, int charge, double modMass, int modPos)
  {
    this.peaks = peaks;
    this.peptide = peptide;
    this.parentMass = pm;
    this.charge = charge;
    this.modMass = modMass;
    this.modPos = modPos;
  }
  
  public Spectrum(Vector peaks, String peptide, String name, double pm, int charge, double modMass, int modPos, int scanNumber)
  {
    this.peaks = peaks;
    this.peptide = peptide;
    this.parentMass = pm;
    this.charge = charge;
    this.modMass = modMass;
    this.modPos = modPos;
    this.spectrumName = name;
    this.scanNumber = scanNumber;
  }
  
  public Spectrum(String spectrumFile, String format)
  {
    if (format.equals("MGF")) {
      readSpectrumFromMGF(spectrumFile);
    }
  }
  
  public Spectrum readSpectrumFromNIST(String fileName)
  {
    return new Spectrum();
  }
  
  public Spectrum(Spectrum s1, Spectrum s2, double scale1, double scale2, double tolerance)
  {
    s1 = s1.toNormVector(BINWIDTH, MINMASS, MAXMASS);
    s2 = s2.toNormVector(BINWIDTH, MINMASS, MAXMASS);
    
    this.peptide = (s1.peptide + " & " + s2.peptide);
    this.modMass = 0.0D;
    this.modPos = 0;
    this.charge = s1.charge;
    this.peaks = new Vector();
    this.parentMass = Math.max(s1.parentMass, s2.parentMass);
    
    this.modMass = (this.parentMass - Math.min(s1.parentMass, s2.parentMass));
    

    int i = 0;int j = 0;
    do
    {
      Peak p1 = (Peak)s1.peaks.get(i);
      Peak p2 = (Peak)s2.peaks.get(j);
      double mz1 = p1.getMass();
      double mz2 = p2.getMass();
      if (Math.abs(mz1 - mz2) <= tolerance)
      {
        this.peaks.add(new Peak(mz1, p1.getIntensity() * scale1 + p2.getIntensity() * scale2));
        i++;
        j++;
      }
      else if (mz1 < mz2)
      {
        this.peaks.add(new Peak(mz1, p1.getIntensity() * scale1));
        i++;
      }
      else
      {
        this.peaks.add(new Peak(mz2, p2.getIntensity() * scale2));
        j++;
      }
      if (i >= s1.peaks.size()) {
        break;
      }
    } while (j < s2.peaks.size());
    while (i < s1.peaks.size())
    {
      Peak p1 = (Peak)s1.peaks.get(i);
      this.peaks.add(new Peak(p1.getMass(), p1.getIntensity() * scale1));
      i++;
    }
    while (j < s2.peaks.size())
    {
      Peak p2 = (Peak)s2.peaks.get(j);
      this.peaks.add(new Peak(p2.getMass(), p2.getIntensity() * scale2));
      j++;
    }
  }
  
  public Spectrum(Spectrum s1, Spectrum s2, double scale1, double scale2, boolean toVector, double tolerance)
  {
    if (toVector)
    {
      s1 = s1.toVector(BINWIDTH, MINMASS, MAXMASS);
      s2 = s2.toVector(BINWIDTH, MINMASS, MAXMASS);
      s1.computePeakRank();
      s2.computePeakRank();
    }
    int commonpeak = 0;
    double mag1 = s1.sumOfPeaks();
    mag1 *= mag1;
    double mag2 = s2.sumOfPeaks();
    mag2 *= mag2;
    System.out.println(s1.peptide + " & " + s2.peptide + " real alpha " + mag1 / s1.magnitude() / (mag2 / s2.magnitude()));
    
    s1.scaleSpectrum(1.0D / mag1);
    s2.scaleSpectrum(1.0D / mag2);
    this.peptide = (s1.peptide + " & " + s2.peptide);
    this.modMass = 0.0D;
    this.modPos = 0;
    s1.charge += s2.charge;
    this.peaks = new Vector();
    this.parentMass = s1.parentMass;
    
    this.modMass = (this.parentMass - Math.min(s1.parentMass, s2.parentMass));
    
    int i = 0;int j = 0;
    do
    {
      Peak p1 = (Peak)s1.peaks.get(i);
      Peak p2 = (Peak)s2.peaks.get(j);
      double mz1 = p1.getMass();
      double mz2 = p2.getMass();
      if (Math.abs(mz1 - mz2) < tolerance)
      {
        this.peaks.add(new Peak(mz1, p1.getIntensity() * scale1 + p2.getIntensity() * scale2));
        if (p1.getIntensity() > p2.getIntensity()) {
          ((Peak)this.peaks.get(this.peaks.size() - 1)).copyRank(p1);
        } else {
          ((Peak)this.peaks.get(this.peaks.size() - 1)).copyRank(p2);
        }
        commonpeak++;
        i++;
        j++;
      }
      else if (mz1 < mz2)
      {
        this.peaks.add(new Peak(mz1, p1.getIntensity() * scale1));
        ((Peak)this.peaks.get(this.peaks.size() - 1)).copyRank(p1);
        i++;
      }
      else
      {
        this.peaks.add(new Peak(mz2, p2.getIntensity() * scale2));
        ((Peak)this.peaks.get(this.peaks.size() - 1)).copyRank(p2);
        j++;
      }
      if (i >= s1.peaks.size()) {
        break;
      }
    } while (j < s2.peaks.size());
    while (i < s1.peaks.size())
    {
      Peak p1 = (Peak)s1.peaks.get(i);
      this.peaks.add(new Peak(p1.getMass(), p1.getIntensity() * scale1));
      ((Peak)this.peaks.get(this.peaks.size() - 1)).copyRank(p1);
      i++;
    }
    while (j < s2.peaks.size())
    {
      Peak p2 = (Peak)s2.peaks.get(j);
      this.peaks.add(new Peak(p2.getMass(), p2.getIntensity() * scale2));
      ((Peak)this.peaks.get(this.peaks.size() - 1)).copyRank(p2);
      j++;
    }
  }
  
  public Spectrum(Spectrum s1, Spectrum s2, double scale1, double scale2, boolean toVector, boolean merge, double tolerance)
  {
    this(s1, s2, scale1, scale2, toVector, tolerance);
    if (merge) {
      mergePeaks(this, 0.25D);
    }
  }
  
  private void mergePeaks(Spectrum mix, double tolerance)
  {
    System.out.println("we start with  " + mix.peaks.size() + " peaks");
    for (int i = 0; i < mix.getPeaks().size(); i++)
    {
      Peak current = (Peak)mix.getPeaks().get(i);
      current.setIntensity(current.getIntensity() + Math.random() * 1.0E-7D);
    }
    TreeMap<Double, Peak> intensityMap = new TreeMap();
    TreeMap<Double, Peak> massMap = new TreeMap();
    List<Peak> newPeakList = new ArrayList();
    for (int i = 0; i < mix.getPeaks().size(); i++)
    {
      Peak current = (Peak)mix.getPeaks().get(i);
      intensityMap.put(Double.valueOf(current.getIntensity()), current);
      massMap.put(Double.valueOf(current.getMass()), current);
    }
    while (intensityMap.size() > 0)
    {
      Map.Entry<Double, Peak> e = intensityMap.lastEntry();
      Peak current = (Peak)e.getValue();
      Double currentKey = Double.valueOf(current.getIntensity());
      SortedMap<Double, Peak> subMap = massMap.subMap(Double.valueOf(current.getMass() - tolerance), Double.valueOf(current.getMass() + tolerance));
      
      List<Peak> toBeRemoved = new ArrayList();
      while (subMap.size() > 1)
      {
        for (Iterator<Peak> it = subMap.values().iterator(); it.hasNext();)
        {
          Peak neigh = (Peak)it.next();
          if (!neigh.equals(current))
          {
            double weight = current.getIntensity() / (current.getIntensity() + neigh.getIntensity());
            current.setMoz(current.getMass() * weight + neigh.getMass() * (1.0D - weight));
            current.setIntensity(current.getIntensity() + neigh.getIntensity());
            toBeRemoved.add(neigh);
          }
        }
        for (int j = 0; j < toBeRemoved.size(); j++)
        {
          intensityMap.remove(Double.valueOf(((Peak)toBeRemoved.get(j)).getIntensity()));
          massMap.remove(Double.valueOf(((Peak)toBeRemoved.get(j)).getMass()));
        }
        subMap = massMap.subMap(Double.valueOf(current.getMass() - tolerance), Double.valueOf(current.getMass() + tolerance));
      }
      newPeakList.add(current);
      intensityMap.remove(currentKey);
    }
    Collections.sort(newPeakList, PeakMassComparator.comparator);
    mix.setPeaks(newPeakList);
    System.out.println("after merging we have " + mix.peaks.size());
  }
  
  public Spectrum(Spectrum s1, Spectrum s2, double tolerance)
  {
    this(s1, s2, 1.0D, 1.0D, tolerance);
  }
  
  public void setPeaks(List<Peak> peaks)
  {
    this.peaks = peaks;
  }
  
  public List<Peak> getPeak()
  {
    return this.peaks;
  }
  
  public void readSpectrumFromMSP(String fileName)
  {
    try
    {
      BufferedReader bf = new BufferedReader(new FileReader(fileName));
      readSpectrumFromMSP(bf);
    }
    catch (IOException ioe)
    {
      System.out.println("Cannot Open MSP file");
      System.out.println(ioe.getMessage());
    }
  }
  
  public boolean readSpectrumFromMSP(BufferedReader bf)
  {
    try
    {
      boolean isPeaks = false;
      
      String line = bf.readLine();
      while ((line != null) && (line.length() != 0))
      {
        if (line.startsWith("Name:"))
        {
          this.charge = Integer.valueOf(line.split("[ ,/]")[2]).intValue();
          String pep = line.split("[ ,/]")[1];
          this.peptide = new String(pep);
        }
        else if (line.startsWith("Comment:"))
        {
          String temp = line.substring(line.indexOf("Protein"));
          
          String prot = temp.split("[=, ]")[1];
          this.protein = new String(this.protein);
          temp = line.substring(line.indexOf("Parent"));
          this.parentMass = Double.parseDouble(temp.split("[=, ]")[1]);
          
          temp = line.substring(line.indexOf("Mods"));
          String[] tokens = temp.split("[,,=, ,/]");
          
          int mod = Integer.parseInt(tokens[1]);
          if (mod == 0)
          {
            this.modPos = -1;
            this.modMass = 0.0D;
          }
          else
          {
            this.modPos = Integer.parseInt(tokens[2]);
            this.modMass = 100.0D;
          }
        }
        else if (line.startsWith("Num peaks:"))
        {
          int numPeaks = Integer.parseInt(line.split(": ")[1]);
          this.peaks = new ArrayList(numPeaks);
          isPeaks = true;
        }
        else if (isPeaks)
        {
          String[] tokens = line.split("\t");
          this.peaks.add(new Peak(Double.valueOf(tokens[0]).doubleValue(), 
            Double.valueOf(tokens[1]).doubleValue()));
        }
        line = bf.readLine();
      }
      if (line == null) {
        return false;
      }
    }
    catch (IOException ioe)
    {
      System.out.println("Cannot Open MSP file");
      System.out.println(ioe.getMessage());
      return false;
    }
    return true;
  }
  
  public boolean readSpectrumFromSplib(BufferedReader bf)
  {
    try
    {
      boolean isPeaks = false;
      
      String line = bf.readLine();
      while ((line != null) && (line.length() != 0))
      {
        if (line.startsWith("Name:"))
        {
          this.charge = Integer.valueOf(line.split("[ ,/]")[2]).intValue();
          String pep = line.split("[ ,/]")[1];
          this.peptide = new String(pep);
        }
        else if (line.startsWith("PrecursorMZ:"))
        {
          this.parentMass = Double.parseDouble(line.split(" ")[1]);
        }
        else if (line.startsWith("Comment:"))
        {
          String temp = line.substring(line.indexOf(" Protein="));
          String prot = temp.split("[=\\s+]")[2];
          if (prot.startsWith("\"")) {
            prot = prot.substring(1);
          }
          this.protein = new String(prot);
          if (line.contains("DECOY")) {
            this.protein = ("DECOY_" + this.protein);
          }
          temp = line.substring(line.indexOf("Mods"));
          String[] tokens = temp.split("[,,=, ,/]");
          
          int mod = Integer.parseInt(tokens[1]);
          if (mod == 0)
          {
            this.modPos = -1;
            this.modMass = 0.0D;
          }
          else
          {
            this.modPos = Integer.parseInt(tokens[2]);
            this.modMass = 100.0D;
          }
        }
        else if (line.startsWith("NumPeaks:"))
        {
          int numPeaks = Integer.parseInt(line.split(": ")[1]);
          isPeaks = true;
          this.peaks = new ArrayList(numPeaks);
        }
        else if (isPeaks)
        {
          String[] tokens = line.split("\t");
          this.peaks.add(new Peak(Double.valueOf(tokens[0]).doubleValue(), 
            Double.valueOf(tokens[1]).doubleValue()));
        }
        line = bf.readLine();
      }
      if (line == null) {
        return false;
      }
    }
    catch (IOException ioe)
    {
      System.out.println("Cannot Open splib file");
      System.out.println(ioe.getMessage());
      return false;
    }
    return true;
  }
  
  public boolean readSpectrumFromMS2(BufferedReader bf)
  {
    try
    {
      boolean isPeaks = false;
      
      String line = bf.readLine();
      while ((line != null) && ((!isPeaks) || (!line.startsWith("S"))))
      {
        String[] tokens = line.split("\\s+");
        if (line.startsWith("S"))
        {
          this.scanNumber = Integer.parseInt(tokens[1]);
          
          this.spectrumName = ("Scan Number: " + this.scanNumber);
          this.parentMass = Double.parseDouble(tokens[3]);
          this.charge = -1;
        }
        else if (line.startsWith("Z"))
        {
          this.charge = Integer.parseInt(tokens[1]);
          if (this.parentMass <= 0.0D) {
            this.parentMass = ((Double.parseDouble(tokens[2]) + Mass.PROTON_MASS * (this.charge - 1)) / this.charge);
          }
        }
        else if (line.startsWith("D"))
        {
          if (tokens[1].equals("seq")) {
            this.peptide = tokens[2];
          }
        }
        else if (Character.isDigit(line.charAt(0)))
        {
          isPeaks = true;
        }
        if (isPeaks)
        {
          tokens = line.split("\\s+");
          this.peaks.add(new Peak(Double.valueOf(tokens[0]).doubleValue(), 
            Double.valueOf(tokens[1]).doubleValue()));
        }
        bf.mark(200);
        line = bf.readLine();
      }
      bf.reset();
      if (line == null) {
        return false;
      }
    }
    catch (IOException ioe)
    {
      System.out.println("Cannot Open ms2 file");
      System.out.println(ioe.getMessage());
      return false;
    }
    return true;
  }
  
  public void readSpectrumFromMGF(String fileName)
  {
    try
    {
      BufferedReader bf = new BufferedReader(new FileReader(fileName));
      readSpectrumFromMGF(bf);
    }
    catch (IOException ioe)
    {
      System.out.println("Cannot Open MGF file");
      System.out.println(ioe.getMessage());
    }
  }
  
  public boolean readSpectrumFromMGF(BufferedReader bf)
  {
    try
    {
      String line;
      do
      {
        line = bf.readLine();
      } while ((line != null) && (!line.equals("BEGIN IONS")));
      boolean isPeaks = false;
      
      line = bf.readLine();
      while ((line != null) && (!line.contains("END IONS")))
      {
        if (line.length() >= 1) {
          if (line.startsWith("PEPMASS"))
          {
            String[] mass = line.split("=");
            mass = mass[1].split("\\s+");
            this.parentMass = Double.valueOf(mass[0]).doubleValue();
          }
          else if (line.startsWith("CHARGE"))
          {
            String charge = line.split("=")[1];
            if (charge.contains("+")) {
              charge = charge.replace("+", "");
            }
            this.charge = Integer.valueOf(charge).intValue();
          }
          else if ((line.startsWith("PEPSEQ")) || (line.startsWith("SEQ")))
          {
            this.peptide = line.split("=")[1];
          }
          else if (line.startsWith("PROTEIN"))
          {
            this.protein = line.split("=")[1];
          }
          else if (line.startsWith("TITLE"))
          {
            String[] tokens = line.split("=");
            if (tokens.length > 1)
            {
              this.spectrumName = line.split("=")[1];
              if (this.spectrumName.contains("Scan Number:"))
              {
                String[] tokens2 = this.spectrumName.split("\\s+");
                this.scanNumber = Integer.parseInt(tokens2[2]);
              }
              if (this.spectrumName.contains("DECOY")) {
                this.protein = ("DECOY_" + this.protein);
              }
            }
          }
          else if (line.startsWith("PEPMOD"))
          {
            this.modMass = 0.0D;
            this.modPos = 0;
          }
          else if (Character.isDigit(line.charAt(0)))
          {
            isPeaks = true;
          }
        }
        if (isPeaks)
        {
          String[] token = line.split("\\s+");
          this.peaks.add(new Peak(Double.valueOf(token[0]).doubleValue(), 
            Double.valueOf(token[1]).doubleValue()));
        }
        line = bf.readLine();
      }
      List<Peak> packed = new ArrayList(this.peaks.size());
      packed.addAll(this.peaks);
      this.peaks = packed;
      if (line == null) {
        return false;
      }
    }
    catch (IOException ioe)
    {
      System.out.println("Cannot Open MGF file");
      System.out.println(ioe.getMessage());
      return false;
    }
    return true;
  }
  
  public int compareTo(Spectrum s)
  {
    if (s.score == this.score) {
      return 0;
    }
    if (s.score < this.score) {
      return 1;
    }
    return -1;
  }
  
  public Spectrum toNormVector(double binWidth, double minMass, double maxMass)
  {
    Spectrum newSpectrum = toVector(binWidth, minMass, maxMass);
    newSpectrum.normalize();
    return newSpectrum;
  }
  
  public Spectrum toNormVector()
  {
    return toNormVector(BINWIDTH, MINMASS, MAXMASS);
  }
  
  public Spectrum toVector(double binWidth, double minMass, double maxMass)
  {
    Vector<Peak> newPeaks = new Vector();
    int i = 0;
    while (i < this.peaks.size()) {
        double mz = ((Peak)this.peaks.get(i)).getMass();
        double intensity = ((Peak)this.peaks.get(i)).getIntensity();
        newPeaks.add(new Peak(mz, intensity));
        i++;
    }

    Spectrum ret = new Spectrum(newPeaks, this.peptide, this.parentMass, this.charge, 
      this.modMass, this.modPos);
    ret.scanNumber = this.scanNumber;
    ret.spectrumName = this.spectrumName;
    ret.protein = this.protein;
    ret.specIndex = this.specIndex;
    return ret;
  }
  
  private void normalize()
  {
    double total = magnitude();
    scaleSpectrum(1.0D / total);
  }
  
  public void toRelIntensity()
  {
    double biggest = getMaxIntensity();
    scaleSpectrum(1.0D / biggest);
  }
  
  public String getPeptide()
  {
    return this.peptide;
  }
  
  public void setPeptide(String peptide)
  {
    this.peptide = peptide;
  }
  
  public List<Peak> getPeaks()
  {
    return this.peaks;
  }
  
  public double getMaxIntensity()
  {
    double max = 0.0D;
    for (int i = 0; i < this.peaks.size(); i++) {
      if (((Peak)this.peaks.get(i)).getIntensity() > max) {
        max = ((Peak)this.peaks.get(i)).getIntensity();
      }
    }
    return max;
  }
  
  public boolean checkOutLiner()
  {
    int i = 0;
    double total = magnitude();
    for (i = 0; i < this.peaks.size(); i++) {
      if (((Peak)this.peaks.get(i)).getIntensity() > 0.9D * total) {
        return true;
      }
    }
    return false;
  }
  
  protected double magnitude()
  {
    double total = 0.0D;
    for (int i = 0; i < this.peaks.size(); i++) {
      total += ((Peak)this.peaks.get(i)).getIntensity() * ((Peak)this.peaks.get(i)).getIntensity();
    }
    total = Math.pow(total, 0.5D);
    return total;
  }
  
  public double sumOfPeaks()
  {
    double total = 0.0D;
    for (int i = 0; i < this.peaks.size(); i++) {
      total += ((Peak)this.peaks.get(i)).getIntensity();
    }
    total = Math.pow(total, 0.5D);
    return total;
  }
  
  public double sumMagnitude()
  {
    double total = 0.0D;
    for (int i = 0; i < this.peaks.size(); i++) {
      total += ((Peak)this.peaks.get(i)).getIntensity();
    }
    return total;
  }
  
  public void scaleSpectrum(double factor)
  {
    for (int i = 0; i < this.peaks.size(); i++) {
      ((Peak)this.peaks.get(i)).scaleIntensity(factor);
    }
  }
  
  public void scaleMass(double factor)
  {
    for (int i = 0; i < this.peaks.size(); i++) {
      ((Peak)this.peaks.get(i)).scaleMass(factor);
    }
  }
  
  public void shiftSpectrum(double shift)
  {
    for (int i = 0; i < this.peaks.size(); i++) {
      ((Peak)this.peaks.get(i)).shiftMass(shift);
    }
  }
  
  public void shiftSpectrumPPM(double shiftPPM)
  {
    for (int i = 0; i < this.peaks.size(); i++) {
      ((Peak)this.peaks.get(i)).shiftMassPPM(shiftPPM);
    }
    this.parentMass += shiftPPM * this.parentMass / 1000000.0D;
  }
  
  public void sqrtSpectrum()
  {
    for (int i = 0; i < this.peaks.size(); i++) {
      ((Peak)this.peaks.get(i)).setIntensity(Math.pow(((Peak)this.peaks.get(i)).getIntensity(), 0.5D));
    }
  }
  
  public double cosine(Spectrum s1)
  {
    return 0.0D;
  }

  public double cosineSim1(Spectrum s1, double tolerance)
  {
      if (s1.peaks.size() == 0 || this.peaks.size() == 0) {
          return 0.0;
      }
      double magnitude = this.sumOfPeaks();
      magnitude *= s1.sumOfPeaks();
      double shift = this.parentMass - s1.parentMass;
      Set<List<Integer>> zeroShiftAlignments = this.findMatchPeaks(s1, 0, tolerance);
      Set<List<Integer>> realShiftAlignments;
      if (Math.abs(shift) > tolerance) {
          //realShiftAlignments = this.findMatchPeaks(s1, shift, tolerance); 
      }
      else {
          realShiftAlignments = new HashSet<List<Integer>>();
      }
      HashMap<Double, Set<List<Integer>>> pairsPerScore = new HashMap<Double, Set<List<Integer>>>();
      updatePairsPerScoreMap2(pairsPerScore, zeroShiftAlignments, s1);
      //updatePairsPerScoreMap2(pairsPerScore, realShiftAlignments, s1);
      ArrayList allScoresSorted = new ArrayList(new TreeSet(pairsPerScore.keySet()));
      Collections.reverse(allScoresSorted);
      Set<Integer> spec0PeakUsed = new HashSet<Integer>();
      Set<Integer> spec1PeakUsed = new HashSet<Integer>();
      double totalScore = 0.0;
      for (int i = 0; i < allScoresSorted.size(); i++) {
          double score = (double) allScoresSorted.get(i);
          Set<List<Integer>> scorePairs = pairsPerScore.get(score);
          Iterator<List<Integer>> it = scorePairs.iterator();
          while (it.hasNext()) {
              List<Integer> pair = it.next();
              int i0 = pair.get(0);
              int i1 = pair.get(1);
              if (!spec0PeakUsed.contains(i0) && !spec1PeakUsed.contains(i1)) {
                  spec0PeakUsed.add(i0);
                  spec1PeakUsed.add(i1);
                  totalScore += score;
              }
          }
      }
      return totalScore / magnitude;
  }  
  
  public double cosineSim2(Spectrum s1, double tolerance)
  {
      if (s1.peaks.size() == 0 || this.peaks.size() == 0) {
          return 0.0;
      }
      double magnitude = this.magnitude();
      magnitude *= s1.magnitude();
      double shift = this.parentMass - s1.parentMass;
      Set<List<Integer>> zeroShiftAlignments = this.findMatchPeaks(s1, 0, tolerance);
      Set<List<Integer>> realShiftAlignments;
      if (Math.abs(shift) > tolerance) {
          //realShiftAlignments = this.findMatchPeaks(s1, shift, tolerance); 
      }
      else {
          realShiftAlignments = new HashSet<List<Integer>>();
      }
      HashMap<Double, Set<List<Integer>>> pairsPerScore = new HashMap<Double, Set<List<Integer>>>();
      updatePairsPerScoreMap2(pairsPerScore, zeroShiftAlignments, s1);
      //updatePairsPerScoreMap2(pairsPerScore, realShiftAlignments, s1);
      ArrayList allScoresSorted = new ArrayList(new TreeSet(pairsPerScore.keySet()));
      Collections.reverse(allScoresSorted);
      Set<Integer> spec0PeakUsed = new HashSet<Integer>();
      Set<Integer> spec1PeakUsed = new HashSet<Integer>();
      double totalScore = 0.0;
      for (int i = 0; i < allScoresSorted.size(); i++) {
          double score = (double) allScoresSorted.get(i);
          Set<List<Integer>> scorePairs = pairsPerScore.get(score);
          Iterator<List<Integer>> it = scorePairs.iterator();
          while (it.hasNext()) {
              List<Integer> pair = it.next();
              int i0 = pair.get(0);
              int i1 = pair.get(1);
              if (!spec0PeakUsed.contains(i0) && !spec1PeakUsed.contains(i1)) {
                  spec0PeakUsed.add(i0);
                  spec1PeakUsed.add(i1);
                  totalScore += score;
              }
          }
      }
      return Math.pow(totalScore, .5) / magnitude;
  }
  
  public double shiftCosineSim(Spectrum s1, double fragmentTolerance)
  {
    Spectrum original = toNormVector();
    shiftSpectrum(0.5D);
    Spectrum rightShift = toNormVector();
    shiftSpectrum(-1.0D);
    Spectrum leftShift = toNormVector();
    shiftSpectrum(0.5D);
    Spectrum vS1 = s1.toNormVector();
    double cosineOriginal = original.cosineSim(vS1, fragmentTolerance);
    double cosineLeftShift = leftShift.cosineSim(vS1, fragmentTolerance);
    double cosineRightShift = rightShift.cosineSim(vS1, fragmentTolerance);
    
    return Math.max(cosineRightShift, 
      Math.max(cosineOriginal, cosineLeftShift));
  }
  
  public double shareSim(Spectrum s1, double fragmentTolerance)
  {
    double s = shareSim(s1, 0, fragmentTolerance);
    return s;
  }
  
  public double shareSim(Spectrum s1, int mode, double tolerance)
  {
      if (s1.peaks.size() == 0 || this.peaks.size() == 0) {
          return 0.0;
      }
      double product = 0;
      double magnitude = this.peaks.size();
      double magnitude2 = s1.peaks.size();
      double shift = this.parentMass - s1.parentMass;
      Set<List<Integer>> zeroShiftAlignments = this.findMatchPeaks(s1, 0, tolerance);
      Set<List<Integer>> realShiftAlignments;
      if (Math.abs(shift) > tolerance) {
          //realShiftAlignments = this.findMatchPeaks(s1, shift, tolerance); 
      }
      else {
          realShiftAlignments = new HashSet<List<Integer>>();
      }
      HashMap<Double, Set<List<Integer>>> pairsPerScore = new HashMap<Double, Set<List<Integer>>>();
      updatePairsPerScoreMap(pairsPerScore, zeroShiftAlignments, s1);
      //updatePairsPerScoreMap(pairsPerScore, realShiftAlignments, s1);
      ArrayList allScoresSorted = new ArrayList(new TreeSet(pairsPerScore.keySet()));
      Collections.reverse(allScoresSorted);
      Set<Integer> spec0PeakUsed = new HashSet<Integer>();
      Set<Integer> spec1PeakUsed = new HashSet<Integer>();
      double totalScore = 0.0;
      for (int i = 0; i < allScoresSorted.size(); i++) {
          double score = (double) allScoresSorted.get(i);
          Set<List<Integer>> scorePairs = pairsPerScore.get(score);
          Iterator<List<Integer>> it = scorePairs.iterator();
          while (it.hasNext()) {
              List<Integer> pair = it.next();
              int i0 = pair.get(0);
              int i1 = pair.get(1);
              if (!spec0PeakUsed.contains(i0) && !spec1PeakUsed.contains(i1)) {
                  spec0PeakUsed.add(i0);
                  spec1PeakUsed.add(i1);
                  totalScore += score;
                  product++;
              }
          }
      }
      if(mode == 0){
          return product/((magnitude+magnitude2)/2);
      }else if(mode == 1){
            return product/magnitude;
      }else{
            return product/magnitude2;
      }
  }
  
  public double projectedCosine(Spectrum s1, double tolerance)
  {
      if (s1.peaks.size() == 0 || this.peaks.size() == 0) {
          return 0.0;
      }
      double magnitude = this.magnitude();
      double projectedNorm = 0.00000001;
      double shift = this.parentMass - s1.parentMass;
      Set<List<Integer>> zeroShiftAlignments = this.findMatchPeaks(s1, 0, tolerance);
      Set<List<Integer>> realShiftAlignments;
      if (Math.abs(shift) > tolerance) {
          //realShiftAlignments = this.findMatchPeaks(s1, shift, tolerance); 
      }
      else {
          realShiftAlignments = new HashSet<List<Integer>>();
      }
      HashMap<Double, Set<List<Integer>>> pairsPerScore = new HashMap<Double, Set<List<Integer>>>();
      updatePairsPerScoreMap(pairsPerScore, zeroShiftAlignments, s1);
      //updatePairsPerScoreMap(pairsPerScore, realShiftAlignments, s1);
      ArrayList allScoresSorted = new ArrayList(new TreeSet(pairsPerScore.keySet()));
      Collections.reverse(allScoresSorted);
      Set<Integer> spec0PeakUsed = new HashSet<Integer>();
      Set<Integer> spec1PeakUsed = new HashSet<Integer>();
      double totalScore = 0.0;
      for (int i = 0; i < allScoresSorted.size(); i++) {
          double score = (double) allScoresSorted.get(i);
          Set<List<Integer>> scorePairs = pairsPerScore.get(score);
          Iterator<List<Integer>> it = scorePairs.iterator();
          while (it.hasNext()) {
              List<Integer> pair = it.next();
              int i0 = pair.get(0);
              int i1 = pair.get(1);
              if (!spec0PeakUsed.contains(i0) && !spec1PeakUsed.contains(i1)) {
                  spec0PeakUsed.add(i0);
                  spec1PeakUsed.add(i1);
                  totalScore += score;
                  projectedNorm += s1.peaks.get(i1).getIntensity() * s1.peaks.get(i1).getIntensity();
              }
          }
      }
      magnitude = magnitude * Math.pow(projectedNorm, 0.5);
      return totalScore / magnitude;
  }
  
  public double projectedCosine1(Spectrum s1, double tolerance)
  {
      if (s1.peaks.size() == 0 || this.peaks.size() == 0) {
          return 0.0;
      }
      double magnitude = this.sumOfPeaks();
      double projectedNorm = 0.00000001;
      double shift = this.parentMass - s1.parentMass;
      Set<List<Integer>> zeroShiftAlignments = this.findMatchPeaks(s1, 0, tolerance);
      Set<List<Integer>> realShiftAlignments;
      if (Math.abs(shift) > tolerance) {
          //realShiftAlignments = this.findMatchPeaks(s1, shift, tolerance); 
      }
      else {
          realShiftAlignments = new HashSet<List<Integer>>();
      }
      HashMap<Double, Set<List<Integer>>> pairsPerScore = new HashMap<Double, Set<List<Integer>>>();
      updatePairsPerScoreMap3(pairsPerScore, zeroShiftAlignments, s1);
      //updatePairsPerScoreMap3(pairsPerScore, realShiftAlignments, s1);
      ArrayList allScoresSorted = new ArrayList(new TreeSet(pairsPerScore.keySet()));
      Collections.reverse(allScoresSorted);
      Set<Integer> spec0PeakUsed = new HashSet<Integer>();
      Set<Integer> spec1PeakUsed = new HashSet<Integer>();
      double totalScore = 0.0;
      for (int i = 0; i < allScoresSorted.size(); i++) {
          double score = (double) allScoresSorted.get(i);
          Set<List<Integer>> scorePairs = pairsPerScore.get(score);
          Iterator<List<Integer>> it = scorePairs.iterator();
          while (it.hasNext()) {
              List<Integer> pair = it.next();
              int i0 = pair.get(0);
              int i1 = pair.get(1);
              if (!spec0PeakUsed.contains(i0) && !spec1PeakUsed.contains(i1)) {
                  spec0PeakUsed.add(i0);
                  spec1PeakUsed.add(i1);
                  totalScore += score;
                  projectedNorm += s1.peaks.get(i1).getIntensity();
              }
          }
      }
      magnitude = magnitude * Math.pow(projectedNorm, 0.5);
      return totalScore / magnitude;
  }
  
  public double alpha(Spectrum a, Spectrum b, double fragmentTolerance)
  {
    double C = dot(a, fragmentTolerance);
    double D = a.dot(b, fragmentTolerance);
    double E = dot(b, fragmentTolerance);
    double A = a.dot(a, fragmentTolerance);
    double B = b.dot(b, fragmentTolerance);
    double alpha = (B * C - D * E) / (A * E - C * D);
    if ((alpha < 0.1D) || (alpha > 10.0D) || (Double.isNaN(alpha))) {
      return 0.1D;
    }
    if (alpha > 1.0D) {
      return 1.0D / alpha;
    }
    return alpha;
  }
  
  public double maxScore(Spectrum a, Spectrum b, double alpha, double fragmentTolerance)
  {
    alpha = Math.pow(alpha, 0.5D);
    Spectrum answerMix1 = new Spectrum(a, b, alpha, 1.0D, fragmentTolerance);
    Spectrum answerMix2 = new Spectrum(b, a, alpha, 1.0D, fragmentTolerance);
    
    Spectrum baseAnswer = new Spectrum(a, b, 1.0D, 0.1D, fragmentTolerance);
    
    double score1 = cosineSim(answerMix1, fragmentTolerance);
    double score2 = cosineSim(answerMix2, fragmentTolerance);
    double score = cosineSim(baseAnswer, fragmentTolerance);
    if ((score1 > score2) && (score1 > score)) {
      return score1;
    }
    if (score2 > score) {
      return score2;
    }
    return score;
  }
  
  public double shiftMaxScore(Spectrum a, Spectrum b, double alpha, double fragmentTolerance)
  {
    Spectrum origA = a.toNormVector();
    a.shiftSpectrum(0.5D);
    Spectrum rightA = a.toNormVector();
    a.shiftSpectrum(-1.0D);
    Spectrum leftA = a.toNormVector();
    a.shiftSpectrum(0.5D);
    Spectrum origB = b.toNormVector();
    b.shiftSpectrum(0.5D);
    Spectrum rightB = b.toNormVector();
    b.shiftSpectrum(-1.0D);
    Spectrum leftB = b.toNormVector();
    b.shiftSpectrum(0.5D);
    Spectrum mix = toNormVector();
    
    double best = 0.0D;
    alpha = mix.alpha(origA, origB, fragmentTolerance);
    best = Math.max(best, mix.maxScore(origA, origB, alpha, fragmentTolerance));
    alpha = mix.alpha(origA, leftB, fragmentTolerance);
    best = Math.max(best, mix.maxScore(origA, leftB, alpha, fragmentTolerance));
    alpha = mix.alpha(origA, rightB, fragmentTolerance);
    best = Math.max(best, mix.maxScore(origA, rightB, alpha, fragmentTolerance));
    
    alpha = mix.alpha(leftA, origB, fragmentTolerance);
    best = Math.max(best, mix.maxScore(leftA, origB, alpha, fragmentTolerance));
    alpha = mix.alpha(leftA, leftB, fragmentTolerance);
    best = Math.max(best, mix.maxScore(leftA, leftB, alpha, fragmentTolerance));
    alpha = mix.alpha(leftA, rightB, fragmentTolerance);
    best = Math.max(best, mix.maxScore(leftA, rightB, alpha, fragmentTolerance));
    
    alpha = mix.alpha(rightA, origB, fragmentTolerance);
    best = Math.max(best, mix.maxScore(rightA, origB, alpha, fragmentTolerance));
    alpha = mix.alpha(rightA, leftB, fragmentTolerance);
    best = Math.max(best, mix.maxScore(rightA, leftB, alpha, fragmentTolerance));
    alpha = mix.alpha(rightA, rightB, fragmentTolerance);
    best = Math.max(best, mix.maxScore(rightA, rightB, alpha, fragmentTolerance));
    
    return best;
  }
  
  public double inverseMaxScore(Spectrum a, Spectrum b, double alpha, double fragmentTolerance)
  {
    Spectrum answerMix = new Spectrum(a, b, fragmentTolerance);
    Spectrum mix = toNormVector();
    
    mix.upscale(a, b, alpha);
    
    return mix.cosineSim(answerMix, fragmentTolerance);
  }
  
  public void upscale(Spectrum s1, Spectrum s2, double alpha)
  {
    int i = 0;int j = 0;int k = 0;
    while ((i < this.peaks.size()) && (j < s2.peaks.size()))
    {
      double mz1 = ((Peak)this.peaks.get(i)).getMass();
      double mz2 = ((Peak)s2.peaks.get(j)).getMass();
      double mz3;
      if (k >= s1.peaks.size()) {
        mz3 = 3000.0D;
      } else {
        mz3 = ((Peak)s1.peaks.get(k)).getMass();
      }
      if (mz1 < mz2)
      {
        i++;
      }
      else if (mz3 < mz2)
      {
        k++;
      }
      else if ((mz1 == mz2) && (mz3 > mz2))
      {
        if (((Peak)s2.peaks.get(j)).getIntensity() > 0.1D) {
          ((Peak)this.peaks.get(i)).scaleIntensity(alpha + 1.0D);
        }
        i++;
        j++;
      }
      else if ((mz1 == mz2) && (mz1 == mz3))
      {
        if (((Peak)s2.peaks.get(j)).getIntensity() > 0.1D)
        {
          double intensity = ((Peak)this.peaks.get(i)).getIntensity();
          intensity -= ((Peak)s1.peaks.get(k)).getIntensity();
          
          intensity *= (alpha + 1.0D);
          if (intensity > 0.0D) {
            ((Peak)this.peaks.get(i)).setIntensity(((Peak)this.peaks.get(i)).getIntensity() + intensity);
          }
        }
        i++;
        j++;
        k++;
      }
      else
      {
        j++;
      }
    }
  }
  
  public double dot(Spectrum s1, double tolerance)
  {
      if (s1.peaks.size() == 0 || this.peaks.size() == 0) {
          return 0.0;
      }
      double shift = this.parentMass - s1.parentMass;
      Set<List<Integer>> zeroShiftAlignments = this.findMatchPeaks(s1, 0, tolerance);
      Set<List<Integer>> realShiftAlignments;
      if (Math.abs(shift) > tolerance) {
          //realShiftAlignments = this.findMatchPeaks(s1, shift, tolerance); 
      }
      else {
          realShiftAlignments = new HashSet<List<Integer>>();
      }
      HashMap<Double, Set<List<Integer>>> pairsPerScore = new HashMap<Double, Set<List<Integer>>>();
      updatePairsPerScoreMap2(pairsPerScore, zeroShiftAlignments, s1);
      //updatePairsPerScoreMap2(pairsPerScore, realShiftAlignments, s1);
      ArrayList allScoresSorted = new ArrayList(new TreeSet(pairsPerScore.keySet()));
      Collections.reverse(allScoresSorted);
      Set<Integer> spec0PeakUsed = new HashSet<Integer>();
      Set<Integer> spec1PeakUsed = new HashSet<Integer>();
      double totalScore = 0.0;
      for (int i = 0; i < allScoresSorted.size(); i++) {
          double score = (double) allScoresSorted.get(i);
          Set<List<Integer>> scorePairs = pairsPerScore.get(score);
          Iterator<List<Integer>> it = scorePairs.iterator();
          while (it.hasNext()) {
              List<Integer> pair = it.next();
              int i0 = pair.get(0);
              int i1 = pair.get(1);
              if (!spec0PeakUsed.contains(i0) && !spec1PeakUsed.contains(i1)) {
                  spec0PeakUsed.add(i0);
                  spec1PeakUsed.add(i1);
                  totalScore += score;
              }
          }
      }
      return totalScore;
  }
  
  public void removeSharePeaks(Spectrum s)
  {
    Iterator<Peak> p1 = this.peaks.iterator();
    Iterator<Peak> p2 = s.getPeak().iterator();
    Peak peak1 = null;Peak peak2 = null;
    if ((p1.hasNext()) && (p2.hasNext()))
    {
      peak1 = (Peak)p1.next();
      peak2 = (Peak)p2.next();
    }
    while ((p1.hasNext()) && (p2.hasNext())) {
      if (peak1.getMass() < peak2.getMass())
      {
        peak1 = (Peak)p1.next();
      }
      else if (peak1.getMass() == peak2.getMass())
      {
        p1.remove();
        peak1 = (Peak)p1.next();
        peak2 = (Peak)p2.next();
      }
      else
      {
        peak2 = (Peak)p2.next();
      }
    }
    if (peak1.getMass() == peak2.getMass()) {
      p1.remove();
    }
  }
  
  public void removePrecursors(double tolerance)
  {
    List<Peak> toBeRemoved = new ArrayList();
    for (int i = 0; i < this.peaks.size(); i++)
    {
      Peak p = (Peak)this.peaks.get(i);
      if (Math.abs(p.getMass() - this.parentMass) < tolerance)
      {
        System.out.println("removing " + p);
        toBeRemoved.add(p);
      }
      if (Math.abs(p.getMass() - (this.parentMass - Mass.WATER / this.charge)) < tolerance)
      {
        System.out.println("removing " + p);
        toBeRemoved.add(p);
      }
      if (Math.abs(p.getMass() - (this.parentMass - Mass.NH3 / this.charge)) < tolerance)
      {
        System.out.println("removing " + p);
        toBeRemoved.add(p);
      }
      if (Math.abs(p.getMass() - (this.parentMass - 2.0D * Mass.WATER / this.charge)) < tolerance)
      {
        System.out.println("removing " + p);
        toBeRemoved.add(p);
      }
      if (Math.abs(p.getMass() - (this.parentMass - (Mass.WATER + Mass.NH3) / this.charge)) < tolerance)
      {
        System.out.println("removing " + p);
        toBeRemoved.add(p);
      }
    }
    this.peaks.removeAll(toBeRemoved);
  }
  
  public double residual(Spectrum s, double tolerance)
  {
        Iterator<Peak> p1 = this.peaks.iterator();
        Iterator<Peak> p2 = s.getPeak().iterator();
        Peak peak1 = null, peak2 = null;
        int exp = 4;
//      if(this.isSqrtTrans){
//          exp = 4;
//      }else{
//          exp = 2;
//      }
        double shareIntensity = 0;
        double residIntensity = 0.001; //avoid div-by-zero 
        if(p1.hasNext() && p2.hasNext()){
            peak1 = p1.next();
            peak2 = p2.next();
        }else{
            return shareIntensity/residIntensity;
        }
        double remain;
        //System.out.println("magnitude is: " + this.magnitude());
        while(p1.hasNext() && p2.hasNext()){
            if(Math.abs(peak1.getMass() - peak2.getMass()) <= tolerance){
                remain = peak1.getIntensity() - peak2.getIntensity();
                if(remain < 0){
                    remain = 0;
                }
                shareIntensity += Math.pow(peak1.getIntensity(),exp);
                //System.out.println("share: " + peak2.getMass() + "\t" + peak2.getIntensity());
                //System.out.println("remaining: " + peak1.getMass() + "\t" + remain);
                //residIntensity += Math.pow(remain,exp);
                peak1 = p1.next();
                peak2 = p2.next();
            }
            else if(peak1.getMass() < peak2.getMass()){
                residIntensity += Math.pow(peak1.getIntensity(),exp);
                //System.out.println("res intensity: " + peak1.getMass() + "\t" + peak1.getIntensity());
                peak1 = p1.next();
            }else{
                //shareIntensity += Math.pow(peak1.getIntensity(), 2);
                peak2 = p2.next();
            }
        }
        
        if(Math.abs(peak1.getMass() - peak2.getMass()) < tolerance){ //take care of last element
                shareIntensity += Math.pow(peak1.getIntensity(), exp);
                remain = peak1.getIntensity() - peak2.getIntensity();
                if(remain < 0){
                    remain = 0;
                }
                //residIntensity += Math.pow(remain, 2);
        }
        if(!p1.hasNext()){
            residIntensity += Math.pow(peak1.getIntensity(), exp);
        }
        
//      if(!p2.hasNext()){
//          shareIntensity += Math.pow(peak2.getIntensity(), 2);
//          //System.out.println("share: " + peak2.getMass() + "\t" + peak2.getIntensity());
//      }
        
        while(p1.hasNext()){
            residIntensity += Math.pow(p1.next().getIntensity(), exp);
        }
        
//      while(p2.hasNext()){
//          shareIntensity += Math.pow(p2.next().getIntensity(), 2);
//      }
            
        //System.out.println("shareIntensity: " + Math.pow(shareIntensity,0.5));
//      System.out.println("s has intensity: " + s.magnitude());
        //System.out.println("resIntensity: " + Math.pow(residIntensity,0.5));
        //double alpha = Math.pow(shareIntensity/residIntensity, 0.5); //initial estimate
        double alpha = Math.pow(residIntensity, 0.5);
        //System.out.println("alpha: " + alpha);
        alpha = alpha*alpha / (1-alpha*alpha);  //reestimate by taking into account that mixture is normalized to one
        alpha =  Math.pow(alpha, 0.5); //thus each component is down-weighted slightly
        if(alpha > 1.0){
            return 1/alpha;
        }else{
            return alpha;
        }
//return alpha;
  }
  
  public void subtractSharePeaks(Spectrum s, double tolerance)
  {
    Iterator<Peak> p1 = this.peaks.iterator();
    Iterator<Peak> p2 = s.getPeak().iterator();
    Peak peak1 = null;Peak peak2 = null;
    if ((p1.hasNext()) && (p2.hasNext()))
    {
      peak1 = (Peak)p1.next();
      peak2 = (Peak)p2.next();
    }
    while ((p1.hasNext()) && (p2.hasNext())) {
        if (Math.abs(peak1.getMass() - peak2.getMass()) <= tolerance)
        {
          peak1.setIntensity(peak1.getIntensity() - peak2.getIntensity());
          if (peak1.getIntensity() < 0.0D) {
            p1.remove();
          }
          peak1 = (Peak)p1.next();
          peak2 = (Peak)p2.next();
        }
       else if (peak1.getMass() < peak2.getMass())
      {
        peak1 = (Peak)p1.next();
      }
      else
      {
        peak2 = (Peak)p2.next();
      }
    }
    if (Math.abs(peak1.getMass() - peak2.getMass()) <= tolerance) {
      p1.remove();
    }
  }
  
  public void filterPeaks(int n)
  {
    Vector<Peak> sortedPeakList = new Vector();
    sortedPeakList.addAll(this.peaks);
    Collections.sort(sortedPeakList, new peakComparator(null));
    int i = 0;
    for (i = 0; i < sortedPeakList.size() - n - 1; i++) {
      this.peaks.remove(sortedPeakList.get(i));
    }
  }
  
  public int explainedIntensity(double percent)
  {
    Vector<Peak> sortedPeakList = new Vector();
    sortedPeakList.addAll(this.peaks);
    Collections.sort(sortedPeakList, new peakComparator(null));
    double mag = magnitude();
    mag *= mag;
    int i = 0;
    
    double current = 0.0D;
    for (i = sortedPeakList.size() - 1; i > 0; i--)
    {
      current = current + ((Peak)sortedPeakList.get(i)).getIntensity() * ((Peak)sortedPeakList.get(i)).getIntensity();
      if (current / mag > percent) {
        return sortedPeakList.size() - i;
      }
    }
    return sortedPeakList.size();
  }
  
  public int intensePeakCount(double threshold)
  {
    int count = 0;
    for (int i = 0; i < this.peaks.size(); i++) {
      if (((Peak)this.peaks.get(i)).getIntensity() > threshold) {
        count++;
      }
    }
    return count;
  }
  
  public void windowFilterPeaks(int n, double deltaM)
  {
    int i = 0;
    
    Vector<Peak> toBeRemove = new Vector();
    Vector<Peak> neighbors = new Vector();
    for (int j = 0; j < this.peaks.size(); j++)
    {
      Peak p = (Peak)this.peaks.get(j);
      
      removeNonNeighborPeaks(deltaM, neighbors, p);
      
      i = addNeighborPeaks(deltaM, neighbors, p, i);
      
      Vector sortedNeighbors = new Vector(neighbors);
      Collections.sort(sortedNeighbors, new peakComparator(null));
      if (getPeakRank(sortedNeighbors, p) > n) {
        toBeRemove.add(p);
      }
    }
    this.peaks.removeAll(toBeRemove);
  }
  
  public void windowFilterPeaks2(int n, double deltaM)
  {
    int current = 0;int left = 0;int right = 0;
    Collection<Peak> toBeRemove = new ArrayList();
    TreeSet<Peak> neighs = new TreeSet(PeakIntensityComparator.comparator);
    while (current < this.peaks.size())
    {
      Peak p = (Peak)this.peaks.get(current);
      for (int i = left; i < current; i++)
      {
        Peak smaller = (Peak)this.peaks.get(i);
        if (p.getMass() - smaller.getMass() > deltaM)
        {
          neighs.remove(smaller);
          left = i;
        }
      }
      for (int i = right + 1; i < this.peaks.size(); i++)
      {
        Peak bigger = (Peak)this.peaks.get(i);
        if (bigger.getMass() - p.getMass() <= deltaM)
        {
          neighs.add(bigger);
        }
        else
        {
          right = i - 1;
          break;
        }
      }
      Iterator<Peak> it = neighs.descendingIterator();
      for (int i = 0; it.hasNext(); i++)
      {
        Peak p2 = (Peak)it.next();
        if ((p2 == p) && (i > n)) {
          toBeRemove.add(p);
        }
      }
      System.out.println();
      current++;
    }
    this.peaks.removeAll(toBeRemove);
  }
  
  public void windowFilterAndRank(int n, double deltaM, int topN)
  {
    int i = 0;
    
    Vector<Peak> toBeRemove = new Vector();
    Vector<Peak> neighbors = new Vector();
    for (int j = 0; j < this.peaks.size(); j++)
    {
      Peak p = (Peak)this.peaks.get(j);
      
      removeNonNeighborPeaks(deltaM, neighbors, p);
      
      i = addNeighborPeaks(deltaM, neighbors, p, i);
      
      Vector sortedNeighbors = new Vector(neighbors);
      Collections.sort(sortedNeighbors, new peakComparator(null));
      if (getPeakRank(sortedNeighbors, p) > n) {
        toBeRemove.add(p);
      }
    }
    this.peaks.removeAll(toBeRemove);
    Collections.sort(toBeRemove, PeakIntensityComparator.comparator);
    computePeakRank();
    int count = 1;
    int toBeInsert = topN - this.peaks.size() + 1;
    for (int k = toBeRemove.size() - 1; (k > toBeRemove.size() - toBeInsert) && (k > 0); k--)
    {
      ((Peak)toBeRemove.get(k)).setRank(this.peaks.size() + count);
      this.peaks.add((Peak)toBeRemove.get(k));
      count++;
    }
    Collections.sort(this.peaks, PeakMassComparator.comparator);
  }
  
  public List<Peak> topWindowPeaks(int n, double deltaM)
  {
    int i = 0;
    
    List<Peak> toBeKept = new Vector();
    List<Peak> neighbors = new Vector();
    for (int j = 0; j < this.peaks.size(); j++)
    {
      Peak p = (Peak)this.peaks.get(j);
      
      removeNonNeighborPeaks(deltaM, neighbors, p);
      
      i = addNeighborPeaks(deltaM, neighbors, p, i);
      
      List<Peak> sortedNeighbors = new Vector(neighbors);
      Collections.sort(sortedNeighbors, new peakComparator(null));
      if (getPeakRank(sortedNeighbors, p) <= n) {
        toBeKept.add(p);
      }
    }
    System.out.println("We start with " + this.peaks.size() + " peaks");
    System.out.println("After window-filtering we have: " + toBeKept.size() + " peaks");
    return toBeKept;
  }
  
  public Map<Peak, Integer> getPeakRank()
  {
    List<Peak> sortedList = new Vector();
    sortedList.addAll(this.peaks);
    Collections.sort(sortedList, new peakComparator(null));
    Map<Peak, Integer> peakRankMap = new HashMap();
    
    int i = 0;
    for (int size = sortedList.size(); i < size; i++)
    {
      Peak p = (Peak)sortedList.get(i);
      peakRankMap.put(p, new Integer(i));
    }
    return peakRankMap;
  }
  
  public void computePeakRank()
  {
    List<Peak> sortedList = new Vector();
    sortedList.addAll(this.peaks);
    Collections.sort(sortedList, new peakComparator(null));
    Map<Peak, Peak> peakMap = new HashMap();
    
    int i = 0;
    for (int size = sortedList.size(); i < size; i++)
    {
      Peak p = (Peak)sortedList.get(i);
      p.setRank(size - i);
    }
  }
  
  public List<Peak> getSortedPeaks()
  {
    List<Peak> sortedList = new Vector();
    sortedList.addAll(this.peaks);
    Collections.sort(sortedList, new peakComparator(null));
    return sortedList;
  }
  
  public List<Peak> getTopPeaks(int N)
  {
    List<Peak> sortedList = getSortedPeaks();
    int begin = sortedList.size() - N;
    begin = begin < 0 ? 0 : begin;
    return sortedList.subList(begin, sortedList.size());
  }
  
  public List<Peak> getTopPeaks(int N, double minMassDiff)
  {
    List<Peak> sortedList = getSortedPeaks();
    List<Peak> ret = new ArrayList();
    int addCount = 0;
    int index = sortedList.size() - 1;
    while ((index >= 0) && (addCount < N))
    {
      Peak toBeAdd = (Peak)sortedList.get(index);
      boolean canAdd = true;
      for (int i = 0; i < ret.size(); i++) {
        if (Math.abs(((Peak)ret.get(i)).getMass() - toBeAdd.getMass()) < minMassDiff)
        {
          canAdd = false;
          break;
        }
      }
      if (canAdd)
      {
        ret.add(toBeAdd);
        addCount++;
      }
      index--;
    }
    return ret;
  }
  
  public List<Peak> getTopPeaks(int N, List<Peak> exclusionList)
  {
    List<Peak> sortedList = getSortedPeaks();
    List<Peak> out = new ArrayList();
    for (int i = sortedList.size() - N; (i > 0) && (out.size() < N); i--)
    {
      Peak p = (Peak)sortedList.get(i);
      if (!exclusionList.contains(p)) {
        out.add(p);
      }
    }
    return out;
  }
  
  private int addNeighborPeaks(double deltaM, Collection<Peak> v, Peak p, int beginInd)
  {
    boolean beginFlag = false;
    for (int i = beginInd; i < this.peaks.size(); i++) {
      if (Math.abs(((Peak)this.peaks.get(i)).getMass() - p.getMass()) <= deltaM)
      {
        v.add((Peak)this.peaks.get(i));
        beginFlag = true;
      }
      else
      {
        if (beginFlag) {
          return i;
        }
        if ((beginFlag) && (i == this.peaks.size() - 1)) {
          return i + 1;
        }
      }
    }
    return beginInd;
  }
  
  private void removeNonNeighborPeaks(double deltaM, Collection<Peak> v, Peak p)
  {
    Iterator<Peak> it = v.iterator();
    while (it.hasNext()) {
      if (Math.abs(((Peak)it.next()).getMass() - p.getMass()) > deltaM) {
        it.remove();
      }
    }
  }
  
  private int getPeakRank(List<Peak> l, Peak p)
  {
    for (int i = 0; i < l.size(); i++) {
      if (((Peak)l.get(i)).equals(p)) {
        return l.size() - i;
      }
    }
    return -1;
  }
  
  public int numberOfPeaks(double percent)
  {
    Vector<Peak> sortedPeakList = new Vector();
    sortedPeakList.addAll(this.peaks);
    Collections.sort(sortedPeakList, new peakComparator(null));
    
    double total = 0.0D;
    int i = 0;int count = 0;
    for (i = sortedPeakList.size() - 1; i >= 0; i--) {
      total += ((Peak)sortedPeakList.get(i)).getIntensity();
    }
    double current = 0.0D;
    for (i = sortedPeakList.size() - 1; i >= 0; i--)
    {
      current += ((Peak)sortedPeakList.get(i)).getIntensity();
      if (current / total >= percent)
      {
        count++;return count;
      }
      count++;
    }
    return count;
  }
  
  public String toString()
  {
    StringBuffer sb = new StringBuffer();
    sb.append("BEGIN IONS\n");
    
    sb.append("TITLE=" + this.spectrumName + "\n");
    
    sb.append("CHARGE=" + this.charge + "\n");
    sb.append("PEPMASS=" + this.parentMass + "\n");
    if ((this.peptide != null) && (this.peptide.contains(" & "))) {
      sb.append("SEQ=" + this.peptide + "\n");
    } else {
      sb.append("SEQ=" + this.peptide + "\n");
    }
    if (this.scanNumber > 0) {
      sb.append("SCAN=" + this.scanNumber);
    }
    if (this.specIndex > 0) {
      sb.append("SPECINDEX=" + this.specIndex);
    }
    for (int i = 0; i < this.peaks.size(); i++) {
      sb.append(this.peaks.get(i) + "\n");
    }
    sb.append("END IONS\n");
    return sb.toString();
  }
  
  public static void main(String[] args)
  {
    testReadMGF();
    testtoVector();
    testVectorRep();
    testMixSpect(.5);
    testRemoveSharePeak();
    testPeakCount();
    
    testWindowFilterPeak();
    testShiftCosine();
  }
  
  public static void testReadMGF()
  {
    String filename = "testspectrum.mgf";
    Spectrum msms = new Spectrum();
    msms.readSpectrumFromMGF(filename);
    System.out.print(msms);
  }
  
  public static void testtoVector()
  {
    System.out.println("generating vectored spectrum");
    String filename = "testspectrum.mgf";
    Spectrum msms = new Spectrum();
    msms.readSpectrumFromMGF(filename);
    Spectrum s = msms.toVector(50.0D, 100.0D, 500.0D);
    System.out.println(s);
  }
  
  public static void testVectorRep()
  {
    String filename = "testspectrum.mgf";
    Spectrum msms = new Spectrum();
    msms.readSpectrumFromMGF(filename);
    Spectrum msms2 = new Spectrum();
    msms2.readSpectrumFromMGF(filename);
    System.out.println("Comparing self to self: " + 
      msms2.cosineSim(msms, .5));
    
    msms = msms.toNormVector(0.5D, 100.0D, 500.0D);
    msms2 = msms2.toNormVector(0.5D, 100.0D, 500.0D);
    System.out.println("Comparing self to self: " + 
      msms2.cosineSim(msms, .5));
    System.out.println("scaling one copy of self");
    msms.scaleSpectrum(0.03D);
    System.out.println("Comparing self to self: " + 
      msms2.cosineSim(msms, .5));
    System.out.println();
  }
  
  public static void testMixSpect(double fragmentTolerance)
  {
    System.out.println("generating mix spectrum");
    String filename = "testspectrum.mgf";
    Spectrum s1 = new Spectrum();
    s1.readSpectrumFromMGF("testspectrum.mgf");
    Spectrum s2 = new Spectrum();
    s2.readSpectrumFromMGF("testspectrum2.mgf");
    Spectrum msms12 = new Spectrum(s1, s2, fragmentTolerance);
    s1 = s1.toNormVector();
    s2 = s2.toNormVector();
    System.out.println(s1);
    System.out.println(s2);
    System.out.println(msms12);
    System.out.println("cosine: " + s2.cosineSim(msms12, .5));
    System.out.println("projected: " + s2.projectedCosine(msms12, fragmentTolerance));
    System.out.println("cosine: " + s1.cosineSim(msms12, .5));
    System.out.println("projected: " + s1.projectedCosine(msms12, fragmentTolerance));
    System.out.println("distinct: " + s1.cosineSim(s2, .5));
    System.out.println("distinct projected: " + s1.projectedCosine(s2, fragmentTolerance));
  }
  
  public static void testRemoveSharePeak()
  {
    String filename = "testspectrum.mgf";
    Spectrum msms = new Spectrum();
    msms.readSpectrumFromMGF(filename);
    Spectrum msms2 = new Spectrum();
    msms2.readSpectrumFromMGF(filename);
    System.out.println("original spectrum");
    System.out.println(msms);
    System.out.println("removing self against self");
    msms.removeSharePeaks(msms2);
    System.out.println(msms);
    msms.readSpectrumFromMGF(filename);
    ((Peak)msms.peaks.get(0)).setMoz(20.0D);
    System.out.println(msms);
    System.out.println("removing all but the first peak");
    msms.removeSharePeaks(msms2);
    System.out.println(msms);
  }
  
  public static void testPeakCount()
  {
    String filename = "testspectrum.mgf";
    Spectrum msms = new Spectrum();
    msms.readSpectrumFromMGF(filename);
    System.out.println(msms);
    System.out.println("total magnitude: " + msms.magnitude());
    for (double i = 0.0D; i < 11.0D; i += 1.0D) {
      System.out.println("number of peaks at" + i * 10.0D + "%: " + msms.numberOfPeaks(i / 10.0D));
    }
  }
  
  public static void testfilterPeak()
  {
    String filename = "testspectrum.mgf";
    Spectrum msms = new Spectrum();
    msms.readSpectrumFromMGF(filename);
    System.out.println(msms);
    System.out.println("total magnitude: " + msms.magnitude());
    for (int i = 20; i > 0; i--)
    {
      msms.filterPeaks(i);
      System.out.println(msms);
    }
  }
  
  public static void testWindowFilterPeak()
  {
    String filename = "testspectrum.mgf";
    Spectrum msms = new Spectrum();
    msms.readSpectrumFromMGF(filename);
    System.out.println("before window filter: ");
    System.out.println(msms);
    msms.windowFilterPeaks(2, 50.0D);
    System.out.println("after window filter: ");
    System.out.println(msms);
  }
  
  public static void testShiftCosine()
  {
    double fragmentTolerance = .5;
    String filename = ".\\mixture_compressed\\new80min.mgf";
    String fileMix = ".\\mixture_compressed\\new2min.mgf";
    SpectrumLib lib1 = new SpectrumLib(filename, "MGF");
    SpectrumLib mixlib = new SpectrumLib(fileMix, "MGF");
    Spectrum s = (Spectrum)lib1.getSpectra("spec_3182.dta..1").get(0);
    Spectrum m = (Spectrum)mixlib.getSpectra("spec_186.dta..1").get(0);
    Spectrum sCopy = s;
    System.out.println("shift cosine: " + s.shiftCosineSim(m, fragmentTolerance));
    System.out.println("shift cosine: " + s.shiftCosineSim(m, fragmentTolerance));
    System.out.println("shift cosine: " + m.shiftCosineSim(s, fragmentTolerance));
    System.out.println("");
    System.out.println("shift cosine with self: " + s.shiftCosineSim(sCopy, fragmentTolerance));
    System.out.println("shift cosine with self: " + sCopy.shiftCosineSim(s, fragmentTolerance));
    System.out.println("");
    s = s.toNormVector();
    sCopy = sCopy.toNormVector();
    System.out.println("cosine: " + s.cosineSim(sCopy, fragmentTolerance));
    System.out.println("cosine: " + sCopy.cosineSim(s, fragmentTolerance));
  }
  
  private class peakComparator
    implements Comparator<Peak>
  {
    private peakComparator(Object object) {}
    
    public int compare(Peak p0, Peak p1)
    {
      if (p0.getIntensity() > p1.getIntensity()) {
        return 1;
      }
      if (p0.getIntensity() == p1.getIntensity()) {
        return 0;
      }
      return -1;
    }
  }
  
  public double cosineSim(Spectrum s1, double tolerance) {
      if (s1.peaks.size() == 0 || this.peaks.size() == 0) {
          return 0.0;
      }
      double magnitude = this.magnitude();
      magnitude *= s1.magnitude();
      double shift = this.parentMass - s1.parentMass;
      Set<List<Integer>> zeroShiftAlignments = this.findMatchPeaks(s1, 0, tolerance);
      Set<List<Integer>> realShiftAlignments;
      if (Math.abs(shift) > tolerance) {
          //realShiftAlignments = this.findMatchPeaks(s1, shift, tolerance); 
      }
      else {
          realShiftAlignments = new HashSet<List<Integer>>();
      }
      HashMap<Double, Set<List<Integer>>> pairsPerScore = new HashMap<Double, Set<List<Integer>>>();
      updatePairsPerScoreMap(pairsPerScore, zeroShiftAlignments, s1);
      //updatePairsPerScoreMap(pairsPerScore, realShiftAlignments, s1);
      ArrayList allScoresSorted = new ArrayList(new TreeSet(pairsPerScore.keySet()));
      Collections.reverse(allScoresSorted);
      Set<Integer> spec0PeakUsed = new HashSet<Integer>();
      Set<Integer> spec1PeakUsed = new HashSet<Integer>();
      double totalScore = 0.0;
      for (int i = 0; i < allScoresSorted.size(); i++) {
          double score = (double) allScoresSorted.get(i);
          Set<List<Integer>> scorePairs = pairsPerScore.get(score);
          Iterator<List<Integer>> it = scorePairs.iterator();
          while (it.hasNext()) {
              List<Integer> pair = it.next();
              int i0 = pair.get(0);
              int i1 = pair.get(1);
              if (!spec0PeakUsed.contains(i0) && !spec1PeakUsed.contains(i1)) {
                  spec0PeakUsed.add(i0);
                  spec1PeakUsed.add(i1);
                  totalScore += score;
              }
          }
      }
      return totalScore / magnitude;
  }
  
  public void updatePairsPerScoreMap(
          HashMap<Double, Set<List<Integer>>> pairsPerScore, Set<List<Integer>> alignments, Spectrum s1
  ) {
      Iterator<List<Integer>> it = alignments.iterator();
      while(it.hasNext()) {
          List<Integer> l = it.next();
          int i0 = l.get(0);
          int i1 = l.get(1);
          double score = this.peaks.get(i0).getIntensity() * s1.peaks.get(i1).getIntensity();
          Set<List<Integer>> currVal;
          if (pairsPerScore.containsKey(score)) {
              currVal = pairsPerScore.get(score);
          }
          else {
              currVal = new HashSet<List<Integer>>();
          }
          currVal.add(l);
          pairsPerScore.put(score, currVal);
      }
  }
  
  public void updatePairsPerScoreMap1(
          HashMap<Double, Set<List<Integer>>> pairsPerScore, Set<List<Integer>> alignments, Spectrum s1
  ) {
      Iterator<List<Integer>> it = alignments.iterator();
      while(it.hasNext()) {
          List<Integer> l = it.next();
          int i0 = l.get(0);
          int i1 = l.get(1);
          double score = Math.pow(this.peaks.get(i0).getIntensity(), .5) * Math.pow(s1.peaks.get(i1).getIntensity(), .5);
          Set<List<Integer>> currVal;
          if (pairsPerScore.containsKey(score)) {
              currVal = pairsPerScore.get(score);
          }
          else {
              currVal = new HashSet<List<Integer>>();
          }
          currVal.add(l);
          pairsPerScore.put(score, currVal);
      }
  }
  
  public void updatePairsPerScoreMap2(
          HashMap<Double, Set<List<Integer>>> pairsPerScore, Set<List<Integer>> alignments, Spectrum s1
  ) {
      Iterator<List<Integer>> it = alignments.iterator();
      while(it.hasNext()) {
          List<Integer> l = it.next();
          int i0 = l.get(0);
          int i1 = l.get(1);
          double score = Math.pow(this.peaks.get(i0).getIntensity(), 2) * Math.pow(s1.peaks.get(i1).getIntensity(), 2);
          Set<List<Integer>> currVal;
          if (pairsPerScore.containsKey(score)) {
              currVal = pairsPerScore.get(score);
          }
          else {
              currVal = new HashSet<List<Integer>>();
          }
          currVal.add(l);
          pairsPerScore.put(score, currVal);
      }
  }
 
  public void updatePairsPerScoreMap3(
          HashMap<Double, Set<List<Integer>>> pairsPerScore, Set<List<Integer>> alignments, Spectrum s1
  ) {
      Iterator<List<Integer>> it = alignments.iterator();
      while(it.hasNext()) {
          List<Integer> l = it.next();
          int i0 = l.get(0);
          int i1 = l.get(1);
          double score = Math.pow(this.peaks.get(i0).getIntensity() * s1.peaks.get(i1).getIntensity(), .5);
          Set<List<Integer>> currVal;
          if (pairsPerScore.containsKey(score)) {
              currVal = pairsPerScore.get(score);
          }
          else {
              currVal = new HashSet<List<Integer>>();
          }
          currVal.add(l);
          pairsPerScore.put(score, currVal);
      }
  }
  
  public Set<List<Integer>> findMatchPeaks(Spectrum s1, double shift, double tolerance) {
     double adjTolerance =  tolerance + 0.000001;
     int low = 0;
     int high = 0;
     Set<List<Integer>> matchPeaks = new HashSet<List<Integer>>();
     for (int i = 0; i < this.peaks.size(); i++) {
         Peak peak = this.peaks.get(i);
         low = s1.peaks.size() - 1;
         while (low > 0 && (peak.getMass() - adjTolerance) < (s1.peaks.get(low).getMass() + shift)) {
            low--;
         }
         while (low < s1.peaks.size() && (peak.getMass() - adjTolerance) > (s1.peaks.get(low).getMass() + shift)){
             low++;
         }
         while (high < s1.peaks.size() && (peak.getMass() + adjTolerance) >= (s1.peaks.get(high).getMass() + shift)) {
             high++;
         }
         //return list of pairs of index of "this" peak, index of s1 peak
         for (int j = low; j < high; j++) {
             List<Integer> pair = new ArrayList<Integer>(2);
             pair.add(i);
             pair.add(j);
             matchPeaks.add(pair);
         }
     }
     return matchPeaks;
  }
  
}
