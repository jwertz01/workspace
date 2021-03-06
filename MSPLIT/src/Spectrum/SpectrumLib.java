package Spectrum;

//a spectrum library that represents a collection of spectrum
//index by their corresponding peptides
//thus each peptide is a key that map to a list of
//spectrums that corresponds to the peptide

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStreamWriter;
import java.util.Hashtable;
import java.util.ArrayList;
import java.util.Vector;
import java.util.Collection;
import java.util.Iterator;
import java.util.Collections;
import java.util.TreeMap;
import java.io.Serializable;
import java.util.Random;
import java.util.GregorianCalendar;
//import org.jgrapht.*;
//import org.jgrapht.graph.*;
//import org.jgrapht.alg.*;
import java.util.Set;
import java.util.List;
import java.util.Map;

import IO.LargeSpectrumLibIterator;
import IO.SortedMZXMLReader;
import IO.SortedSpectrumReader;

/**
 * This class represent a spectral library, it contains two views on the spectral library, 
 * a Map view that index the spectra by its peptide and a List view that contains all the spectra
 * 
 * @author Jian
 *
 */
public class SpectrumLib implements Iterable, Serializable{
    public static final long serialVersionUID = 1L; //for serializable interface, see doc
    public static boolean DETAIL = true;
    public static boolean NODETAIL = false;
    private static boolean DEBUG = false;
    private double parentMassTolerance = 2000;
    private double fragmentTolerance = 0.5;
    private int topSpectraKept = 500;
    private int maxPairs=15000;
    private Map<String, List<Spectrum>> spectrumLibrary; 
    private String Filename = "";
    private String QueryFilename ="";
    private String outputFile ="";
    private BufferedWriter out=null;
    
    private List<Spectrum> spectrumList;
    public SpectrumLib(){
        this.spectrumLibrary = new Hashtable();
        this.spectrumList = null;
        this.out = new BufferedWriter(new OutputStreamWriter(System.out)); //set default out to stdout
    }
        
    public SpectrumLib(Map<String, List<Spectrum>> lib){
        this.spectrumLibrary = lib;
        this.spectrumList = this.getAllSpectrums();
    }
    
    public SpectrumLib(String file){
        this(file, file.substring(file.lastIndexOf(".")+1));        
    }
    
    public SpectrumLib(String file, String format){
        this();
        System.out.println("library format: " + format);
        format = format.toUpperCase();
        this.Filename = file;
        if(format.equals("MGF")){
            this.readSpectrumsFromMGF(file);
        }else if(format.equals("MSP")){
            this.readSpectrumsFromMSP(file);
        }else if(format.equals("SPTXT")){
            this.readSpectrumsFromSplib(file);
        }else if(format.equals("MS2")){
            this.readSpectrumsFromMS2(file);
        }else{
            throw new IllegalArgumentException("Unrecognizable library format");
        }
        this.spectrumList = this.getAllSpectrums();
        System.out.println("Read in spectra: " + this.spectrumList.size());
    }
    
    public void initOutPut(String outfile){
        try{
            this.outputFile = outfile;
            this.out = new BufferedWriter(new FileWriter(this.outputFile));
            
        }catch(IOException ioe){
            System.err.println("Error creating output file");
            System.err.println(ioe.getMessage());
        }
    }
    //read a list of spectrums from MGF file format
    public SpectrumLib readSpectrumsFromMGF(String fileName){
        try{
            BufferedReader bf = new BufferedReader(new FileReader(fileName));
            Spectrum s = new Spectrum();
            boolean success = s.readSpectrumFromMGF(bf);
            Vector v;
            int count = 1;
            while(success){
                //System.out.println(s.spectrumName + "\t" + s.peptide);
                if(spectrumLibrary.containsKey(s.peptide)){
                    v = (Vector)spectrumLibrary.get(s.peptide);
                }else{
                    v = new Vector();
                }
                s.specIndex = count++;
                v.add(s);
                spectrumLibrary.put(s.peptide, v);
                s = new Spectrum();
                success = s.readSpectrumFromMGF(bf);
                //System.out.println(count);
                if(count % 10000 == 0){
                    System.out.println("read in spectra: " +  count);
                }
            }        
        }catch(IOException ioe){
            System.out.println("Cannot Open MGF file");
            System.out.println(ioe.getMessage());
        }
        System.out.println("Done reading spectrum File");
        return this;
    }
    
    public SpectrumLib readDecoySpectrumsFromMGF(String fileName){
        try{
            BufferedReader bf = new BufferedReader(new FileReader(fileName));
            Spectrum s = new Spectrum();
            boolean success = s.readSpectrumFromMGF(bf);
            Vector v;
            while(success){
                s.shiftSpectrum(20);
                s.peptide = "X"+s.peptide;
                if(spectrumLibrary.containsKey(s.peptide)){
                    v = (Vector)spectrumLibrary.get(s.peptide);
                }else{
                    v = new Vector();
                }
                v.add(s);
                spectrumLibrary.put(s.peptide, v);
                s = new Spectrum();
                success = s.readSpectrumFromMGF(bf);
            }        
            
        }catch(IOException ioe){
            System.out.println("Cannot Open MGF file");
            System.out.println(ioe.getMessage());
        }
        return this;
    }
    
    public SpectrumLib readSpectrumsFromSplib(String fileName){
        try{
            BufferedReader bf = new BufferedReader(new FileReader(fileName));
            Spectrum s = new Spectrum();
            boolean success = s.readSpectrumFromSplib(bf);
            Vector v;
            int count = 1;
            while(success){
                if(spectrumLibrary.containsKey(s.peptide)){
                    v = (Vector)spectrumLibrary.get(s.peptide);
                }else{
                    v = new Vector();
                }
                s.specIndex=count++;
                v.add(s);
                spectrumLibrary.put(s.peptide, v);
                s = new Spectrum();
                success = s.readSpectrumFromSplib(bf);
            }        
            
        }catch(IOException ioe){
            System.out.println("Cannot Open splib file");
            System.out.println(ioe.getMessage());
        }
        return this;
    }
    
    
    public SpectrumLib readSpectrumsFromMS2(String fileName){
        try{
            BufferedReader bf = new BufferedReader(new FileReader(fileName));
            Spectrum s = new Spectrum();
            boolean success = s.readSpectrumFromMS2(bf);
            Vector v;
            int count=1;
            while(success){
                if(spectrumLibrary.containsKey(s.peptide)){
                    v = (Vector)spectrumLibrary.get(s.peptide);
                }else{
                    v = new Vector();
                }
                s.specIndex = count++;
                v.add(s);
                spectrumLibrary.put(s.peptide, v);
                s = new Spectrum();
                success = s.readSpectrumFromMS2(bf);
            }        
            
        }catch(IOException ioe){
            System.out.println("Cannot Open MS2 file");
            System.out.println(ioe.getMessage());
            System.out.println(ioe.getStackTrace());
        }
        return this;
    }
    
    public SpectrumLib readSpectrumsFromMSP(String fileName){
        try{
            BufferedReader bf = new BufferedReader(new FileReader(fileName));
            Spectrum s = new Spectrum();;
            boolean success = s.readSpectrumFromMSP(bf);
            //System.out.println("reading is: " + success);
            Vector v;
            int count=1;
            while(success){
                //System.out.println(s) ;
                if(spectrumLibrary.containsKey(s.peptide)){
                    v = (Vector)spectrumLibrary.get(s.peptide);
                }else{
                    v = new Vector();
                }
                s.specIndex=count++;
                v.add(s);
                spectrumLibrary.put(s.peptide, v);
                s = new Spectrum();
                success = s.readSpectrumFromMSP(bf);
                //System.out.println("reading is: " + success);
            }        
            
        //  System.out.println("***********************************************************") ;
        }catch(IOException ioe){
            System.err.println("Cannot Open MSP file");
            System.err.println(ioe.getMessage());
        }
        return this;
    }
    
    /**
     * Provide a list view of the spectral library
     * @return
     */
    public List<Spectrum> getAllSpectrums() {
        //being a little careless here, no strong typing
        Vector  spects = new Vector<Spectrum>();
        Iterator it = this.spectrumLibrary.values().iterator();
        while(it.hasNext()){
            spects.addAll((List<Spectrum>)it.next());
        }
        return spects;
 
    }
    
    /**
     * Add a spectrum into spectral library
     * @param s
     */
    public void addSpectrum(Spectrum s){
        Vector v;
        if(spectrumLibrary.containsKey(s.peptide)){
            v = (Vector)spectrumLibrary.get(s.peptide);
        }else{
            v = new Vector();
        }
        v.add(s);
        spectrumLibrary.put(s.peptide, v);
        this.spectrumList.add(s);
    }
    
    public Iterator<Spectrum> iterator(){
        if(this.spectrumList == null){
            this.spectrumList = this.getAllSpectrums();
        }
        return this.spectrumList.iterator();
    }
    
    //the following few methods basically just apply some method to every
    //spectrum in the library, to be more elegant we should have implement
    //this as some higher order function, but for simplicity we will just 
    //implement this one by one
    
    //transform every spectrum in this libary to a normalized vector
    public void toNormVector(double binWidth, double minMass, double maxMass){
        Iterator<Spectrum> it = this.iterator();
        Spectrum s, v;
        while(it.hasNext()){
            s = (Spectrum)it.next();
            v = s.toNormVector(binWidth, minMass, maxMass);
            s.setPeaks(v.getPeak());
        }
    }
    
    //transform every spectrum in this libary to a vectorized one
    public void toVector(double binWidth, double minMass, double maxMass){
        Iterator<Spectrum> it = this.iterator();
        Spectrum s, v;
        while(it.hasNext()){
            s = (Spectrum)it.next();
            v = s.toVector(binWidth, minMass, maxMass);
            s.setPeaks(v.getPeak());
        }
    }
    
    //normalize the intensity by taking its square root
    public void normIntensity(){
        Iterator<Spectrum> it = this.iterator();
        while(it.hasNext()){
            it.next().sqrtSpectrum();
        }
    }
    
    //shift all the spectrum by a mass of delta
    public void shiftSpectrum(double delta){
        Iterator<Spectrum> it = this.iterator();
        while(it.hasNext()){
            it.next().shiftSpectrum(delta);
        }
    }
    
    public void scaleSpectrumMass(double scale){
        Iterator<Spectrum> it = this.iterator();
        while(it.hasNext()){
            it.next().scaleMass(scale);
        }
    }

    public String toString(){
        StringBuffer sb = new StringBuffer();    
        Iterator it = this.iterator();
        Spectrum spect;
        while(it.hasNext()){
            spect = (Spectrum)it.next();
            sb.append(spect.toString());
            sb.append("\n");
        }
        return sb.toString();   
    }
    
    //printing varous info for debuging
    public static void printDetail(String s){
        if(SpectrumLib.DEBUG){
            System.out.println(s);
        }
    }
    
    /**
     * Print out statistics about the spectral library
     * @param detail
     */
    public void printStat(boolean detail){
        StringBuffer sb = new StringBuffer();    
        Iterator it = this.spectrumLibrary.keySet().iterator();
        List <Spectrum> spects;
        String pep;
        int singleCount = 0;
        int multipleCount = 0;
        //counting spectrums for each peptide, number
        //of peptide with single and multiple spectrums
        while(it.hasNext()){
            pep = (String)it.next();
            spects = spectrumLibrary.get(pep);
            sb.append("peptide: " + pep 
                    + " # spectrums:\t" +  spects.size() 
                    + "\n");
            if(spects.size() > 1){
                multipleCount++;
            }else{
                singleCount++;
            }       
        }
        
        if(detail){
            System.out.println(sb.toString());
        }
        System.out.println("Total Number of Peptide: " + (singleCount + multipleCount));
        System.out.println("Peptide with single spectrum : " + singleCount);
        System.out.println("Peptide with multiple spctra : " + multipleCount);
    }
    
    /**
     * Remove spectra whose annotation contain modification
     */
    public void removeModSpectra(){
        Collection <List <Spectrum>> values = this.spectrumLibrary.values();
        Iterator <List <Spectrum>> it = values.iterator();
        List <Spectrum>  v;
        int index = 0;
        //iterate over each values, which is a list of spectrums
        while(it.hasNext()){
            v = it.next();
            for(index = 0; index < v.size(); index++){
                if(v.get(index).modMass > 0){
                    v.remove(index);
                    index--; //since we just remove one element, backtrack one pointer
                }
                                
            }
            //if this peptide has no more peptide remove it as well
            if(v.size() == 0){
                it.remove();
            }
        }
        //the "naive" way to update the vector
        this.spectrumList = this.getAllSpectrums();
        
    }
    
    
    public void removeSpectrum(String id){
        List<Spectrum> list = this.spectrumLibrary.get(id);
        this.spectrumLibrary.remove(id);
        if(this.spectrumList != null){
            this.spectrumList.removeAll(list);
        }
    }
    
    /**
     * Remove all spectral library entry with only one spectrum
     */
    public void removeSingle(){
        Collection <List <Spectrum>> values = this.spectrumLibrary.values();
        Iterator <List <Spectrum>> it = values.iterator();
        //iterate over each values, which is a list of spectrums
        List <Spectrum> v;
        int index = 0;
        while(it.hasNext()){
            v = it.next();
            if(v.size() == 1){
                it.remove();
            }
        }
        //the "naive" way to update the vector
        this.spectrumList = this.getAllSpectrums();
    }
    
//  public void removeLowComplexity(){
//      Collection <List <Spectrum>> values = this.spectrumLibrary.values();
//      Iterator <List <Spectrum>> it = values.iterator();
//      //iterate over each values, which is a list of spectrums
//      List <Spectrum> v;
//      int index = 0;
//      while(it.hasNext()){
//          v = it.next();
//          Spectrum s = v.get(0);
//          s.toRelIntensity();
//          if(s.intensePeakCount(0.1) < 10){
//              it.remove();
//          }
//      }
//      //the "naive" way to update the vector
//      this.spectrumList = this.getAllSpectrums();
//  }
    
    //only retain top n peaks in the spectrums in this library
    public void filterPeaks(int top){
        Iterator<Spectrum> it = this.iterator();
        while(it.hasNext()){
            it.next().filterPeaks(top);
        }
    }
    
    public void windowFilterPeaks(int top, double winSize){
        Iterator<Spectrum> it = this.iterator();
        while(it.hasNext()){
            it.next().windowFilterPeaks(top, winSize);
        }
    }
    
    public void toRelativeIntensity(){
        Iterator<Spectrum> it = this.iterator();
        while(it.hasNext()){
            it.next().toRelIntensity();
        }
    }
    
    public void computeRank(){
        Iterator<Spectrum> it = this.iterator();
        while(it.hasNext()){
            it.next().computePeakRank();
        }
    }
    
    /**
     * We divide the dataset into two sets 
     * one to simulate mix spectrum, the other
     * to perform the search, note this method modify the
     * original spectLib so those spectra that are separted
     * into one set is removed from this Spectrumlib
     * @return
     */
    public SpectrumLib Divide(){
        Hashtable <String, List <Spectrum>> newTable = new Hashtable();
        Collection <List <Spectrum>> values = this.spectrumLibrary.values();
        Iterator <List <Spectrum>> it = values.iterator();
        //iterate over each values, which is a list of spectrums
        List <Spectrum> v, vNew;
        while(it.hasNext()){
            v = it.next();
            //if the peptide has more than one peptide, we take the first one out
            //to create a different spectrum library, otherwise we don't use it
            //cause there is no way we can identitfy mix spectrum that consist of that peptide
            if(v.size() > 1){
                vNew = new Vector <Spectrum>();
                vNew.add(v.get(0));
                newTable.put(v.get(0).peptide, vNew);
                v.remove(0);
            }
        }
        this.spectrumList = this.getAllSpectrums();
        return new SpectrumLib(newTable);
    }
    
    /**
     * Creating a new spectrum library where we mix a pair
     * of spectrum in this library to simulated a mixture spectrum
     * @param size
     * @return
     */
    public SpectrumLib createMix(int size, double fragmentTolerance){
        return createMix(size, 1, 1, fragmentTolerance);
    }
    
    //edited: for now, if we see spectrum where there is one very high peak
    //and all very low peaks we ignore it for now to see performance of cosineSim function
    public SpectrumLib createMix(int size, double scale, double scale2, double tolerance){
        Hashtable <String, List <Spectrum>> newTable = new Hashtable();
        Vector <List <Spectrum>> v = new Vector(this.spectrumLibrary.values());
        List <Spectrum> vNew; 
        Spectrum mixture;
        for(int i = 0; i < size; i++){
            for(int j = i+1; j < size; j++){
                    mixture = new Spectrum(v.get(i).get(0), v.get(j).get(0), scale, scale2, fragmentTolerance);
                    vNew = (new Vector());
                    vNew.add(mixture);
                    newTable.put(mixture.peptide, vNew);
                
            }
        }
        return new SpectrumLib(newTable);
    }
    
    //create mix with cosineSim filter
    public SpectrumLib createMix(int size, double scale, double scale2, double minSim, double maxSim, double fragmentTolerance){
        Hashtable <String, List <Spectrum>> newTable = new Hashtable();
        List <List <Spectrum>> v = new Vector(this.spectrumLibrary.values());
        List <Spectrum> vNew; 
        Spectrum mixture, s1, s2;
        int counts = 0;
        for(int i = 0; i < this.spectrumLibrary.keySet().size(); i++){
            for(int j = i+1; j < this.spectrumLibrary.keySet().size(); j++){
                s1 = v.get(i).get(0);
                s2 = v.get(j).get(0);
                if((s1.cosineSim(s2, fragmentTolerance) >= minSim) && (s1.cosineSim(s2, fragmentTolerance) < maxSim)){
                    mixture = new Spectrum(v.get(i).get(0), v.get(j).get(0), scale, scale2, fragmentTolerance);
                    mixture.score = s1.cosineSim(s2, fragmentTolerance);
                    vNew = (new Vector());
                    vNew.add(mixture);
                    newTable.put(mixture.peptide, vNew);
                    counts++;
                }
                if(counts > size){
                    return new SpectrumLib(newTable);
                }
                
            }
        }
        return new SpectrumLib(newTable);
    }
    
    //with mass difference filters
    public SpectrumLib createMix(int size, double scale, double scale2, double minSim, double maxSim, double massDiff, double fragmentTolerance){
        Hashtable <String, List <Spectrum>> newTable = new Hashtable();
        List <List <Spectrum>> v = new Vector(this.spectrumLibrary.values());
        List <Spectrum> vNew; 
        Spectrum mixture, s1, s2;
        int counts = 0;
        for(int i = 0; i < this.spectrumLibrary.keySet().size(); i++){
            for(int j = i+1; j < this.spectrumLibrary.keySet().size(); j++){
                s1 = v.get(i).get(0);
                s2 = v.get(j).get(0);
                if((s1.cosineSim(s2, fragmentTolerance) >= minSim) && (s1.cosineSim(s2, fragmentTolerance) < maxSim)  //also require same charge for now
                        && (Math.abs((s1.parentMass - s2.parentMass))  < massDiff)){
                    mixture = new Spectrum(v.get(i).get(0), v.get(j).get(0), scale, scale2, fragmentTolerance);
                    mixture.score = s1.cosineSim(s2, fragmentTolerance);
                    vNew = (new Vector());
                    vNew.add(mixture);
                    newTable.put(mixture.peptide, vNew);
                    counts++;
                }
                if(counts > size){
                    return new SpectrumLib(newTable);
                }
            }
        }
        return new SpectrumLib(newTable);
    }
    
    //when lib itself is very large,  iterating through each spectra to create mix is not a good way to 
    //sample the space of all possible mix so we choose to use a randomize scheme, doing it this way may lead
    //to duplicate mix spectra, but the chance is small if lib is big, so it is okay, since we are hashing using peptide as the key
    public SpectrumLib createRandomMix(int size, double scale, double scale2, double minSim, double maxSim, double massDiff, double fragmentTolerance){
        Hashtable <String, List <Spectrum>> newTable = new Hashtable();
        Vector <Vector <Spectrum>> v = new Vector(this.spectrumLibrary.values());
        Vector <Spectrum> vNew; 
        Spectrum mixture, s1, s2;
        int counts = 0, i = 0, j = 0;
        while(counts < size){
                i = (int) (Math.random()*this.spectrumLibrary.keySet().size());
                j =  (int)(Math.random()*this.spectrumLibrary.keySet().size());
                s1 = v.get(i).get(0);
                s2 = v.get(j).get(0);
                if((s1.cosineSim1(s2, fragmentTolerance) >= minSim) && (s1.cosineSim1(s2, fragmentTolerance) < maxSim)  //also require same charge for now
                        && (Math.abs((s1.parentMass - s2.parentMass))  < massDiff)
                        &&  i != j){  //we do not mix two same spectrum
                    mixture = new Spectrum(v.get(i).get(0), v.get(j).get(0), scale, scale2, fragmentTolerance);
                    mixture.score = s1.cosineSim1(s2, fragmentTolerance);
                    //System.out.println(mixture.peptide + "\t" + mixture.score);
                    vNew = (new Vector());
                    vNew.add(mixture);
                    newTable.put(mixture.peptide, vNew);
                    counts++;
                }       
        }   
        return new SpectrumLib(newTable);
    }
    
    public SpectrumLib createRandomMix(SpectrumLib lib, int size, double scale, double scale2, double minSim, double maxSim, double massDiff, double fragmentTolerance){
        Hashtable <String, List <Spectrum>> newTable = new Hashtable();
        Vector <Vector <Spectrum>> v = new Vector(this.spectrumLibrary.values());
        Vector <Vector <Spectrum>> w = new Vector(lib.spectrumLibrary.values());
        Vector <Spectrum> vNew; 
        Spectrum mixture, s1, s2;
        int counts = 0, i = 0, j = 0;
        while(counts < size){
                i = (int) (Math.random()*this.spectrumLibrary.keySet().size());
                j =  (int)(Math.random()*lib.spectrumLibrary.keySet().size());
                s1 = v.get(i).get(0);
                s2 = w.get(j).get(0);
                if((s1.cosineSim1(s2, fragmentTolerance) >= minSim) && (s1.cosineSim1(s2, fragmentTolerance) < maxSim)  //also require same charge for now
                        && (Math.abs((s1.parentMass - s2.parentMass))  < massDiff)
                        &&  i != j){  //we do not mix two same spectrum
                    mixture = new Spectrum(v.get(i).get(0), w.get(j).get(0), scale, scale2, fragmentTolerance);
                    mixture.score = s1.cosineSim1(s2, fragmentTolerance);
                    //System.out.println(mixture.peptide + "\t" + mixture.score);
                    vNew = (new Vector());
                    vNew.add(mixture);
                    newTable.put(mixture.peptide, vNew);
                    counts++;
                }       
        }   
        return new SpectrumLib(newTable);
    }
    
    public SpectrumLib createRandomMix(SpectrumLib lib, int size, double minScale, double minSim, double maxSim, double massDiff, boolean toVector, double fragmentTolerance){
        Hashtable <String, List <Spectrum>> newTable = new Hashtable();
        Vector <Vector <Spectrum>> v = new Vector(this.spectrumLibrary.values());
        Vector <Vector <Spectrum>> w = new Vector(lib.spectrumLibrary.values());
        Vector <Spectrum> vNew; 
        Spectrum mixture, s1, s2;
        int counts = 0, i = 0, j = 0;
        while(counts < size){
                i = (int) (Math.random()*this.spectrumLibrary.keySet().size());
                j =  (int)(Math.random()*lib.spectrumLibrary.keySet().size());
                s1 = v.get(i).get(0);
                s2 = w.get(j).get(0);
                double scale = Math.random();
                scale = scale*(1.0-minScale)+minScale;
                if((s1.cosineSim1(s2, fragmentTolerance) >= minSim) && (s1.cosineSim1(s2, fragmentTolerance) < maxSim)  //also require same charge for now
                        && (Math.abs((s1.parentMass - s2.parentMass))  < massDiff)
                        &&  i != j){  //we do not mix two same spectrum
                    mixture = new Spectrum(v.get(i).get(0), w.get(j).get(0), 1.0, scale, toVector, fragmentTolerance);
                    mixture.score = s1.cosineSim1(s2, fragmentTolerance);
                    //System.out.println(mixture.peptide + "\t" + mixture.score);
                    vNew = (new Vector());
                    vNew.add(mixture);
                    newTable.put(mixture.peptide, vNew);
                    counts++;
                }       
        }   
        return new SpectrumLib(newTable);
    }
    
    public SpectrumLib createRandomMix(int size, double scale, double scale2, double minSim, double maxSim, double massDiff, boolean toVector, double fragmentTolerance){
        Hashtable <String, List <Spectrum>> newTable = new Hashtable();
        Vector <Vector <Spectrum>> v = new Vector(this.spectrumLibrary.values());
        Vector <Spectrum> vNew; 
        Spectrum mixture, s1, s2;
        int counts = 0, i = 0, j = 0;
        while(counts < size){
                i = (int) (Math.random()*this.spectrumLibrary.keySet().size());
                j =  (int)(Math.random()*this.spectrumLibrary.keySet().size());
                s1 = v.get(i).get(0);
                s2 = v.get(j).get(0);
                if((s1.cosineSim1(s2, fragmentTolerance) >= minSim) && (s1.cosineSim1(s2, fragmentTolerance) < maxSim)  //also require same charge for now
                        && (Math.abs((s1.parentMass - s2.parentMass))  < massDiff)
                        &&  i != j
                        && s1.charge <= 3
                        && s2.charge <= 3){  //we do not mix two same spectrum
                    mixture = new Spectrum(v.get(i).get(0), v.get(j).get(0), scale, scale2, toVector, fragmentTolerance);
                    mixture.score = s1.cosineSim1(s2, fragmentTolerance);
                    //System.out.println(mixture.peptide + "\t" + mixture.score);
                    vNew = (new Vector());
                    vNew.add(mixture);
                    newTable.put(mixture.peptide, vNew);
                    counts++;
                }       
        }   
        return new SpectrumLib(newTable);
    }

    public SpectrumLib createMix(int size, double scale, double scale2, double minSim, double maxSim, double massDiff, boolean toVector, double fragmentTolerance){
        Hashtable <String, List <Spectrum>> newTable = new Hashtable();
        List <List <Spectrum>> v = new Vector(this.spectrumLibrary.values());
        List <Spectrum> vNew; 
        Spectrum mixture, s1, s2;
        int counts = 0;
        for(int i = 0; i < this.spectrumLibrary.keySet().size(); i++){
            for(int j = i+1; j < this.spectrumLibrary.keySet().size(); j++){
                s1 = v.get(i).get(0);
                s2 = v.get(j).get(0);
                if((s1.cosineSim(s2, fragmentTolerance) >= minSim) && (s1.cosineSim(s2, fragmentTolerance) < maxSim)  //also require same charge for now
                        && (Math.abs((s1.parentMass - s2.parentMass))  < massDiff)){
                    mixture = new Spectrum(v.get(i).get(0), v.get(j).get(0), scale, scale2, toVector, fragmentTolerance);
                    mixture.score = s1.cosineSim(s2, fragmentTolerance);
                    vNew = (new Vector());
                    vNew.add(mixture);
                    newTable.put(mixture.peptide, vNew);
                    counts++;
                }
                if(counts > size){
                    return new SpectrumLib(newTable);
                }
                
            }
        }
        return new SpectrumLib(newTable);
    }
    
    public SpectrumLib createMix(int size, double minScale, double minSim, double maxSim, double massDiff, boolean toVector, double fragmentTolerance){
        Hashtable <String, List <Spectrum>> newTable = new Hashtable();
        List <List <Spectrum>> v = new Vector(this.spectrumLibrary.values());
        List <Spectrum> vNew; 
        Spectrum mixture, s1, s2;
        int counts = 0;
        for(int i = 0; i < this.spectrumLibrary.keySet().size(); i++){
            for(int j = i+1; j < this.spectrumLibrary.keySet().size(); j++){
                s1 = v.get(i).get(0);
                s2 = v.get(j).get(0);
                if((s1.cosineSim(s2, fragmentTolerance) >= minSim) && (s1.cosineSim(s2, fragmentTolerance) < maxSim)  //also require same charge for now
                        && (Math.abs((s1.parentMass - s2.parentMass))  < massDiff)
                        && (s1.charge >=2 && s1.charge <=3 ) 
                        && (s2.charge >= 2 && s2.charge<=3)){
                    double scale = Math.random();
                    scale = scale*(1.0-minScale)+minScale;
                    mixture = new Spectrum(v.get(i).get(0), v.get(j).get(0), 1.0, scale, toVector, fragmentTolerance);
                    mixture.score = s1.cosineSim(s2, fragmentTolerance);
                    vNew = (new Vector());
                    vNew.add(mixture);
                    newTable.put(mixture.peptide, vNew);
                    counts++;
                }
                if(counts > size){
                    return new SpectrumLib(newTable);
                }
                
            }
        }
        return new SpectrumLib(newTable);
    }
    
    public SpectrumLib createRandomMix(int size, double minScale, double minSim, double maxSim, double massDiff, boolean toVector, double fragmentTolerance){
        Hashtable <String, List <Spectrum>> newTable = new Hashtable();
        Vector <Vector <Spectrum>> v = new Vector(this.spectrumLibrary.values());
        Vector <Spectrum> vNew; 
        Spectrum mixture, s1, s2;
        int counts = 0, i = 0, j = 0;
        while(counts < size){
                i = (int) (Math.random()*this.spectrumLibrary.keySet().size());
                j =  (int)(Math.random()*this.spectrumLibrary.keySet().size());
                s1 = v.get(i).get(0);
                s2 = v.get(j).get(0);
                if((s1.cosineSim1(s2, fragmentTolerance) >= minSim) && (s1.cosineSim1(s2, fragmentTolerance) < maxSim)  //also require same charge for now
                        && (Math.abs((s1.parentMass - s2.parentMass))  < massDiff)
                        &&  i != j
                        && (s1.charge >= 2 && s1.charge <= 3)
                        && (s2.charge >= 2 && s2.charge <=3)){  //we do not mix two same spectrum
                    double scale = Math.random();
                    scale = scale*(1.0-minScale)+minScale;
                    System.out.println("creating mix " + s1.getPeptide() + " & " +  s2.getPeptide() + " with alpha: " + scale); 
                    mixture = new Spectrum(v.get(i).get(0), v.get(j).get(0), 1.0, scale, toVector, fragmentTolerance);
                    mixture.score = s1.cosineSim1(s2, fragmentTolerance);
                    //System.out.println(mixture.peptide + "\t" + mixture.score);
                    vNew = (new Vector());
                    vNew.add(mixture);
                    newTable.put(mixture.peptide, vNew);
                    counts++;
                }       
        }   
        return new SpectrumLib(newTable);
    }
    
    //create mixture from specified peptide pair in a file
    public SpectrumLib createMix(String file, int size, double scale1, double scale2, double minSim, double maxSim, double massDiff, double fragmentTolerance){
        Spectrum s1, s2, mixture;
        SpectrumLib lib = null;
        Hashtable <String, List <Spectrum>> newTable = new Hashtable();
        int counts = 0;
        List <Spectrum> vNew; 
        try{
            BufferedReader bf = new BufferedReader(new FileReader(file));
            String[] peps = null;
            String line = bf.readLine();
            while(line != null){
                peps = line.split(" & ");
                s1 = this.spectrumLibrary.get(peps[0]).get(0);
                s2 = this.spectrumLibrary.get(peps[1]).get(0);
                //s1.filterPeaks(75);
                //s2.filterPeaks(25);
                mixture = new Spectrum(s1, s2, scale1, scale2, fragmentTolerance);
                mixture.score = s1.cosineSim(s2, fragmentTolerance);
                if((mixture.score >= minSim) && (mixture.score < maxSim)  //also require same charge for now
                        && (Math.abs((s1.parentMass - s2.parentMass))  < massDiff)){  
                    vNew = new Vector<Spectrum>();
                    vNew.add(mixture);
                    newTable.put(mixture.peptide, vNew);
                    //System.out.println("" + counts);
                    //System.out.println(mixture.peptide + "\t" + mixture.score);
                    counts++;
                    if(counts > size){
                        return new SpectrumLib(newTable);
                    }
                }
                line = bf.readLine();
            }
            return new SpectrumLib(newTable);
        }catch(IOException ioe){
            System.out.println(ioe.getMessage());
        }catch(NullPointerException e){
            System.out.println("specified peptide is not found in the library");
            System.out.println(e.getMessage());
            System.out.println(e.getStackTrace());
        }
        return lib;
    }
    
    public SpectrumLib createMix(String file, int size, double scale1, double scale2, double minSim, double maxSim, double massDiff, boolean toVector, double fragmentTolerance){
        Spectrum s1, s2, mixture;
        SpectrumLib lib = null;
        Hashtable <String, List <Spectrum>> newTable = new Hashtable();
        int counts = 0;
        List <Spectrum> vNew; 
        try{
            BufferedReader bf = new BufferedReader(new FileReader(file));
            String[] peps = null;
            String line = bf.readLine();
            while(line != null){
                peps = line.split(" & ");
                s1 = this.spectrumLibrary.get(peps[0]).get(0);
                s2 = this.spectrumLibrary.get(peps[1]).get(0);
                //s1.filterPeaks(75);
                //s2.filterPeaks(25);
                mixture = new Spectrum(s1, s2, scale1, scale2, toVector, fragmentTolerance);
                mixture.score = s1.cosineSim(s2, fragmentTolerance);
                if((mixture.score >= minSim) && (mixture.score < maxSim)  //also require same charge for now
                        && (Math.abs((s1.parentMass - s2.parentMass))  < massDiff)){  
                    vNew = new Vector<Spectrum>();
                    vNew.add(mixture);
                    newTable.put(mixture.peptide, vNew);
                    //System.out.println("" + counts);
                    //System.out.println(mixture.peptide + "\t" + mixture.score);
                    counts++;
                    if(counts > size){
                        return new SpectrumLib(newTable);
                    }
                }
                line = bf.readLine();
            }
            return new SpectrumLib(newTable);
        }catch(IOException ioe){
            System.out.println(ioe.getMessage());
        }catch(NullPointerException e){
            System.out.println("specified peptide is not found in the library");
            System.out.println(e.getMessage());
            System.out.println(e.getStackTrace());
        }
        return lib;
    }
    
    public SpectrumLib createMix(String file, int size, double minSim, double maxSim, double massDiff, boolean toVector, double fragmentTolerance){
        Spectrum s1, s2, mixture;
        SpectrumLib lib = null;
        Hashtable <String, List <Spectrum>> newTable = new Hashtable();
        int counts = 0;
        List <Spectrum> vNew;
        double alpha = 0;
        try{
            BufferedReader bf = new BufferedReader(new FileReader(file));
            String[] tokens = null;
            String line = bf.readLine();
            while(line != null){
                tokens = line.split("\t");
                if(this.spectrumLibrary.containsKey(tokens[0]) && this.spectrumLibrary.containsKey(tokens[1])){
                    s1 = this.spectrumLibrary.get(tokens[0]).get(0);
                    s2 = this.spectrumLibrary.get(tokens[1]).get(0);
                }else{
                    line = bf.readLine();
                    continue;
                }
                alpha = Double.parseDouble(tokens[2]);
                mixture = new Spectrum(s1, s2, 1.0, alpha, toVector, false, fragmentTolerance);     
                mixture.score = s1.cosineSim(s2, fragmentTolerance);
                if((mixture.score >= minSim) && (mixture.score < maxSim)  //also require same charge for now
                        && (Math.abs((s1.parentMass - s2.parentMass))  < massDiff)){  
                    vNew = new Vector<Spectrum>();
                    vNew.add(mixture);
                    newTable.put(mixture.peptide, vNew);
                    //System.out.println("" + counts);
                    //System.out.println(mixture.peptide + "\t" + mixture.score);
                    counts++;
                    if(counts > size){
                        return new SpectrumLib(newTable);
                    }
                }
                line = bf.readLine();
            }
            return new SpectrumLib(newTable);
        }catch(IOException ioe){
            System.out.println(ioe.getMessage());
        }catch(NullPointerException e){
            System.out.println("specified peptide is not found in the library");
            System.out.println(e.getMessage());
            System.out.println(e.getStackTrace());
        }
        return lib;
    }
    
    
    public List<Spectrum> getSpectra(String peptide){
        return this.spectrumLibrary.get(peptide);
    }
    
    public Spectrum getSpectrumById(String id){
        for(int i = 0; i < this.spectrumList.size(); i++){
            Spectrum curr = this.spectrumList.get(i);
            if(curr.spectrumName.equals(id)){
                return curr;
            }
        }
        throw new IllegalArgumentException("cannot find specified spectrum in library");
    }
    
    /**since multiple spectrum can correspond to same
     * peptide, this method try to see how similar are
     * the spectrum corresponds to same peptide
     **/
    public void getSpectrumSimilarity(int size, double fragmentTolerance){
        Iterator<List<Spectrum>> it = this.spectrumLibrary.values().iterator();
        Spectrum first, next;
        List<Spectrum> spects;
        double maxSim = 0;
        int i = 0, j = 0;
        int count = 0;
        while(it.hasNext()){
            spects = it.next();
            //only can calcuate similarity if there is more than one peptide
            if(spects.size() > 1){
                maxSim = 0;
                for(i = 0; i < spects.size(); i++){
                    for(j = i+1; j < spects.size(); j++){
                        first = spects.get(i).toNormVector(1, 0.5, 2000);
                        next = spects.get(j).toNormVector(1, 0.5, 2000);
                        if(maxSim < next.cosineSim(first, fragmentTolerance)){
                            maxSim = next.cosineSim(first, fragmentTolerance);
                        }
                        //System.out.println(first);
                        //System.out.println(next);
                        System.out.println("" + first.peptide + "\t" + next.cosineSim(first, fragmentTolerance));
                        count++;
                        if(count > size){
                            return;
                        }
                    
                    }
                }
                
            }
        }
    }   
    
    
    //we compute the inter-group spectra similarity, it serves
    //as a negative control for measuring similarity for those
    //spectra that belong to the same peptide with same charge
    public void getSpectrumDifference(int size, double fragmentTolerance){
        int count = 0;
        String[] peptides = this.spectrumLibrary.keySet().toArray(new String[10]);
        //because all possible pair in the negative set
        //is usually very large we, set a limit on the total
        //number of scores to generate
        Spectrum s1, s2;
        //for each peptide we randomly select one spectrum to compare
        for(int i = 0; i < peptides.length; i++){
            for(int j = i+1; j < peptides.length; j++){
                s1 = this.getRandomSpectrum(peptides[i]).toNormVector(1, 0.5, 2000);
                s2 = this.getRandomSpectrum(peptides[j]).toNormVector(1, 0.5, 2000);
                System.out.println("" + peptides[i] + " | " + peptides[j] + "\t" + s1.cosineSim(s2, fragmentTolerance));
                if(count == size){
                    return;
                }
                count++;
            }
        }
        
    }
    
    //given a candidate spectrum, score it against all other
    //spectrum in this library and return the rank of similarity
    //of the target peptide
    public int psimilarityRank(Spectrum m, String targetPeptide, double fragmentTolerance){
        int rank1 = 1, rank2 = 1;
        int prank1 = 1, prank2 = 1; //peptide rank
        int rank = 1, prank = 1; //the rank after we remove some peaks
        int count1 = 0, count2 = 0, count = 0;
        int i = 0;
        int topN = 5;  //consider the top N element as from the "better" component of the mix
        Iterator<List<Spectrum>> it = this.spectrumLibrary.values().iterator();
        Iterator<Spectrum> curr;
        Spectrum bestSpect1, bestSpect2, currSpect, bestSpect;
        double best, best2;
        String[] peps = targetPeptide.split(" & ");
        
        //find the highest psimilarity (p for projectedCosine) of the answer
        bestSpect1 = bestPSim(m, peps[0], fragmentTolerance);
        best = bestSpect1.score;
        if(peps.length >= 2){
            bestSpect2 = bestPSim(m, peps[1], fragmentTolerance);
            best2 = bestSpect2.score;
            System.out.println("best score: " + best + "\t" + best2);
        }else{
            bestSpect2 = bestSpect1; //if this is not a mix only one first rank is meaningful 
            best2 = bestSpect2.score;//so make both of them the same
        }
        System.out.println("best score: " + best + "\t" + best2);
        //find the psimiliarity for everything else in the library
        while(it.hasNext()){
            curr = it.next().iterator();
            count1=0;
            count2=0;
            count = 0;
            while(curr.hasNext()){
                currSpect = curr.next();
                //currSpect.score = currSpect.projectedCosine(m);
                currSpect.score = currSpect.cosineSim(m, fragmentTolerance);
                count1 = currSpect.score > best ?  count1+1 : count1;
                count2 = currSpect.score > best2 ?  count2+1 : count2;
                if(currSpect.score > best && false){
                    System.out.println(currSpect.peptide + "\t" + currSpect.score + ": better than 1st peptide");
                }
                if(currSpect.score > best2 && false){
                    System.out.println(currSpect.peptide + "\t" + currSpect.score + ": better than 2nd peptide");
                }
                count++;
            }
            //System.out.println("total count: " + count);
            //if(count1 > 1 || count2 > 1){
            //  System.out.println("count: " + count1 + "\t" + count2);
            //}
            //Here prank stands for peptide ranks rather than spectrum ranks
            rank1 = rank1 + count1;
            rank2 = rank2 + count2;
            prank1 = count1 > 0 ? prank1+1 : prank1;
            prank2 = count2 > 0 ? prank2+1 : prank2;
            if(count1 > 1 || count2 > 1){
                //System.out.println("r:  " + rank1 + "\t" + rank2);
                //System.out.println("pr: " + prank1 + "\t" + prank2);
            }
        }
        
        Collections.sort(this.spectrumList); //we sort the spectrums by scores
        //remove  the peaks and then rerank the 2nd peptide
        for(i = 1; i <= 1; i++){
            currSpect = this.spectrumList.get(spectrumList.size()-i);
            //System.out.println("removing peaks on: "  + currSpect.peptide);
            m.subtractSharePeaks(currSpect, fragmentTolerance);
        }
        if(best > best2){
            bestSpect = bestPSim(m, bestSpect2.peptide, fragmentTolerance);         
        }else{
            bestSpect = bestPSim(m, bestSpect1.peptide, fragmentTolerance);
        }
        best = bestSpect.score;
        //System.out.println("2nd round best:\t" + best);
        it = this.spectrumLibrary.values().iterator();
        while(it.hasNext()){
            curr = it.next().iterator();
            count = 0;
            while(curr.hasNext()){
                currSpect = curr.next();
                currSpect.score = currSpect.cosineSim(m, fragmentTolerance);
                count = currSpect.score > best ?  count+1 : count;
            }
            rank = rank + count;
            prank = count > 0 ? prank+1 : prank;
        }
        
        System.out.println("Psimilarity Rank: " + targetPeptide +  "\t" + m.score + "\t" + rank1 + "\t" +  rank2 + "\t" + rank
                        + "\t" + prank1 + "\t" +  prank2 + "\t" + prank);
        return Math.max(rank1, rank2);
        
    }
    
    //given a spectrum and a pepeitde, return the  spectrum
    //that belong to the peptide that has best pconsine score
    private Spectrum bestPSim(Spectrum m, String peptide, double fragmentTolerance){
        List<Spectrum> targetSpects = this.getSpectra(peptide);
        Spectrum curr, bestSpect = targetSpects.get(0);
        double score, best = 0;
        //find the psimilarity (p for projectedCosine) of the answer
        for(int i = 0; i <  targetSpects.size(); i++){
            curr = targetSpects.get(i);
//          score = curr.projectedCosine(m);
            score = curr.cosineSim(m, fragmentTolerance);
            bestSpect = score > best ? curr : bestSpect;
            best = score > best ? score : best;
        }
        bestSpect.score = best;
        return bestSpect;
    }
    
    //we score the mixspectrum against a combined spectrum of the
    //two peptide given as argument, this way we try to evaluate
    //how similar is a mixed spectrum and a combined spectrum in the library
    public Spectrum mixSpectrumSimilarity(Spectrum mix, String peptide1, String peptide2, double fragmentTolerance){             //NEW CHANGE 
        double score  = 0, bestScore = 0, alpha = 0;
        List<Spectrum> v1, v2;
        int i = 0, j = 0;
        Spectrum targetMix, best = new Spectrum();
        v1 = this.getSpectra(peptide1);
        v2 = this.getSpectra(peptide2);
        for(i = 0; i < v1.size(); i++){
            for(j = 0; j < v2.size(); j++){
                //System.out.println(mix);
                //System.out.println(targetMix);
                //score = mix.cosineSim(targetMix);
                alpha = mix.alpha(v1.get(i), v2.get(j), fragmentTolerance);
                targetMix = mix.toNormVector();
                //targetMix.upscale(v1.get(i), v2.get(j), alpha);
                score = mix.maxScore(v1.get(i), v2.get(j), alpha, fragmentTolerance);
                //score = mix.inverseMaxScore(v1.get(i), v2.get(j), alpha);
                if(score > bestScore){
                    bestScore = score;
                    best = targetMix;
                    best.score = score;
                }
                //System.out.println("" + v1.get(i));
                //System.out.println("" + v2.get(j));
                //System.out.println("" + peptide1 + ": " + v1.get(i).cosineSim(targetMix) + "\t" + mix.shareSim(v1.get(i),1));
                //System.out.println("" + peptide2 + ": " + v2.get(j).cosineSim(targetMix) + "\t" + mix.shareSim(v2.get(j),1));
                System.out.println("" + peptide1 + ": " + v1.get(i).cosineSim(mix, fragmentTolerance) + "\t" + mix.shareSim(v1.get(i),1));
                System.out.println("" + peptide2 + ": " + v2.get(j).cosineSim(mix, fragmentTolerance) + "\t" + mix.shareSim(v2.get(j),1));
                System.out.println("" + targetMix.peptide + "\t" + mix.score + "\t" + score); 
            }
        }
        //System.out.println("" + mix.peptide + " bestscore: " + bestScore);
        return best;                                                                            //NEW CHANGE END 
    } 
    
    //we get the distribution of the negative set
    public void mixSpectrumDifference(Spectrum mix, String peptide1, String peptide2, int size, double fragmentTolerance){             //NEW CHANGE 
        Spectrum curr = null;
        int i = 0;
        for(i = 0; i < size; i++){
            curr = this.spectrumList.get(i);
            if(!curr.peptide.equals(peptide1) 
                    && !curr.peptide.endsWith(peptide2)){
                System.out.println("" + curr.peptide + ": " + curr.projectedCosine(mix, fragmentTolerance));
            }
        }
    }                                                                                           //NEW CHANGE END
    
    /**
     * Given a mix spectrum directly search the spectral library using branch and bound method
     * @param mix
     * @return
     */
    public Spectrum searchAndBoundLib(Spectrum mix, double fragmentTolerance){
        Iterator<Spectrum> it =  this.iterator();
        Spectrum curr, best, top, best1, best2;
        Spectrum original = mix; //save the original mixture spectrum
        double mean=0, meanIncrease=0;
        //mix = mix.toNormVector(1, 0.5, 2000); //we make sure we normalize the mixture spectrum, since the bound makes this assumption
        double bestscore=0, nextscore = 0, upperbound = 0;
        double alpha = 1, S1=0, S2=0;
        boolean bSearchFlag = true;
        int i=0, j=0, k=0, iLimit=0, jLimit=0;
        int counts = 0; //we keep track of combo we need to create in order to track the efficacy of our bound
        //scoring all 
        it = this.iterator();
        while(it.hasNext()){
            curr = it.next();
            if(Math.abs(curr.parentMass - mix.parentMass) > this.parentMassTolerance){
                curr.score = 0.00000001;
            }else{
                curr.score = curr.cosineSim(mix, fragmentTolerance); //when spectrum are not expected to overlap much
            }                                  //direclty using cosine as criteria to filter is not that bad
        }
        Collections.sort(this.spectrumList);
        top = this.spectrumList.get(this.spectrumList.size()-1);
        //System.out.println("Spectrum: " + mix.peptide + " top single-peptide spectra: " + top.peptide + "\t" + top.score);
        i = spectrumList.size()-1;
        j = i - 1;
        k = (i + 1) / 2;
        alpha = mix.alpha(spectrumList.get(i), spectrumList.get(j), fragmentTolerance);
        //alpha = mix.residual(spectrumList.get(i));
        printDetail("new alpha: " + alpha);
        curr = new Spectrum(spectrumList.get(i), spectrumList.get(j), 1, alpha, fragmentTolerance);  //according to convention alpha < 1 thus it is used to scale the 2nd peptide
        bestscore = mix.maxScore(spectrumList.get(i), spectrumList.get(j), alpha, fragmentTolerance);//actually to be more accurated we need to do like max score check both
        //bestscore = curr.cosineSim(mix);
        best = curr;
        best1 = spectrumList.get(i);
        best2 = spectrumList.get(j);

        //complicated shit, need to draw out the actual array to get a feel of how this bound works
        //basically we update current score so far, for each of such update, we do 
        //binary search for the reduce the size of the vector that we need to search
        
        while((i >= iLimit && j >= jLimit) && i > 0 && j > 0 && i-k >= 1 && k >= 0){
            curr = new Spectrum(spectrumList.get(i), spectrumList.get(k), fragmentTolerance);
            nextscore = curr.cosineSim(mix, fragmentTolerance);
            //upperbound = (spectrumList.get(i).score + spectrumList.get(k).score) / 1.4142;  //root of two
            S1 = spectrumList.get(i).score;
            S2 = spectrumList.get(k).score;
            upperbound = (S1*S1 + S2*S2) / Math.pow((S1*S1 + S2*S2), 0.5);
            printDetail("" + i + "\t" + j + "\t" + k + "\t" + iLimit + "\t" + jLimit);
            printDetail("" + bestscore + "\t" + nextscore + "\t" +  upperbound + "\t" 
                    + spectrumList.get(i).peptide + "\t" + spectrumList.get(k).peptide);
            
            if(bestscore > upperbound){ //now best score is better than upper bound, we do not need to consider spectrum with lower score than this
                iLimit = k;
                jLimit = k;
                k = k + (j - k + 1)/2;  //reducing the search space by a factor
                
            }else{ //cannot bound the search space, try again
                k = k  - (k - iLimit + 1)/2;
            }
            
            if(k == jLimit || k == j || !bSearchFlag){      //reach searching limit or did not find a better estimate of optimal soln, try another branch
                j--;
                k = j;
                alpha = mix.alpha(spectrumList.get(i), spectrumList.get(k), fragmentTolerance); 
                //alpha = mix.residual(spectrumList.get(i));
                printDetail("new alpha: " + alpha);
                curr = new Spectrum(spectrumList.get(i), spectrumList.get(k), 1, alpha, fragmentTolerance);
                nextscore = mix.maxScore(spectrumList.get(i), spectrumList.get(k), alpha, fragmentTolerance);
                mean += nextscore;
                meanIncrease += nextscore - spectrumList.get(i).score;
                //nextscore = curr.cosineSim(mix);
                if(nextscore <= bestscore){
                    bSearchFlag = false;
                }else{
                    bSearchFlag = true;
                }
                best = nextscore > bestscore ? curr : best;
                if(nextscore > bestscore){
                    best1 = spectrumList.get(i);
                    best2 = spectrumList.get(k);
                }
                bestscore = nextscore > bestscore ? nextscore : bestscore;
            }   
            
            if(j == jLimit ){     //loop through all the j, advance i and try again
                i--;
                j = i - 1;
                k = j;
                alpha = mix.alpha(spectrumList.get(i), spectrumList.get(j), fragmentTolerance);
                //alpha = mix.residual(spectrumList.get(i));
                printDetail("new alpha: " + alpha);
                curr = new Spectrum(spectrumList.get(i), spectrumList.get(j), 1, alpha, fragmentTolerance);
                nextscore = mix.maxScore(spectrumList.get(i), spectrumList.get(j), alpha, fragmentTolerance);
                mean += nextscore;
                meanIncrease += nextscore - spectrumList.get(i).score;
                //nextscore = curr.cosineSim(mix);
                if(nextscore > bestscore){
                    bSearchFlag = true;
                }
                best = nextscore > bestscore ? curr : best;
                if(nextscore > bestscore){
                    best1 = spectrumList.get(i);
                    best2 = spectrumList.get(j);
                }
                bestscore = nextscore > bestscore ? nextscore : bestscore;
            }
            if(counts > 15000){ //we impose a limit, so we don't search everything  
                break;
            }
            counts++;
        }       
        //System.out.println("bestscore: " + bestscore);
        double bestalpha = mix.alpha(best1, best2, fragmentTolerance);
//      if(bestalpha < 1){
//          bestalpha = 1/bestalpha;
//      }
        Spectrum normBest = best.toNormVector();
        Spectrum normMix = mix.toNormVector();
        Spectrum normBest1 = best1.toNormVector();
        Spectrum normBest2 = best2.toNormVector();
        //System.out.println("cosine is: " + mix.cosineSim2(normBest) + "\t" + mix.cosineSim2(best1) + "\t" + mix.cosineSim2(best2));
        double simBias = normMix.cosineSim2(normBest, fragmentTolerance) / (bestscore+0.01);
        double score1 = normMix.cosineSim(best1, fragmentTolerance);
        double score2 = normMix.cosineSim(best2, fragmentTolerance);
        double simBias1 = normMix.cosineSim2(normBest1, fragmentTolerance) / (score1+0.01);
        double simBias2 = normMix.cosineSim2(normBest2, fragmentTolerance) / (score2+0.01);
        //System.out.println("numbers of attempts: " + counts);
        if(counts==0){
            counts=1;
        }
        try{
            //String peptide1 = new Peptide(best1.peptide+"."+best1.charge).toString();
            //String peptide2 = new Peptide(best2.peptide+"."+best2.charge).toString(); 
            String peptide1 = best1.peptide;
            String peptide2 = best2.peptide;
            //protein field in NIST and sptxt format has extra " in the beginning
            //don't want to output them
            if(best1.protein.contains("\"")){
                best1.protein = best1.protein.replaceAll("\"", "");
            }
            
            if(best2.protein.contains("\"")){
                best2.protein = best2.protein.replaceAll("\"", "");
            }
            this.out.write(this.QueryFilename + "\t" + mix.scanNumber + "\t" + peptide1 + " \t" + peptide2 + "\t"   
                    +  best1.protein + "\t" + best2.protein +"\t" + best1.charge + "\t" + best2.charge 
                    + "\t" + bestscore + "\t" +  score1 + "\t" + score2 + "\t" + best1.cosineSim(best2, fragmentTolerance)
                    + "\t" + bestalpha + "\t" + mix.residual(best1, fragmentTolerance)
                    + "\t" + mix.explainedIntensity(0.85) + "\t" + simBias
                    + "\t" + simBias1  + "\t" + simBias2 
                    + "\t" + mix.projectedCosine(best, fragmentTolerance) + "\t" + mix.projectedCosine(best1, fragmentTolerance) + "\t" + mix.projectedCosine(best2, fragmentTolerance)
                    + "\t" + mean/counts + "\t" + meanIncrease/counts
                    + "\t" + mix.parentMass + "\t" + best1.parentMass + "\t" + best2.parentMass +"\t" 
                    + mix.specIndex + "\t" + best1.specIndex + "\t" + best2.specIndex + "\n"
                );
        }catch(IOException ioe){
            System.err.println("error writing to outputFile");
            System.err.println(ioe.getMessage());
            return null;
        }
        best.score = bestscore;
        return best;
    }
    
    /**
     * Same as searchAndBoundLib, but instead of using a binary search method to search for the bound
     * it simply scan the list of library spectra by decreasing cosine score, it is straightforward,
     * but slightly less efficient
     * @param mix
     * @return
     */
    public Spectrum searchAndBoundLib2(Spectrum mix, double fragmentTolerance){  //with custom candidate pairs generations
        Iterator<Spectrum> it =  this.iterator();
        Spectrum curr, curr1, curr2, best, top, best1, best2;
        Spectrum original = mix; //save the original mixture spectrum
        // mix = mix.toNormVector(1, 0.5, 2000); //we make sure we normalize the mixture spectrum, since the bound makes this assumption
        double bestscore=0, nextscore = 0, upperbound = 0;
        double alpha = 1, S1=0, S2=0;
        double bestalpha = 0;
        int counts = 0; //we keep track of combo we need to create in order to track the efficacy of our bound
        //preprocessing steps, find top score pair, record upperbound
        //scoring all 
        it = this.iterator();
        while(it.hasNext()){
            curr = it.next();
            curr.score = curr.cosineSim(mix, fragmentTolerance); //when spectrum are not expected to overlap much
                                              //direclty using cosine as criteria to filter is not that bad
        }
        Collections.sort(this.spectrumList);
        top = this.spectrumList.get(this.spectrumList.size()-1);
        it = this.iterator();
        while(it.hasNext()){
            curr = it.next();
            S1 = top.score;
            S2 = curr.score;
            upperbound = (S1*S1 + S2*S2) / Math.pow((S1*S1 + S2*S2), 0.5);
            curr.upperBound = upperbound;  //recorded upperbound
            //System.out.println("upperboud: " + curr.upperBound);
        }
        
        //candidate selection using custom filter
        it = this.iterator();
        while(it.hasNext()){
            curr = it.next();
            curr.score = curr.projectedCosine(mix, fragmentTolerance); //sort using custom filter
        }
        Collections.sort(this.spectrumList);
        
        List<Spectrum> candidates = new ArrayList();
        for(int i = this.spectrumList.size()-1; i > spectrumList.size()-this.topSpectraKept; i--){
            candidates.add(this.spectrumList.get(i));
        }
        System.out.println(original.peptide);
        it = this.iterator();
        //put score back to compute upper bound
        while(it.hasNext()){
            curr = it.next();
            curr.score = curr.cosineSim(mix, fragmentTolerance); //sort using custom filter
        }
        
        //start searching sort candidates by cosine
        Collections.sort(candidates);
//      System.out.println("best score: " + candidates.get(candidates.size()-1).score);
//      System.out.println("worst score: " + candidates.get(0).score);
        //candidates = this.spectrumList;
        best1 = null; best = null; best2 = null;
        for(int i = candidates.size()-1; i > 0; i--){
            curr1 = candidates.get(i);
//          if(curr1.upperBound < bestscore){
//              break;
//          }
            if(counts > 50000){
                break;
            }
            for(int j = i-1; j > 0; j--){
                curr2 = candidates.get(j);
                if(counts > 50000){
                    break;
                }
                S1 = curr1.score;
                S2 = curr2.score;
                upperbound = (S1*S1 + S2*S2) / Math.pow((S1*S1 + S2*S2), 0.5);
                curr2.upperBound = upperbound;  //recorded upperbound
                if(curr2.upperBound > bestscore){
                    alpha = mix.alpha(curr1, curr2, fragmentTolerance);
                    //System.out.println("alpha is: " + alpha);
                    nextscore = mix.maxScore(curr1, curr2, alpha, fragmentTolerance);
                    //System.out.println("current score: " + nextscore);
                    if(nextscore > bestscore){
                        bestscore = nextscore;
                        best1 = candidates.get(i);
                        best2 = candidates.get(j);
                        bestalpha = alpha;
                        best = new Spectrum(candidates.get(i), candidates.get(j), 1, alpha, fragmentTolerance);
                    }
                    //System.out.println("i: " + i + " j: " + j);
                    counts++;
                }
            }
        }
        if(bestalpha < 1){
            bestalpha = 1/bestalpha;
        }
        System.out.println("numbers of attempts: " + counts);
        System.out.println("best answer: " + mix.peptide + "\t" + best1.peptide + " and " + best2.peptide + "\t"   
            +   "\t" + mix.parentMass + "\t" + best1.parentMass + "\t" + best2.parentMass
            + "\t" + bestscore + "\t" +  best1.score + "\t" + best2.score + "\t" + best1.cosineSim(best2, fragmentTolerance)
            + "\t" + bestalpha + "\t" + mix.residual(best1, fragmentTolerance));
        best.score = bestscore;
        return best;
    }
    
    
    private Spectrum findAnswerInSortedSpectrumList(Spectrum mix, String mixedpeptides){
        return findAnswerInSortedSpectrumList(mix, mixedpeptides, this.spectrumList);
    }
    
    //find where the answer lies in our sorted list
    private Spectrum findAnswerInSortedSpectrumList(Spectrum mix, String mixedpeptides, List<Spectrum> sortedList){
        String[] peps = mixedpeptides.split(" & ");
        Spectrum curr, s1= new Spectrum(), s2=s1;
        int index = 0;
        if(peps.length == 1){
            return new Spectrum(); //we return a dummy when mix is not a mix
        }
        for(index = sortedList.size()-1; index >=0; index--){
            curr = sortedList.get(index);
            if(curr.peptide.equals(peps[0])){
                s1 = curr;
                System.out.println(peps[0] + " @ pos " + index + " score: " +  curr.score);
                System.out.println(peps[0] + " sharing " + mix.shareSim(s1, 1)*100 + "% peaks with mixture");
                printDetail(curr.toString());
            }
            if(curr.peptide.equals(peps[1])){
                s2 = curr;
                System.out.println(peps[1] + " @ pos " + index + " score: " + curr.score);
                System.out.println(peps[1] + " sharing " + mix.shareSim(s2, 1)*100 + "% peaks with mixture");
                printDetail(curr.toString());
            }
            
        }
        //curr = this.mixSpectrumSimilarity(mix, peps[0], peps[1]);
//      double alpha = mix.alpha(s1, s2);
//      System.out.println("real estimate of alpha is: "  + alpha);
        curr = new Spectrum(s1, s2, 10, 1, fragmentTolerance);
//      curr= curr.toNormVector(1, 0.5, 2000);
//      System.out.println(curr);
//      System.out.println("the correct mixed: " + curr.cosineSim(mix));
        return curr;
    }
    
    public void maxSolution(Spectrum mix, double fragmentTolerance) {   
        List<Spectrum> spectra = getAllSpectrums() ;
        TreeMap <Double, Spectrum> sortedSpectra = new TreeMap() ;
        Hashtable answer = new Hashtable() ;
        Vector <Spectrum> v;
        Spectrum temp;
        SpectrumLib newLib, fin ;
        String realPeptide[] = mix.peptide.split(" & ") ;
        String sol = realPeptide[1] + " & " + realPeptide[0] ;
        int i = 0;
        long size, rank = 1;
        boolean flag1, flag2 ;      
        size = spectra.size();
                
        flag1 = flag2 = false ;
        
        // Filtering
        
        // calculate the projected cosines and store it
        // in a treemap (ordered)
        double score1 = 0, score2 = 0;                      
        for(int j = 0; j < spectra.size(); j++) {
            
            temp = (spectra.get(j)).toNormVector(1,0.5,2000) ; 

            if (temp.projectedCosine(mix, fragmentTolerance) > 0) {

                sortedSpectra.put(new Double(temp.projectedCosine(mix, fragmentTolerance)), temp) ; 
        
                if (temp.peptide.equals(realPeptide[0])){
                    score1 = temp.cosineSim(mix, fragmentTolerance);
                    flag1 = true ;
                    
                }
                    
                if (temp.peptide.equals(realPeptide[1])) {
                    score2 = temp.cosineSim(mix, fragmentTolerance);
                    flag2 = true;
                    
                }
            
            }
        }
        
        if (flag1 && flag2)
            System.out.print("2\t") ;
        else if (flag1 || flag2)
            System.out.println("1\t0") ;
        else
            System.out.println("0\t-1") ;
        
        
        // Get best scores (above threshold)
        
        if (flag1 && flag2) {
            
            spectra = new Vector() ;    
            
            for (int k = 0 ; k < 500 && sortedSpectra.size() > 0; k++)
                    spectra.add(sortedSpectra.remove(sortedSpectra.lastKey())) ; 
                    
            sortedSpectra = null ;
        
            Spectrum A = null ;
            Spectrum B = null ;
            
            TreeMap solutions = new TreeMap() ;
            double alpha ;
            
            for (int k = 0; k < spectra.size();k++) {
                A = spectra.get(k) ;
                
                for(int l = k + 1; l < spectra.size();l++){
                    
                    if ( !(((spectra.get(k)).peptide).equals(((spectra.get(l)).peptide))) ) {
                        
                        B = spectra.get(l) ;
                        //alpha = mix.alpha(A, B) ;
                        alpha = 1;
                        solutions.put(new Double(mix.maxScore(A, B, alpha, fragmentTolerance)), A.peptide + " & " +B.peptide) ;
                        //solutions.put(new Double(mix.inverseMaxScore(A, B, alpha)), A.peptide + " & " + B.peptide) ;
                    }
                                    
                }
            }

            String pep ;
            Double key;
            size = solutions.size() ;
                        
            if (!solutions.isEmpty()) {
                key = (Double)solutions.lastKey();
                pep = (String)solutions.remove(solutions.lastKey()) ;
                
                while (!mix.peptide.equals(pep) && !sol.equals(pep) && !solutions.isEmpty()) {
                    key = (Double)solutions.lastKey();      
                    pep = (String)solutions.remove(solutions.lastKey()) ;
                    rank++ ;
            
                }   
                    
                if (mix.peptide.equals(pep) || sol.equals(pep)){
                    
                    System.out.print(rank) ;
                    System.out.println("\t" + key.toString() + "\t" + score1 + "\t" + score2);
                }
                else
                    System.out.println("0\t0\t0\t") ;
            }
                
            
        }
        
    } // end of function

/**
 * First filtered spectral library by projected-cosine then find the max pair
 * @param mix
 * @return
 */
public Spectrum filterAndSearch(Spectrum mix, double fragmentTolerance) {
    
    List <Spectrum> spectra = getAllSpectrums() ;
    TreeMap <Double, Spectrum> sortedSpectra = new TreeMap() ;
    Hashtable answer = new Hashtable() ;
    Vector <Spectrum> v;
    Spectrum temp;
    int i = 0, size = spectra.size();
    Spectrum top = this.topSpectrum(mix, fragmentTolerance);
    System.out.println("top single-peptide: " + top.peptide + "\t" + top.score);
    // Filtering
    
    // calculate the projected cosines and store it
    // in a treemap (ordered)
    double score1 = 0, score2 = 0;                      
    for(int j = 0; j < spectra.size(); j++) {
        
        temp = (spectra.get(j));
        
        if (temp.projectedCosine(mix, fragmentTolerance) > 0) {

            sortedSpectra.put(new Double(temp.projectedCosine(mix, fragmentTolerance)), temp) ; 
        }
    }
    
    // Get best scores (above threshold)
    spectra = new Vector() ;    
    for (int k = 0 ; k < 100 && sortedSpectra.size() > 0; k++)
                spectra.add(sortedSpectra.remove(sortedSpectra.lastKey())) ; 
                
    sortedSpectra = null ;
    Spectrum A = null ;
    Spectrum B = null ;
    Spectrum best = new Spectrum();
    Double bestscore = null, currscore = null;
    TreeMap solutions = new TreeMap() ;
    double alpha = 1, best1 =0, best2 =0;
    double bestoverlap = 0;
    for (int k = 0; k < spectra.size();k++) {
        A = spectra.get(k) ;
        for(int l = k + 1; l < spectra.size();l++){
            if ( !(((spectra.get(k)).peptide).equals(((spectra.get(l)).peptide))) ) {       
                B = spectra.get(l) ;
                alpha = mix.alpha(A, B, fragmentTolerance) ;
                //alpha = 1;
                currscore = new Double(mix.maxScore(A, B, alpha, fragmentTolerance));
                solutions.put(currscore, A.peptide + " & " +B.peptide) ;
                //solutions.put(new Double(mix.inverseMaxScore(A, B, alpha)), A.peptide + " & " + B.peptide) ;
                bestscore = (Double)solutions.lastKey();
                if(bestscore.equals(currscore)){
                    best = new Spectrum(A, B, alpha, 1, fragmentTolerance);
                    best.score = bestscore.doubleValue();
                    best1 = A.score;
                    best2 = B.score; //the score were calculated when we calculate top single spectrum match to the mixture spectrum
                    bestoverlap = A.cosineSim(B, fragmentTolerance);
                }
            }
                                
        }
    }
    System.out.println(mix.peptide + "\tPutative Pair: " + best.peptide + "\t" + best1 + "\t" + best2 + "\toverlapping: " + bestoverlap + 
            "\tcombined: " + best.score + "\tMassdiff: " + best.modMass);
    return best;
} 
    /**We iteratively identify one component of the mixture 
     * substract that component and repeat 
     * @param mix
     * @param iters
     * @return
     */
    public Vector<Spectrum> peals(Spectrum mix, int iters, double fragmentTolerance){
        Vector<Spectrum> answers = new Vector();
        Vector<Spectrum> tempLibrary = new Vector();
        double[] bestscores = new double[iters];
        Iterator<Spectrum> it =  null;
        Spectrum curr, curr2, top;
        top = this.topSpectrum(mix, fragmentTolerance);
        System.out.println("top single-peptide: " + top.peptide + "\t" + top.score);
        tempLibrary.addAll(this.spectrumList);
        for(int i = 0; i < iters; i++){
            it =  this.iterator();
            while(it.hasNext()){
                curr = it.next();
                curr.score = curr.cosineSim(mix, fragmentTolerance); //when spectrum are not expected to overlap much
                //curr.score = curr.projectedCosine(mix);                                  //direclty using cosine as criteria to filter is not that bad
            }
            Collections.sort(tempLibrary);
            bestscores[i] = tempLibrary.lastElement().score;
            mix.subtractSharePeaks(tempLibrary.lastElement(), fragmentTolerance);
            answers.add(tempLibrary.lastElement());
            findAnswerInSortedSpectrumList(mix, mix.peptide, tempLibrary);
            //tempLibrary.remove(tempLibrary.lastElement()); //remove the last elements

        }
        
        //we reinsert the proper scores back to the top hits
        for(int i = 0; i < iters; i++){
            answers.get(i).score = bestscores[i];
        }   
        return answers;
    }
    
    /**J
     * Search the best matched spectrum by using cosine similarity
     * @param mix
     * @return
     */
    public Spectrum topSpectrum(Spectrum mix, double fragmentTolerance){
        Iterator<Spectrum> it =  null;
        Spectrum curr, curr2;
        it =  this.iterator();
        while(it.hasNext()){
            curr = it.next();
            if(Math.abs(curr.parentMass - mix.parentMass) > this.parentMassTolerance){
                curr.score = 0.00000001;
            }else{
                curr.score = curr.cosineSim(mix, fragmentTolerance); //when spectrum are not expected to overlap much
            }                                  //direclty using cosine as criteria to filter is not that bad
        }
        Collections.sort(this.spectrumList);
        System.out.println(mix.peptide + ":\t"  + spectrumList.get(spectrumList.size()-1).peptide);
        return this.spectrumList.get(this.spectrumList.size()-1);
    }
    
    //for a given peptide randomly selected one of the spectrum
    //correspond to this peptide
    public Spectrum getRandomSpectrum(String peptide){
        List <Spectrum> v = this.spectrumLibrary.get(peptide);
        return v.get((int)(Math.random() * (v.size() - 1)));
    }
    
    public Spectrum getRandomSpectrum(){
        //System.out.println("size is: " + spectrumList.size());
        int ind = (int)(Math.random() * (spectrumList.size()-1));
        //System.out.println("index is: " + ind);
        return spectrumList.get(ind);
    }
    public void printLib(){
        Iterator it = this.iterator();
        while(it.hasNext()){
            System.out.println(it.next());
            System.out.println();
        }
    }
    
    //we print out the number of peaks in each spectrum in this library
    //so we can see the number of peaks distribution
    public double getPeakCounts(double fraction){
        Iterator<Spectrum> it = this.iterator();
        Spectrum s = null;
        double average = 0;
        while(it.hasNext()){
            s = it.next();
            average += s.numberOfPeaks(fraction);
            System.out.println("peptide: " + s.peptide.length() 
                    + "\tpeaks: " + s.numberOfPeaks(fraction));
        }
        return average / (this.spectrumList.size());
    }
    

    public void removeDuplicate(double overlap){
        Spectrum s1, s2;
        List<Spectrum> l = this.getAllSpectrums();
        List<Spectrum> copy = new Vector();
        for(int i = 0; i < l.size(); i++){
            s1 = l.get(i);
            for(int j = i+1; j < l.size(); j++){
                s2 = l.get(j);
                if(Math.abs(s1.parentMass - s2.parentMass) < 1.8 && s1.cosineSim(s2, fragmentTolerance) > overlap){
                    if(s1.sumOfPeaks() > s2.sumOfPeaks()){
                        copy.add(s2);
                    }else{
                        copy.add(s1);
                    }
                }
            }
        }
        for(int i = 0; i < copy.size(); i++){
            this.spectrumLibrary.remove(copy.get(i).peptide);
        }
        this.spectrumList = this.getAllSpectrums();
        
    }
    
    //this is a reverse-direction search, given two spectrum that are likely to formed a mixture,
    //this method search the library of mixture to identify potential mixture spectrum that are
    //made of the two single-peptide spectrum
    public Spectrum reverseSearchMix2(Spectrum s1, Spectrum s2, double fragmentTolerance){
        Iterator<Spectrum> it = this.iterator();
        Spectrum mix = null, best = null;
        double alpha, bestscore = 0;
        while(it.hasNext()){
            mix = it.next();
            mix = mix.toNormVector();
            //alpha = 1.0; //just dummy alpha, shiftMaxScore is not using it in current version
            alpha = mix.alpha(s1, s2, fragmentTolerance);
            mix.score = mix.shiftMaxScore(s1, s2, alpha, fragmentTolerance);
            //mix.score = mix.maxScore(s1, s2, alpha);
            //System.out.println("mix: " + mix.peptide);
            //System.out.println("alpha is: " + alpha);
            //System.out.println("max score is: " + mix.score);
            if(mix.score > bestscore){
                bestscore = mix.score;
                best = mix;
            }
        }
//      if(best.shiftCosineSim(s1) - s1.shiftCosineSim(best) > 0.001){
//          System.out.println("shiftCosine seems not communitive");
//      }
//      System.out.println("best: " + best.peptide + " " + s1.peptide + " & " + s2.peptide +  
//              " parentmass: " + best.parentMass + " combined: " + best.score + 
//              " single: " + best.shiftCosineSim(s1) + " & " + best.shiftCosineSim(s2) +
//              " overlap: " + s1.shiftCosineSim(s2)
//      );
        Spectrum ns1, ns2;
        ns1 = s1.toNormVector();
        ns2 = s2.toNormVector();
        alpha = mix.alpha(ns1, ns2, fragmentTolerance);
        System.out.println("best answer: " + best.peptide + " " + s1.peptide + " and " + s2.peptide +  "\t"
                + best.parentMass + "\t" + s1.parentMass + "\t" + s2.parentMass +  "\t" 
                + best.score + "\t" + s1.shiftCosineSim(best, fragmentTolerance) + " \t " + s2.shiftCosineSim(best, fragmentTolerance) + "\t " + s1.shiftCosineSim(s2, fragmentTolerance) + "\t"
                + alpha + "\t" + 1/best.residual(ns1, fragmentTolerance)
        );
        return best;
    }
    public Spectrum reverseSearchMix(Spectrum s1, Spectrum s2, double fragmentTolerance){
        Iterator<Spectrum> it = this.iterator();
        Spectrum mix = null, best = null;
        double alpha, bestscore = 0, bestalpha =0;
        while(it.hasNext()){
            mix = it.next();
            //mix = mix.toNormVector();
            alpha = mix.alpha(s1, s2, fragmentTolerance);
            mix.score = mix.maxScore(s1, s2, alpha, fragmentTolerance);
            if(mix.score > bestscore){
                bestscore = mix.score;
                best = mix;
                bestalpha = alpha;
            }
        }
        
        System.out.println("best answer: " + best.peptide + " " + s1.peptide + " and " + s2.peptide +  "\t"
                + best.parentMass + "\t" + s1.parentMass + "\t" + s2.parentMass +  "\t" 
                + best.score + "\t" + s1.cosineSim(best, fragmentTolerance) + " \t " + s2.cosineSim(best, fragmentTolerance) + "\t " + s1.shiftCosineSim(s2, fragmentTolerance) + "\t"
                + bestalpha + "\t" + 1/best.residual(s1, fragmentTolerance)
        );
        return best;
    }
    
    /**
     * A static method to write a spectrum lib to a file
     * @param outfile
     * @param s
     */
    //note this write the object itself not the string  of the lib
    public static void writeLibToFile(String outfile, SpectrumLib s){
        try{
            BufferedOutputStream bo = new BufferedOutputStream(new FileOutputStream(outfile));
            ObjectOutputStream oo = new ObjectOutputStream(bo);
            oo.writeObject(s);
            oo.flush();
            oo.close();
        }catch(IOException ioe){
            System.out.println(ioe.getMessage());
            
        }
    }
    
    public void printLibToFile(String outfile, SpectrumLib s){
        try{
            BufferedWriter bo = new BufferedWriter(new FileWriter(outfile));
            Iterator<Spectrum> it = this.iterator();
            Spectrum curr;
            while(it.hasNext()){
                bo.write(it.next().toString());
                bo.write("\n");
            }
            bo.flush();
            bo.close();
        }catch(IOException ioe){
            System.out.println(ioe.getMessage());
            
        }
    }
    
    public void printLibToDTAFile(String outfile, SpectrumLib s){
        try{
            Iterator<Spectrum> it = this.iterator();
            Spectrum curr;
            while(it.hasNext()){
                curr = it.next();
                String filename = curr.spectrumName.replaceAll("\\W+", "_");
                String out = outfile+filename+".dta";
                System.out.println("printing to file " + out);
                BufferedWriter bo = new BufferedWriter(new FileWriter(out));
                bo.write("<precursor>\n");
                bo.write((curr.parentMass*curr.charge-(curr.charge-1)*Mass.PROTON_MASS)
                    +" 1000 "+curr.charge+"\n"); //default intensity
                bo.write("</precursor>\n");
                List<Peak> pList = curr.getPeak();
                for(int i = 0; i < pList.size();i++){
                    Peak p = pList.get(i);
                    bo.write(p.getMass() + " " + p.getIntensity()+"\n");
                }
                bo.flush();
                bo.close();
            }
        }catch(IOException ioe){
            System.out.println(ioe.getMessage());
            
        }
    }

    //read the library back from file
    public static SpectrumLib readLibFromFile(String infile){
        SpectrumLib s = null;
        try{
            BufferedInputStream bi = new BufferedInputStream(new FileInputStream(infile));
            ObjectInputStream oi = new ObjectInputStream(bi);
            s = (SpectrumLib)oi.readObject();
            oi.close();
        }catch(Exception ioe){
            System.out.println(ioe.getMessage());
            System.out.println(ioe.getCause());
        }
        return s;
    }
    
    public static void annoateSpectrumFromInspectFile(String file){
        //to be implememtned
    }
    
    public static String parseFileFormat(String filename){
        return filename.split("\\.")[1];
    }
    
    public static void runSearch(String libFile, String mixtureFile, double massTolerance, String resultFile, double parentMassTolerance){
        SpectrumLib lib1 = new SpectrumLib(libFile, parseFileFormat(libFile));
        SpectrumLib mixlib = new SpectrumLib(mixtureFile, parseFileFormat(mixtureFile));
        lib1.toNormVector(massTolerance, -1*massTolerance/2, 2000);
        lib1.normIntensity();
        mixlib.scaleSpectrumMass(0.9995);
        mixlib.toNormVector(massTolerance, -1*massTolerance/2, 2000);
        mixlib.normIntensity();
        lib1.parentMassTolerance = parentMassTolerance;
        for(Iterator<Spectrum> mixIter = mixlib.iterator(); mixIter.hasNext();){
            Spectrum mixture = mixIter.next();
            lib1.searchAndBoundLib2(mixture, massTolerance);
        }
    }
    
    public static void runSearch(SpectrumLib lib, SpectrumLib mixtureLib, double massTolerance, String resultFile, double parentMassTolerance){
        lib.toNormVector(massTolerance*2, massTolerance,  2000);
        lib.normIntensity();
        lib.parentMassTolerance = massTolerance;
        //lib.toNormVector(massTolerance*2, massTolerance,  2000);
        mixtureLib.windowFilterPeaks(10, 25);
        mixtureLib.toNormVector(massTolerance*2, massTolerance, 2000);
        mixtureLib.normIntensity();
        mixtureLib.toNormVector(massTolerance*2, massTolerance, 2000);
        for(Iterator<Spectrum> mixIter = mixtureLib.iterator(); mixIter.hasNext();){
            Spectrum mixture = mixIter.next();
            lib.searchAndBoundLib(mixture, massTolerance);
        }
    }

    
    public double getParentMassTolerance() {
        return parentMassTolerance;
    }

    public void setParentMassTolerance(double parentMassTolerance) {
        this.parentMassTolerance = parentMassTolerance;
    }
    
    //create a dummy spectrum 
    public static Spectrum createDummySpect(){
        Spectrum dummy = new Spectrum();
        dummy.getPeak().add(new Peak(300, 10));
        dummy.getPeak().add(new Peak(500, 10));
        dummy.getPeak().add(new Peak(1000, 10));
        dummy.peptide = "AAAAA";
        dummy.protein = "DECOY_RANDOM_HITS";
        return dummy;
    }
    /**
     * Search a query spectral file against a spectral library
     * This is the main entry point to the M-SPLIT tool
     * @param lib
     * @param query spect File
     * @param parentMassTolerance
     * @param massTolerance
     * @param resultFile
     * @return
     */
    public static boolean runSearch(SpectrumLib lib, String mixtureLibFile, double parentMassTolerance, double massTolerance, String resultFile, int guessCharges){
        lib.scaleSpectrumMass(0.9995);
        lib.toNormVector(massTolerance*2, massTolerance,  2000);
        lib.normIntensity();
        lib.parentMassTolerance = parentMassTolerance;
        lib.toNormVector(massTolerance*2, massTolerance,  2000);
        Iterator it;
        Spectrum.BINWIDTH = massTolerance*2;
        Spectrum.MINMASS = massTolerance;
        
        try{
                it = new SortedSpectrumReader(mixtureLibFile);
        }catch(Exception e){
            System.err.println("error reading query spectrum file");
            System.err.println(e.getMessage());
            e.printStackTrace();
            return false;
        }
        lib.QueryFilename = mixtureLibFile;
        int count = 1, success=0;
        long start = (new GregorianCalendar()).getTimeInMillis();
        
        String libraryObjectFile = lib.Filename;
        CandidateSpectrumGenerator gen = new CandidateSpectrumGenerator();
        gen.loadLibraryFromFile(libraryObjectFile);

        lib.initOutPut(resultFile);
        try{
            lib.out.write("#SpectrumFile\tScan#\tAnnotation1\tAnnotation2\tProtein1\tProtein2\tCharge1\tCharge2\t" +
                "cosine(M, A+B)\tcosine(M,A)\tcosine(M,B)\tcosine(A,B)\talpha\tres-alpha" +
                "\t#Peak-0.85Intensity\tsimBias(A)\tsimBias(B)" +
                "\tprojCos(M,A+B)\tprojCos(M,A)\tprojCos(M,B)\tmeanCos\tmeanDeltaCos\tPrecursor(M)\tPrecursor(A)\tPrecursor(B)\tindex\tlibInd1\tlibIndex2\n");
        }catch(IOException ioe){
            System.err.println("error writing to outputfile");
            System.err.println(ioe.getMessage());
            return false;
        }
        for(Iterator<Spectrum> mixIter = it; mixIter.hasNext();){
            try{
                Spectrum mixture = mixIter.next();
                //if spectrum has too few peaks, we skip it as it has
                //small chance of finding any significant identification
                if(mixture.getPeak().size() < 5){
                    continue; 
                }
                mixture.scaleMass(0.9995);
                mixture.windowFilterPeaks(10, 25);
                mixture = mixture.toNormVector(massTolerance*2, massTolerance, 2000);
                mixture.sqrtSpectrum();
                mixture = mixture.toNormVector(massTolerance*2, massTolerance, 2000);
                if (guessCharges == 1) {
                    lib.spectrumList = gen.getSpectraByMass(mixture, lib.parentMassTolerance);
                }
                else {
                		lib.spectrumList = gen.getSpectraByMassAndCharge(mixture, lib.parentMassTolerance);
                }
                if(lib.spectrumList.size() < 1){
                    //System.out.println(mixture.spectrumName + "\t" + mixture.parentMass + "\t" + mixture.charge + "number of candidates: " + lib.spectrumList.size());
                    count++;
                    continue;
                }else if(lib.spectrumList.size() == 1){
                    Spectrum dummy = createDummySpect();
                    dummy.parentMass = mixture.parentMass;
                    lib.addSpectrum(dummy);
                }
                mixture.spectrumName = ""+count;
                Spectrum best = lib.searchAndBoundLib(mixture, massTolerance);
                if(best != null) success++;
                count++;
            }catch(Exception e){
                System.err.println("error searching spectrum");
                System.err.println(e.getMessage());
                e.printStackTrace();
            }
            if(count % 1000 == 0){
                System.out.println("Finish searching: " + count);
            }
        }
        //finalizing output
        try{
            lib.out.flush();
            lib.out.close();
        }catch(IOException ioe){
            System.err.println("error occur when finalizing output");
            return false;
        }
        System.out.println("matching " + count + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
        return success > -1;
    }
    
    public static boolean runSearch(SpectrumLib lib, String mixtureLibFile, double parentMassTolerance, String resultFile, int guessCharges){
        return runSearch(lib, mixtureLibFile, parentMassTolerance, 0.5, resultFile, guessCharges);
    }
    
    public static void runSearch(SpectrumLib lib, List<Spectrum> mixtureLib, double massTolerance, String resultFile, double parentMassTolerance){
        lib.toNormVector(massTolerance*2, massTolerance,  2000);
        lib.normIntensity();
        lib.toNormVector(massTolerance*2, massTolerance,  2000);
        lib.parentMassTolerance = parentMassTolerance;
        for(Iterator<Spectrum> mixIter = mixtureLib.iterator(); mixIter.hasNext();){
            Spectrum mixture = mixIter.next();
            mixture.windowFilterPeaks(10, 25);
            mixture = mixture.toNormVector(massTolerance*2, massTolerance, 2000);
            mixture.sqrtSpectrum(); 
            mixture = mixture.toNormVector(massTolerance*2, massTolerance, 2000);
            lib.searchAndBoundLib(mixture, massTolerance);
            //lib.topSpectrum(mixture);
        }
    }   
    

    public static void main(String[] args){
        if(args.length < 5){
            //System.err.println("java -jar MSPLIT.jar <library file> <query spectrum file> <precursor mass tolerance> <outputfile>");
            throw new IllegalArgumentException("java -jar MSPLIT.jar <library file> <query spectrum file> <precursor mass tolerance> <ion tolerance> <outputfile> [-guessCharges 0/1]");
            
        }else{
        //  args[0]= "../mixture_linked/MSPLIT_v1.0/test_Blibospec_libs/Acidiphilium_cryptum_JF-5.ms2";
        //  args[1]= "../mixture_linked/yeast_data/klc_010908p_yeast-digest.mzXML";
        //  args[2] = "3";
        //  args[3] = "../mixture_linked/out.txt";
            String filename = args[0];
            String fileMix = args[1];
            String fileout = args[4];
            int guessCharges = 1;
            if (args.length > 5) {
            		String opt = args[5];
            		if (opt.equalsIgnoreCase("-guessCharges") || opt.equalsIgnoreCase("-c")) {
            			guessCharges = Integer.parseInt(args[6]);
            		}
            		else {
            			throw new IllegalArgumentException("Unrecognized argument: " + args[5]);
            		}
            }
            SpectrumLib lib1 = null;
            double parentmassTol = 3;
            double fragmentTolerance = 0.5;
            try{
                if(filename.contains(".map")){
                    lib1 = new SpectrumLib();
                    lib1.Filename = filename; 
                }else {
                    lib1 = new SpectrumLib();
                    lib1.Filename= CandidateSpectrumGenerator.indexSpectralLibrary(filename, fragmentTolerance);
                }
            }catch(Exception e){
                System.err.println("Error loading spectral library: " + filename);
                System.err.println(e.getMessage());
                e.printStackTrace();
                System.exit(1); 
            }
            try{
                parentmassTol = Double.parseDouble(args[2]);
            }catch(NumberFormatException ne){
                System.err.println("precursor tolerance not valide: " + args[2]);
                System.err.println(ne.getMessage());
                System.exit(1);
            }
            try{
                fragmentTolerance = Double.parseDouble(args[3]);
            }catch(NumberFormatException ne){
                System.err.println("fragment tolerance not valid: " + args[3]);
                System.err.println(ne.getMessage());
                System.exit(1);
            }
            //lib1.scaleSpectrumMass(0.9995);
            lib1.outputFile = fileout;
            if(!runSearch(lib1, fileMix, parentmassTol, fragmentTolerance, fileout, guessCharges)){
                System.err.println("Error occur during M-SPLIT searchess");
                System.exit(1);
            }
        }
    }   
    
    
    private static double[] cosineStat(Vector v1, Vector v2, double fragmentTolerance){
        double mean = 0, var = 0;
        double[] cosines = new double[v1.size()*v2.size()];
        int k = 0;
        Spectrum s1, s2;
        for(int i = 0; i < v1.size(); i++){
            s1 = (Spectrum)v1.get(i);
            for(int j = 0; j < v2.size(); j++){
                s2 = (Spectrum)v2.get(j);
                cosines[k] = s1.cosineSim(s2, fragmentTolerance);
                mean += cosines[k];
                k++;
            }
        }
        mean = mean / k;
        for(int i = 0; i < cosines.length; i++){
            var += Math.pow(cosines[i]-mean, 2);
        }
        double[] ret = {mean, var};
        return ret;
    }

    public List<Spectrum> getSpectrumList() {
        return spectrumList;
    }

    public void setSpectrumList(Vector<Spectrum> spectrumList) {
        this.spectrumList = spectrumList;
    }

    public Map<String, List<Spectrum>> getSpectrumLibrary() {
        return spectrumLibrary;
    }

    public void setSpectrumLibrary(
            Hashtable<String, List<Spectrum>> spectrumLibrary) {
        this.spectrumLibrary = spectrumLibrary;
    }
    
}

