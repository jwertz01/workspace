package misc;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.PrintStream;
import java.math.BigDecimal;
import java.math.MathContext;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Random;

import org.systemsbiology.jrap.stax.MSXMLSequentialParser;

import msdbsearch.DBScanner;
import msgf.MSGFDBResultGenerator;
import msutil.AminoAcid;
import msutil.AminoAcidSet;
import msutil.Composition;
import msutil.CompositionFactory;
import msutil.Enzyme;
import msutil.IonType;
import msutil.Peptide;
import msutil.SpecKey;
import msutil.SpectraIterator;

import parser.BufferedLineReader;
import parser.BufferedRandomAccessLineReader;
import parser.MgfSpectrumParser;
import parser.MzXMLSpectraIterator;
import parser.MzXMLSpectraMap;

public class Chores {
	public static void main(String argv[]) throws Exception
	{
//		System.out.println(Enzyme.TRYPSIN.getProbCleavageSites());
//		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
//		DBScanner.setAminoAcidProbabilities("/home/sangtaekim/Research/Data/IPI/IPI_human_3.79.fasta", aaSet);
//		System.out.println(Enzyme.TRYPSIN.getProbCleavageSites());
//		System.out.format("%.2f", 10.346f);
//		System.out.println(new IonType.PrefixIon(1,19));
//		long time = System.currentTimeMillis();
//		float value = 1;
//		for(int i=0; i<1000000000; i++)
//			for(int j=0; j<1000000000; j++)
//				value += Math.log(2.3f);
//		System.out.println(System.currentTimeMillis()-time);
//		System.out.println(IonType.X.getOffset());
//		System.out.println(Character.isUpperCase(']'));
//		System.out.println(Float.MIN_NORMAL);
//		System.out.println(Double.MIN_NORMAL);
//		System.out.println(Enzyme.LysN.getPeptideCleavageEfficiency());
//		System.out.println(Enzyme.LysN.getNeighboringAACleavageEffiency());
//		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
//		aaSet.registerEnzyme(Enzyme.LysC);
//		System.out.println(aaSet.getNeighboringAACleavageCredit());
//		System.out.println(aaSet.getNeighboringAACleavagePenalty());
//		System.out.println(aaSet.getPeptideCleavageCredit());
//		System.out.println(aaSet.getPeptideCleavagePenalty());
//		efdrTest();
//		System.out.println((Composition.N15-Composition.N)*2+(Composition.C13-Composition.C)*6);
//		System.out.println((Composition.N15-Composition.N)*4+(Composition.C13-Composition.C)*6);
//		System.out.println(Composition.OFFSET_Y);
//		System.out.println(Runtime.getRuntime().availableProcessors());
//		System.out.println(Composition.getMass("H-2O-1"));
//		printMaxScanNum();
//		aaSet.registerEnzyme(Enzyme.TRYPSIN);
//		System.out.println(aaSet.getNeighboringAACleavageCredit()+" "+aaSet.getNeighboringAACleavagePenalty());
//		System.out.println(aaSet.getPeptideCleavageCredit()+" "+aaSet.getPeptideCleavagePenalty());
//		testSpecKey();
//		String pep = "-.asdf.AA";
//		System.out.println(pep.matches(".\\..+\\.."));
		System.out.println((byte)200);
	}

	public static void printMaxScanNum() throws Exception
	{
		String mzXMLFileName = "/home/sangtaekim/Test/CCMS/3297b97db35241ba8547906b22377869/spectrum/00000.mzXML";
		MzXMLSpectraMap map = new MzXMLSpectraMap(mzXMLFileName);
		System.out.println(map.getMaxScanNumber());
	}
	
	public static void compareFiles() throws Exception
	{
		String mzXMLFileName = "/home/sangtaekim/Test/CCMS/3297b97db35241ba8547906b22377869/spectrum/00000.mzXML";
		String mgfFileName = "/home/sangtaekim/Test/CCMS/3297b97db35241ba8547906b22377869/spectrum/00000.mgf";
		
		SpectraIterator itr = new SpectraIterator(mgfFileName, new MgfSpectrumParser());
		HashSet<Integer> mgfScanNum = new HashSet<Integer>();
		while(itr.hasNext())
		{
			mgfScanNum.add(itr.next().getScanNum());
		}
		
		MzXMLSpectraIterator itr2 = new MzXMLSpectraIterator(mzXMLFileName);
		HashSet<Integer> mzXMLScanNum = new HashSet<Integer>();
		while(itr2.hasNext())
		{
			mzXMLScanNum.add(itr2.next().getScanNum());
		}
		
		int mzXMLOnly = 0;
		for(int scanNum : mzXMLScanNum)
		{
			if(!mgfScanNum.contains(scanNum))
			{
				System.out.println("MzXMLOnly: " + scanNum);
				mzXMLOnly++;
			}
		}

		int mgfOnly = 0;
		for(int scanNum : mgfScanNum)
		{
			if(!mzXMLScanNum.contains(scanNum))
			{
				System.out.println("MgfOnly: " + scanNum);
				mgfOnly++;
			}
		}
		
		System.out.println("MzXMLOnly: " + mzXMLOnly);
		System.out.println("MgfOnly: " + mgfOnly);
		
	}
	public static void efdrTest() throws Exception
	{
		double[] specProb = new double[136964];
		Random rand = new Random();
		for(int i=0; i<specProb.length; i++)
		{
			specProb[i] = rand.nextDouble();
//			Arrays.sort(specProb[i]);
		}
		
		long time = System.currentTimeMillis();
		Arrays.sort(specProb);
		System.out.println("Sorting: " + (System.currentTimeMillis()-time));
		time = System.currentTimeMillis();
		double cumP = 0;
		for(int i=0; i<specProb.length; i++)
		{
			double probCorr = 1.-specProb[i];
			double pValue = MSGFDBResultGenerator.DBMatch.getPValue(specProb[i], 3000000);
			cumP += pValue;
			double eTD = i+1-cumP;
			double eDD = cumP;
			for(int j=1; j<specProb.length; j++)
			{
				eDD += specProb[j];
			}
			double eFDR = eDD/eTD;
			double dummy = eFDR;
		}
		System.out.println("Time: " + (System.currentTimeMillis()-time));
		
	}
	
	public static void combination(int n, int r) throws Exception
	{
		for(LinkedList<Integer> ins : getCombinations(n,r))
			System.out.println(ins);
	}

	public static ArrayList<LinkedList<Integer>> getCombinations(int n, int r) throws Exception
	{
		ArrayList<LinkedList<Integer>> nHr = new ArrayList<LinkedList<Integer>>();
		
		if(r==1)
		{
			for(int i=0; i<n; i++)
			{
				LinkedList<Integer> newIns = new LinkedList<Integer>();
				newIns.add(i);
				nHr.add(newIns);
			}
		}
		else if(r>1)
		{
			ArrayList<LinkedList<Integer>> nHrMinus1 = getCombinations(n, r-1);
			for(LinkedList<Integer> ins : nHrMinus1)
			{
				int prevLargest = ins.getLast();
				for(int i=prevLargest; i<n; i++)
				{
					LinkedList<Integer> newIns = new LinkedList<Integer>(ins);
					newIns.add(i);
					nHr.add(newIns);
				}
			}
		}
		
		return nHr;
	}
	
	public static void mzXMLLoadingTest() throws Exception
	{
		long time = System.currentTimeMillis();
		String fileName = "/home/sangtaekim/Research/Data/HeckWhole/Spectra/090121_NM_Trypsin_20.mzXML";
		MzXMLSpectraMap map = new MzXMLSpectraMap(fileName);
		System.out.println("LoadingTime: " + (System.currentTimeMillis()-time));
		map.getSpectrumBySpecIndex(100);
		System.out.println("Time: " + (System.currentTimeMillis()-time));
	}
	
	public static void extractRECOMBCPEmails() throws Exception
	{
		String fileName = System.getProperty("user.home")+"/Documents/RECOMB-CP/participantsEmails.txt";
		BufferedLineReader in = new BufferedLineReader(fileName);
		String s;
		while((s=in.readLine()) != null)
		{
			if(s.matches(".+@.+\\..+"))
				System.out.println(s.split("\\s+")[0]);
		}
	}
	public static void extractRECOMBBEAbstracts() throws Exception
	{
		String fileName = System.getProperty("user.home")+"/Desktop/abstracts.txt";
		BufferedLineReader in = new BufferedLineReader(fileName);
		String s;
		int lineNum=0;
		Hashtable<String, String> entries = new Hashtable<String,String>();
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("ID"))
				continue;
			lineNum++;
			String[] token = s.split("\t");
			if(token.length < 12)
				continue;
			String email = token[7];
			if(!email.contains("sangtae") && !email.contains("chorta"))
				entries.put(email, s);
		}
		
		ArrayList<String> list = new ArrayList<String>();
		for(String e : entries.values())
			list.add(e);
		Collections.sort(list);
		for(String e : list)
		{
			String[] token = e.split("\t");
			System.out.println("****************************");
			System.out.println("Date Submitted: " + token[1]);
			System.out.println("Category: " + token[2]);
			System.out.println("Title: " + token[3]);
			System.out.println("Author: " + token[4] + " " + token[5]);
			System.out.println("Co-Authors: " + token[6]);
			System.out.println("Email: " + token[7]);
			System.out.println("Affiliation: " + token[8]);
			System.out.println("Topic: " + token[9]);
			System.out.println("Keywords: " + token[10]);
			System.out.println();
			System.out.println("Abstract:");
			for(int i=11; i<token.length; i++)
				System.out.print(token[i]+"\t");
			System.out.println();
			System.out.println();
		}
	}
	
	public static void extractEmails() throws Exception
	{
		String fileName = "/Users/sangtaekim/Documents/RECOMB-CP/asms09.txt";
		BufferedLineReader in = new BufferedLineReader(fileName);
		String s;
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\\s+");
			for(String t : token)
			{
				if(t.contains("@") && t.contains("."))
				{
					System.out.println(t);
				}
			}
		}
	}
}
