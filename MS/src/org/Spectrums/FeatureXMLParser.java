/**
 * Parse output from feature finder of openMS
 */
package org.Spectrums;
import java.io.IOException;
import java.util.ArrayList;
import java.util.GregorianCalendar;
import java.util.List;
import java.util.Map;

import javax.xml.parsers.*;
import org.apache.commons.collections.map.MultiValueMap;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

public class FeatureXMLParser {
	private Document dom;
	private List<MSFeature> featureList;
	
	public Document getDom() {
		return dom;
	}

	public void setDom(Document dom) {
		this.dom = dom;
	}

	public List<MSFeature> getFeatureList() {
		return featureList;
	}

	public void setFeatureList(List<MSFeature> featureList) {
		this.featureList = featureList;
	}

	public FeatureXMLParser(String filename){
		parseXmlFile(filename);
		parseDocument();
	}
	
	private void parseXmlFile(String filename){
		//get the factory
		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
		
		try {

			//Using factory get an instance of document builder
			DocumentBuilder db = dbf.newDocumentBuilder();

			//parse using builder to get DOM representation of the XML file
			dom = db.parse(filename);
			//dom = db.parse("..//mixture_linked//employee.xml");


		}catch(ParserConfigurationException pce) {
			pce.printStackTrace();
		}catch(SAXException se) {
			se.printStackTrace();
		}catch(IOException ioe) {
			ioe.printStackTrace();
		}
	}
	private void parseDocument(){
		//get the root element
		Element docEle = dom.getDocumentElement();
		//System.out.println("encoding is: " + dom.getInputEncoding());
		//get a nodelist of 
		
		NodeList nl = docEle.getElementsByTagName("feature");
		System.out.println("There are total of features: " + nl.getLength());
		this.featureList = new ArrayList();
		if(nl != null && nl.getLength() > 0) {
			int count = 0;
			for(int i = 0 ; i < nl.getLength(); i++) {
				//get the employee element
				Element el = (Element)nl.item(i);
				//System.out.println("feature id is: " + el.getAttribute("id"));			
				Node parent = el.getParentNode();
				if(parent.getNodeName().equals("featureList")){
					MSFeature feature = new MSFeature();
					double[] timeSpan = parseHullPoint(el);
					double[] location = getPositionValue(el, "position");
					feature.setId(el.getAttribute("id"));
					feature.setIntensity(getDoubleValue(el, "intensity"));
					feature.setMz(location[0]);
					feature.setQuality(getDoubleValue(el, "overallquality"));
					feature.setCharge(getIntValue(el, "charge"));
					//feature.setFeatureList(parseFeature(el));
					feature.setMinRT(timeSpan[0]);
					feature.setMaxRT(timeSpan[1]);
					feature.setRt(0.5*(feature.getMaxRT()-feature.getMinRT())+feature.getMinRT());
					feature.setMz(location[1]);
					System.out.print(el.getAttribute("id") + "\t" 
							+ getTextValue(el, "intensity") + "\t"
							+ getTextValue(el, "overallquality") + "\t"
							+ location[0] + "\t" + location[1] + "\t" 
							+ getTextValue(el, "charge") + "\t"
							+ timeSpan[0] +"\t" + timeSpan[1] +"\t"
							+ "\n");
					
					featureList.add(feature);
					count++;
				}
				
			}
			System.out.println("there are total of " + count + " main features");
//			for(int i =0; i < featureList.size(); i++){
//				System.out.println(featureList.get(i));
//				
//			}
		}
	}
	
	private MultiValueMap parseFeature(Element el){
		MultiValueMap map = new MultiValueMap();
		List<Element> hullPoints = getTopLevelHullPoints(el);
		System.out.println("extracted hullpoints: " + hullPoints.size());
		for(int i = 0; i < hullPoints.size(); i++){
				Element hullpoint = hullPoints.get(i);
				double[] positions = parseHullPoint(hullpoint);
				map.put(positions[0], positions[1]);
		}
		return map;
	}
	/**
	 * Only extract top level hull points, ignore those in subordinate 
	 * @param el
	 * @return
	 */
	private List<Element> getTopLevelHullPoints(Element el){
		NodeList nl = el.getElementsByTagName("hullpoint");
		List<Element> hullPoints = new ArrayList();
		for(int i = 0; i < nl.getLength(); i++){
			Element hullpoint = (Element)nl.item(i);
			if(hullpoint.getParentNode().getParentNode()
					.getParentNode().getNodeName().equals("featureList")){
				hullPoints.add(hullpoint);
			}
		}
		return hullPoints;
	}
	
	//parse the hullpoint object, extract begin and end elution time of this peptide feature
	private double[] parseHullPoint(Element hullpoint){
		double[] positions = new double[2];
		double beginTime = 0;
		double endTime = 0;
		NodeList nl = hullpoint.getElementsByTagName("convexhull");
		//System.out.println("convexhull size: " + nl.getLength());
		for(int i = 0; i < nl.getLength(); i++){
			Element hullLine = (Element)nl.item(i);
			NodeList pts = hullLine.getElementsByTagName("pt");
			//System.out.println("number of points in hull " + hullLine.getElementsByTagName("pt").getLength());
			for(int j = 0; j < pts.getLength(); j++){
				//Element pt = (Element)pts.item(j);
				//System.out.println("point " + pt.getAttribute("x") + "\t" + pt.getAttribute("y"));
			}
			if(i==0){ //want to get widest time-span, presumably it is the first window width
				beginTime = Double.parseDouble(((Element)pts.item(0)).getAttribute("x"));
				endTime = Double.parseDouble(((Element)pts.item(1)).getAttribute("x"));
			}
		}
		positions[0] = beginTime;
		positions[1] = endTime;
		return positions;
	}
	
	private double[] parsePosition(Element feature){
		double[] positions = new double[2];
		NodeList nl = feature.getElementsByTagName("position");
		positions[0] = Double.parseDouble(nl.item(0).getFirstChild().getNodeValue());
		positions[1] = Double.parseDouble(nl.item(1).getFirstChild().getNodeValue());
		System.out.println("position is: " + positions[0] + "\t" + positions[1]);
		return positions;
	}
	
	private String getTextValue(Element ele, String tagName) {
		String textVal = null;
		NodeList nl = ele.getElementsByTagName(tagName);
		if(nl != null && nl.getLength() > 0) {
			Element el = (Element)nl.item(0);
			textVal = el.getFirstChild().getNodeValue();
		}

		return textVal;
	}
	
	private double[] getPositionValue(Element ele, String tagName) {
		String textVal = null, textVal2 = null;
		NodeList nl = ele.getElementsByTagName(tagName);
		if(nl != null && nl.getLength() > 0) {
			Element el = (Element)nl.item(0);
			textVal = el.getFirstChild().getNodeValue();
			Element el2 = (Element)nl.item(1);
			textVal2 = el2.getFirstChild().getNodeValue();
		}

		return new double[]{Double.parseDouble(textVal), 
				Double.parseDouble(textVal2)};
	}
	
	private int getIntValue(Element ele, String tagName) {
		return Integer.parseInt(getTextValue(ele,tagName));
	}

	private double getDoubleValue(Element ele, String tagName) {
		//in production application you would catch the exception
		return Double.parseDouble(getTextValue(ele,tagName));
	}
	
	public void getFeaturePair(double expectedMassDiff, double massTolerance, double maxRTDiff, int mode){
		for(int i = 0; i < this.featureList.size(); i++){
			MSFeature feature1 = featureList.get(i);
			for(int j = i+1; j < this.featureList.size(); j++){
				MSFeature feature2 = featureList.get(j);
				if(feature1.getCharge() == feature2.getCharge()
						//&& (feature1.getMz() - feature2.getMz()) * (feature1.getRt() - feature2.getRt()) > 0   //we assume the lighter form should elude first, not neccessarily this requirement 
						&& SWATHUtils.checkMass(Math.abs(feature1.getMz() - feature2.getMz()), expectedMassDiff/feature1.getCharge(), massTolerance, mode)
						&& Math.abs(feature1.getMinRT() - feature2.getMinRT()) < maxRTDiff){
					System.out.println("Isotopic pairs:\t" + feature1 + "\t" + feature2);
				}
			}
		}
	}
	
	public void mapScanFromRT(String spectrumLibFile){
		MZXMLReader reader = new MZXMLReader(spectrumLibFile);
		List<Spectrum> MS2 = reader.readAllMS2Spectra();
		Map<Integer, Double> rtMap = reader.getRTScanMapping();
		Map<Double, Integer> rtMap2 = reader.getRTScanMappingReverse();
		for(int j = 0; j < this.featureList.size(); j++){
			MSFeature feature = this.featureList.get(j);
			feature.setScan(rtMap2.get(feature.getRt()));
			System.out.println(this.featureList.get(j) );					
		}
	}
	
	public static void testFeatureParser(String featureFile){
		FeatureXMLParser parser = new FeatureXMLParser(featureFile);
	}	
	
	
	public static void testMapMS2ToFeature(){
		String spectrumLibFile = "..\\mixture_linked\\msdata\\UPS_Ecoli\\40fmol_UPS1_IDA_2.mzXML";
		String searchResult = "..\\mixture_linked\\ACG_swathdevelopment_40fmolUPS1_1pepFDR_allPSMs_msgfdb.txt";
		FeatureXMLParser parser = new FeatureXMLParser("..//mixture_linked/openMS/feature_detection/TOPPAS_out/FeatureFinderIDA/40fmol_UPS1_IDA_2.featureXML");
		List<String> results =Utils.FileIOUtils.createListFromFile(searchResult);
		List<MSFeature> featureList = parser.getFeatureList();
		MZXMLReader reader = new MZXMLReader(spectrumLibFile);
		Map<Integer, Double> rtMap = reader.getRTScanMapping();
		int scanNumInd = 1;
		int pepInd = 7;
		int protInd = 8;
		//parser.mapScanFromRT(spectrumLibFile);
		System.out.println("finish mapping");
		for(int i = 0; i < results.size(); i++){
			String result = results.get(i);
			String[] tokens = result.split("\\t");
			int scan = Integer.parseInt(tokens[scanNumInd]);
			Spectrum s = reader.getSpectrum(scan);
			s.peptide = tokens[pepInd];
			double RT = rtMap.get(s.scanNumber);
			for(int j = 0; j < featureList.size(); j++){
				MSFeature feature = featureList.get(j);
				if(feature.isWithnFeature(s.parentMass, s.charge, RT, 0.03)){
					System.out.println("Spectrum " + s.scanNumber +  "\t" + s.parentMass + "\t" + s.charge  +"\t" + RT  + "\t" + tokens[pepInd] + "\t" + tokens[protInd] +"\t"
							+ " is within feature " + feature.getId() +"\t" + feature.getRt() + "\t" + feature.getMz() +"\t" + feature.getCharge() + feature.getIntensity());
					break;
				}
				if(j == featureList.size()-1){
					System.out.println("Spectrum " + s.scanNumber +  "\t" + s.parentMass + "\t" + s.charge  + "\t" + RT + "\t" + tokens[pepInd] + "\t" + tokens[protInd] +"\t"
							+ " is not within any feature ");
				}
			}
		}
		
	}
	
	public static void testMapMSFeatureToPeptide(){
		CombinatoryPeptides combPeps = new CombinatoryPeptides("[TN][PG][AY]K[EF][IQ][DS]R");
		String spectrumLibFile = "..\\mixture_linked\\linked_peptide_library\\Crosslinked_peps_lib-2_cen.mzXML";
		//CombinatoryPeptides combPeps = new CombinatoryPeptides("[NP][ET][GY]K[SQ][IF][AD]R");
		List<String> peptides = combPeps.generateAllPeptides();
		//peptides = addMod(peptides, 42.010564686);
		peptides = combPeps.cutPeptides(peptides);
		System.out.println("we have peptides: " + peptides.size());
		//peptides = addMod(peptides);
		CandidateSpectrumLibFactory factory = 	
			CandidateSpectrumLibFactory.createFactoryFromPeptide(peptides);
		factory.setMinCharge(1);
		factory.setMaxCharge(1);
		factory.indexPeptideByParentMass(0.1);
		//factory.reIndexPeptideByParentMass(0.1);
		//factory.insertPTM(Mass.DSSDANGLE_MASS, new char[]{'K'}, 1);
		//factory.insertPTM(Mass.DSSDANGLE_MASS_D12, new char[]{'K'}, 1);
		//factory.insertPTM(0.9847, new char[]{'N', 'Q'}, 2);
		factory.crossLinkAllPeptides(2, 4);
		FeatureXMLParser parser = new FeatureXMLParser("..//mixture_linked//test.xml");
		parser.mapScanFromRT(spectrumLibFile);
		List<Peptide> pepList = factory.getAllPeptide();
		List<MSFeature> featureList = parser.getFeatureList();
		for(int i = 0; i < pepList.size(); i++){
			Peptide current = pepList.get(i);
			Peptide currentD12 = new Peptide(current);
			currentD12.setParentmass(currentD12.getParentmass()+12.0759/currentD12.getCharge());
			for(int j = 0; j < featureList.size(); j++){
				MSFeature feature = featureList.get(j);
				if(feature.isMatchFeature(current, 0.03)){
					System.out.println("Peptide: " + current + "\t" + current.getParentmass() + "\t" + current.getCharge() 
							+ " match to feature " + feature);
				}
				if(feature.isMatchFeature(currentD12, 0.03)){
					System.out.println("Peptide-Heavy: " + current + "\t" + current.getParentmass() + "\t" + current.getCharge() 
							+ " match to feature " + feature);
				}
			}
		}
	}
	
	public static void testMapPeptideToMSFeature(){
		CombinatoryPeptides combPeps = new CombinatoryPeptides("[TN][PG][AY]K[EF][IQ][DS]R");
		List<String> peptides = combPeps.generateAllPeptides();
		//peptides = addMod(peptides, 42.010564686);
		peptides = combPeps.cutPeptides(peptides);
		System.out.println("we have peptides: " + peptides.size());
		CandidateSpectrumLibFactory factory = 	
			CandidateSpectrumLibFactory.createFactoryFromPeptide(peptides);
		factory.setMinCharge(1);
		factory.setMaxCharge(1);
		factory.indexPeptideByParentMass(0.1);
		factory.crossLinkAllPeptides(2, 4);
		//String spectrumLibFile = "..\\mixture_linked\\linked_peptide_library\\20100215_nLCFTMSMS_CLPL_2pmol_1841_MSonly.mzXML";
		String spectrumLibFile = "..\\mixture_linked\\linked_peptide_library\\Crosslinked_peps_lib-2_cen.mzXML";
		MZXMLReader reader = new MZXMLReader(spectrumLibFile);
		PrecursorMassChecker checker = new PrecursorMassChecker(spectrumLibFile);
		long start = (new GregorianCalendar()).getTimeInMillis();
		FeatureXMLParser parser = new FeatureXMLParser("..//mixture_linked//test.xml");
		List<MSFeature> featureList = parser.getFeatureList();
		parser.mapScanFromRT(spectrumLibFile);
		Map<Integer, Double> rtMap = reader.getRTScanMapping();
		List<Peptide> pepList = factory.getAllPeptide();
		for(int i = 0; i < pepList.size(); i++){
			Peptide current = pepList.get(i);
			Peptide currentD12 = new Peptide(current);
			currentD12.setParentmass(currentD12.getParentmass()+12.0759/currentD12.getCharge());
			//System.out.println("checking : " + current + "\t" + current.getCharge());
			//int matched = checker.matchPeptidePrecursorProfile2(current, 10);
			List[] matched = checker.matchPeptidePrecursorProfilePairDetail(current, 70, (Mass.DEUTERIUM-Mass.PROTON_MASS)*12);
			for(int k = 0; k < matched[2].size(); k++){
				System.out.println("checking1 : " + matched[2].get(k) + "\t" + current + "\t" + current.getParentmass());
				for(int j = 0; j < featureList.size(); j++){
					MSFeature feature = featureList.get(j);
					double RT = rtMap.get(matched[2].get(k));
					if(feature.isMatchFeature(current, 0.03, RT)){
						System.out.println(current +  "\t" + current.getParentmass() + "\t" + current.getCharge()
								+ " with isotopic labels match to feature " + feature);
						break;
					}
					if(feature.isMatchFeature(currentD12, 0.03, RT)){
						System.out.println(current +  "\t" + currentD12.getParentmass() + "\t" + currentD12.getCharge()
								+ " with isotopic labels match to Heavy-feature " + feature);
						break;
					}
				}				
			}
		}		
		System.out.println("matching spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
	public static void testFindFeaturePair(){
		String featureFile = "..//mixture_linked//openMS//feature_detection\\TOPPAS_out\\005-FeatureFinderCentroided\\aleitner_M1108_136.featureXML";
		FeatureXMLParser parser = new FeatureXMLParser(featureFile);
		parser.getFeaturePair(12.0759, 50, 60, 2);
	}
	
	
	
	public static void main(String[] args){
		//String featureFile = "..\\mixture_linked\\SWATH\\FeaturesDetect\\400fmolUPS1_1ugHumanLysate_IDA.featureXMl";
		String featureFile = args[0];
		//testFeatureParser(featureFile);
		//testMapMS2ToFeature();
		//testMapMSFeatureToPeptide();
		//testMapPeptideToMSFeature();
		testFindFeaturePair();
	}
}

