package misc;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.util.Calendar;

import msscorer.NewRankScorer;
import msscorer.ScoringParameterGeneratorWithErrors;
import msutil.ActivationMethod;
import msutil.AminoAcidSet;
import msutil.Enzyme;
import msutil.InstrumentType;
import msutil.Modification;

public class TrainScoringParameters {
	
	private static final String PARAM_DIR = System.getProperty("user.home")+"/Developments/MS_Java_Dev/src/resources/ionstat";
	private static final String BACKUP_DIR = System.getProperty("user.home")+"/Developments/ionstat";
	private static final String SPEC_DIR = System.getProperty("user.home")+"/Research/Data/IonStat/SpectraForTraining";
	
	public static void main(String argv[]) throws Exception
	{
//		backup();
		createParamFiles();
		testParamFiles();
	}
	
	public static void backup() throws Exception
	{
		File paramDir = new File(PARAM_DIR);
		boolean paramExists = false;
		for(File paramFile : paramDir.listFiles())
		{
			if(paramFile.getName().endsWith(".param"))
				paramExists = true;
		}
		if(!paramExists)
		{
			System.out.println("No param file to backup.");
			return;
		}
		Calendar calendar = Calendar.getInstance();
		String dateStr = calendar.get(Calendar.MONTH) + "_" + calendar.get(Calendar.DAY_OF_MONTH) + "_" + calendar.get(Calendar.YEAR);
		String backupDirName = "ParamBackup_" + dateStr;
		File backupDir = new File(BACKUP_DIR+"/"+backupDirName);
		if(backupDir.exists())
		{
			System.out.println("Backup directory already exists: " + backupDir.getPath());
			System.exit(-1);
		}
		backupDir.mkdir();
		System.out.println(backupDir.getPath()+ " is created.");
		
		boolean backupSuccess = true;
		for(File paramFile : paramDir.listFiles())
		{
			if(paramFile.getName().endsWith(".param"))
			{
				File newFile = new File(backupDir, paramFile.getName());
				boolean isBackupSuccessful = paramFile.renameTo(newFile);
				System.out.println("Moving " + paramFile.getPath() + " to " + newFile.getPath() + (isBackupSuccessful ? " succeeded." : " failed."));
				if(!isBackupSuccessful)
				{
					backupSuccess = false;
					break;
				}
			}
		}
		if(backupSuccess)
			System.out.println("Backup complete.");
		else
		{
			backupDir.delete();
			System.out.println(backupDir.getPath()+ " is deleted.");
			System.out.println("Backup failed.");
			System.exit(0);
		}
	}
	
	public static void createParamFiles() throws Exception
	{
		File specDir = new File(SPEC_DIR);
		if(!specDir.exists())
		{
			System.err.println("Training spectra directory doesn't exist:" + specDir.getPath());
			System.exit(-1);
		}
		
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		for(File specFile : specDir.listFiles())
		{
			String specFileName = specFile.getName();
			if(specFileName.endsWith(".mgf"))
			{
				String id = specFileName.substring(0, specFileName.lastIndexOf('.'));
				String[] token = id.split("_");
				if(token.length != 3 && token.length != 4)
				{
					System.err.println("Wrong file name: " + specFile.getName());
					System.exit(-1);
				}
				String actMethodStr = token[0];
				String instTypeStr = token[1];
				String enzymeStr = token[2];
				String modStr = null;
				if(token.length == 4)
					modStr = token[3];
				
				ActivationMethod actMethod = ActivationMethod.get(actMethodStr);
				InstrumentType instType = InstrumentType.get(instTypeStr);
				Enzyme enzyme = Enzyme.getEnzymeByName(enzymeStr);
				Modification mod = null;
				if(modStr != null)
					mod = Modification.get(modStr);
				
				if(actMethod == null || instType == null || enzyme == null)
				{
					System.err.println("Wrong file name: " + specFile.getName());
					System.exit(-1);
				}
				
				String paramFileName = id+".param";
				File outputFile = new File(PARAM_DIR, paramFileName);
				System.out.println("Generating " + outputFile.getPath());
				int numSpecsPerPeptide = 1;
				int errorScalingFactor = 0;
				boolean deconvoluteSpectra = false;
				if(instType == InstrumentType.HIGH_RESOLUTION_LTQ)
				{
					numSpecsPerPeptide = 3;
					errorScalingFactor = 100;
					deconvoluteSpectra = true;
				}
				else if(instType == InstrumentType.TOF)
				{
					errorScalingFactor = 100;
					deconvoluteSpectra = true;
				}
				
				boolean considerPhosLoss = false;
				if(mod == Modification.get("Phosphorylation"))
				{
					numSpecsPerPeptide = 3;
					considerPhosLoss = true;
				}
				if(enzyme == Enzyme.ALP)
					numSpecsPerPeptide = 3;
				
				if(actMethod == ActivationMethod.HCD)
				{
					deconvoluteSpectra = true;
					errorScalingFactor = 0;
					ScoringParameterGeneratorWithErrors.generateParameters(specFile, actMethod, instType, enzyme, numSpecsPerPeptide, errorScalingFactor, considerPhosLoss, deconvoluteSpectra, outputFile, aaSet, false, false);
				}
			}
		}
		System.out.println("Successfully generated parameters!");
	}
	
	public static void testParamFiles() throws Exception
	{
		for(File f : new File(PARAM_DIR).listFiles())
		{
			if(f.getName().endsWith(".param"))
			{
				System.out.println("Reading " + f.getName());
				InputStream is = new BufferedInputStream(new FileInputStream(f));
				NewRankScorer scorer = new NewRankScorer(new BufferedInputStream(is));
				System.out.println(scorer.getActivationMethod().getName()+"_"+scorer.getInstrumentType().getName()+"_"+scorer.getEnzyme().getName());
			}
		}
		System.out.println("Read Success");
	}
}
