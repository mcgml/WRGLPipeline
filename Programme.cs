﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace WRGLPipeline
{
    class Programme
    {
        public const double WRGLversion = 1.7;

        private static void Main(string[] args)
        {
            if (args.Length != 1) {
                throw new ArgumentNullException(@"ERROR: Enter path to MiSeq alignment folder");
            }

            //get run parameters
            string suppliedDir = args[0];
            string localFastqDir = AuxillaryFunctions.GetFastqDir(suppliedDir);
            string localRootRunDir = AuxillaryFunctions.GetRootRunDir(suppliedDir);
            string localMiSeqAnalysisDir = AuxillaryFunctions.GetLocalAnalysisFolderDir(suppliedDir);
            string sampleSheetPath = suppliedDir + @"\SampleSheetUsed.csv";
            string runID = AuxillaryFunctions.GetRunID(suppliedDir);
            string localLogFilename = localRootRunDir + @"\WRGLPipeline.log";
            string networkRootRunDir = AuxillaryFunctions.makeNetworkOutputDir(@"\\GENSTORE-ED\Genetics\Illumina\MiSeqOutput\" + runID);
            ProgrammeParameters parameters = new ProgrammeParameters(); //read ini file

            //Parse samplesheet
            ParseSampleSheet sampleSheet = new ParseSampleSheet(sampleSheetPath, localFastqDir, localLogFilename, parameters);

            //write variables to log
            AuxillaryFunctions.WriteLog(@"Run identifier: " + runID, localLogFilename, 0, true, parameters);
            AuxillaryFunctions.WriteLog(@"Local root run directory: " + localRootRunDir, localLogFilename, 0, false, parameters);
            AuxillaryFunctions.WriteLog(@"Local FASTQ directory: " + localFastqDir, localLogFilename, 0, false, parameters);
            AuxillaryFunctions.WriteLog(@"Supplied directory: " + suppliedDir, localLogFilename, 0, false, parameters);
            AuxillaryFunctions.WriteLog(@"Local MiSeq analysis directory: " + localMiSeqAnalysisDir, localLogFilename, 0, false, parameters);
            AuxillaryFunctions.WriteLog(@"Local output directory: " + localRootRunDir, localLogFilename, 0, false, parameters);
            AuxillaryFunctions.WriteLog(@"Network output directory: " + networkRootRunDir, localLogFilename, 0, false, parameters);
            AuxillaryFunctions.WriteLog(@"Local SampleSheet path: " + sampleSheetPath, localLogFilename, 0, false, parameters);
            AuxillaryFunctions.WriteLog(@"Experiment name: " + sampleSheet.getExperimentName, localLogFilename, 0, false, parameters);
            AuxillaryFunctions.WriteLog(@"Investigator name: " + sampleSheet.getInvestigatorName, localLogFilename, 0, false, parameters);

            //compute MD5 for fastqs; copy fastqs, metrics, samplesheet and MD5 to network
            FileManagement.BackupFiles(localFastqDir, suppliedDir, localRootRunDir, networkRootRunDir, localLogFilename, parameters);

            if (sampleSheet.getAnalyses.Count > 0)
            {
                Dictionary<string, Tuple<string, string>> fastqFileNames = new Dictionary<string, Tuple<string, string>>(GetFASTQFileNames(sampleSheet, localFastqDir, localLogFilename, parameters));

                //AnalyseSamples
                if (sampleSheet.getAnalyses.ContainsKey("G"))
                {
                    GenotypingPipelineWrapper genotypingAnalysis = new GenotypingPipelineWrapper(sampleSheet, localFastqDir, localLogFilename, runID, parameters, networkRootRunDir, fastqFileNames);
                }

                if (sampleSheet.getAnalyses.ContainsKey("P"))
                {
                    PanelPipeline panelAnalysis = new PanelPipeline(sampleSheet, localFastqDir, runID, localLogFilename, parameters, networkRootRunDir, fastqFileNames);
                }
            }

            //delete local oldest run
            if (parameters.getDeleteOldestLocalRun == true) { FileManagement.DeleteOldRuns(localLogFilename, localMiSeqAnalysisDir, parameters); }

            //copy local files to network
            File.Copy(localLogFilename, networkRootRunDir + @"\" + Path.GetFileName(localLogFilename)); //logfile
        }

        private static Dictionary<string, Tuple<string, string>> GetFASTQFileNames(ParseSampleSheet sampleSheet, string localFastqDir, string logFilename, ProgrammeParameters parameters)
        {
            Dictionary<string, Tuple<string, string>> fastqFileNames = new Dictionary<string, Tuple<string, string>>();
            string[] read1Files, read2Files;

            foreach (SampleRecord record in sampleSheet.getSampleRecords)
            {
                //get FASTQ filenames
                if (record.Sample_Name == "")
                {
                    read1Files = Directory.GetFiles(localFastqDir, record.Sample_ID + @"_*_R1_*.fastq.gz");
                    read2Files = Directory.GetFiles(localFastqDir, record.Sample_ID + @"_*_R2_*.fastq.gz");
                }
                else
                {
                    read1Files = Directory.GetFiles(localFastqDir, record.Sample_Name + @"_*_R1_*.fastq.gz");
                    read2Files = Directory.GetFiles(localFastqDir, record.Sample_Name + @"_*_R2_*.fastq.gz");
                }

                if (read1Files.Length == 0 || read2Files.Length == 0) //no paired-end FASTQs found for this Sample_ID
                {
                    AuxillaryFunctions.WriteLog(@"Paired-end FASTQ file(s) not found for " + record.Sample_ID, logFilename, -1, false, parameters);
                    throw new FileNotFoundException();
                }
                else
                {
                    Tuple<string, string> fileNames = new Tuple<string, string>(read1Files[0], read2Files[0]);
                    fastqFileNames.Add(record.Sample_ID, fileNames);
                }
            }

            return fastqFileNames;
        }
    }
}

