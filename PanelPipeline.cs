﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using WinSCP;
using System.IO;
using System.Threading;
using System.Text.RegularExpressions;

namespace WRGLPipeline
{
    class PanelPipeline
    {
        const double PanelPipelineVerison = 2.0;

        //tunnel connection settings
        private string scratchDir;
        private string localFastqDir, localAnalysisDir, networkAnalysisDir, runID, logFilename, localReportFilename, networkReportFilename, winscpLogPath, networkRootRunDir;
        ProgrammeParameters parameters;
        ParseSampleSheet sampleSheet;
        Dictionary<string, HashSet<string>> failedAmplicons = new Dictionary<string, HashSet<string>>();
        ParseBED coreBEDRecords;
        ParseBED targetBEDRecords;
        Dictionary<string, Tuple<string, string>> fastqFileNames;

        private SessionOptions iridis4a = new SessionOptions
        {
            Protocol = Protocol.Sftp,
            HostName = @"iridis4_a.soton.ac.uk",
        };
        
        private SessionOptions iridis4b = new SessionOptions
        {
            Protocol = Protocol.Sftp,
            HostName = @"iridis4_b.soton.ac.uk",
        };

        private SessionOptions iridis4c = new SessionOptions
        {
            Protocol = Protocol.Sftp,
            HostName = @"iridis4_c.soton.ac.uk",
        };

        public PanelPipeline(ParseSampleSheet _sampleSheet, string _localFastqDir, string _runID, string _logFilename, ProgrammeParameters _parameters, 
            string _networkRootRunDir, Dictionary<string, Tuple<string, string>> _fastqFileNames)
        {
            this.localFastqDir = _localFastqDir;
            this.runID = _runID;
            this.logFilename = _logFilename;
            this.sampleSheet = _sampleSheet;
            this.localAnalysisDir = localFastqDir + @"\Panel_" + PanelPipelineVerison;
            this.networkAnalysisDir = _networkRootRunDir + @"\Panel_" + PanelPipelineVerison;
            this.localReportFilename = localAnalysisDir + @"\" + runID + @"_Panel_" + PanelPipelineVerison + ".report";
            this.networkReportFilename = networkAnalysisDir + @"\" + runID + @"_Panel_" + PanelPipelineVerison + ".report";
            this.winscpLogPath = localAnalysisDir + @"\" + runID + @"_WinSCP_Transfer.log";
            this.parameters = _parameters;
            this.scratchDir = @"/scratch/WRGL/";
            this.networkRootRunDir = _networkRootRunDir;
            this.coreBEDRecords = new ParseBED(parameters.getCoreBedFile, logFilename, parameters);
            this.targetBEDRecords = new ParseBED(sampleSheet.getAnalyses[@"P"], logFilename, parameters);
            this.fastqFileNames = _fastqFileNames;

            ExecutePanelPipeline();
        }

        private void ExecutePanelPipeline()
        {
            AuxillaryFunctions.WriteLog(@"Starting panel pipeline...", logFilename, 0, false, parameters);

            //configure connection
            iridis4a.UserName = parameters.getSotonUserName;
            iridis4a.Password = ProgrammeParameters.ToInsecureString(parameters.getSotonPassWord);
            iridis4a.SshHostKeyFingerprint = parameters.getIridisHostKey;
            iridis4a.AddRawSettings(@"Tunnel", @"1");
            iridis4a.AddRawSettings(@"TunnelHostName", @"ssh.soton.ac.uk");
            iridis4a.AddRawSettings(@"TunnelUserName", parameters.getSotonUserName);
            iridis4a.AddRawSettings(@"TunnelPasswordPlain", ProgrammeParameters.ToInsecureString(parameters.getSotonPassWord));
            iridis4a.AddRawSettings(@"TunnelHostKey", parameters.getSSHHostKey);

            iridis4b.UserName = parameters.getSotonUserName;
            iridis4b.Password = ProgrammeParameters.ToInsecureString(parameters.getSotonPassWord);
            iridis4b.SshHostKeyFingerprint = parameters.getIridisHostKey;
            iridis4b.AddRawSettings(@"Tunnel", @"1");
            iridis4b.AddRawSettings(@"TunnelHostName", @"ssh.soton.ac.uk");
            iridis4b.AddRawSettings(@"TunnelUserName", parameters.getSotonUserName);
            iridis4b.AddRawSettings(@"TunnelPasswordPlain", ProgrammeParameters.ToInsecureString(parameters.getSotonPassWord));
            iridis4b.AddRawSettings(@"TunnelHostKey", parameters.getSSHHostKey);

            iridis4c.UserName = parameters.getSotonUserName;
            iridis4c.Password = ProgrammeParameters.ToInsecureString(parameters.getSotonPassWord);
            iridis4c.SshHostKeyFingerprint = parameters.getIridisHostKey;
            iridis4c.AddRawSettings(@"Tunnel", @"1");
            iridis4c.AddRawSettings(@"TunnelHostName", @"ssh.soton.ac.uk");
            iridis4c.AddRawSettings(@"TunnelUserName", parameters.getSotonUserName);
            iridis4c.AddRawSettings(@"TunnelPasswordPlain", ProgrammeParameters.ToInsecureString(parameters.getSotonPassWord));
            iridis4c.AddRawSettings(@"TunnelHostKey", parameters.getSSHHostKey);

            //create local output analysis directory
            try { Directory.CreateDirectory(localAnalysisDir); } catch (Exception e) {
                AuxillaryFunctions.WriteLog(@"Could not create local ouput directory: " + e.ToString(), logFilename, -1, false, parameters);
                throw;
            }

            //create network output analysis directory
            try { Directory.CreateDirectory(networkAnalysisDir); } catch (Exception e) {
                AuxillaryFunctions.WriteLog(@"Could not create network ouput directory: " + e.ToString(), logFilename, -1, false, parameters);
                throw;
            }

            //write variables file
            WriteVariablesFile();

            //upload and execute pipeline
            UploadAndExecute();

            //wait before checking download
            AuxillaryFunctions.WriteLog(@"Pipeline idle. Going to sleep...", logFilename, 0, false, parameters);
            Thread.Sleep(1000 * 60 * 180); //ms

            //poll IRIDIS4 for run completion file
            for (int k = 0; k < 15; ++k)
            {
                AuxillaryFunctions.WriteLog(@"Download data attempt " + (k + 1), logFilename, 0, false, parameters);

                if (GetData() == false) //run pending
                {
                    AuxillaryFunctions.WriteLog(@"Pipeline idle. Going to sleep...", logFilename, 0, false, parameters);
                    Thread.Sleep(1000 * 60 * 10); //ms wait 10 mins before checking again
                }
                else
                {
                    AuxillaryFunctions.WriteLog(@"Files downloaded sucessfully", logFilename, 0, false, parameters);

                    WritePanelReport();

                    //copy file to network
                    File.Copy(localReportFilename, networkReportFilename); //report
                    File.Copy(localReportFilename, parameters.getPanelRepo + @"\" + runID + @"_Panel_" + PanelPipelineVerison + ".report");

                    AuxillaryFunctions.WriteLog(@"Variant report path is " + localReportFilename, logFilename, 0, false, parameters);
                    AuxillaryFunctions.SendRunCompletionEmail(logFilename, parameters.getPanelRepo + @"\" + Path.GetFileName(localReportFilename), sampleSheet, @"Panel_" + PanelPipelineVerison, runID, parameters);
                    return;
                }
            }
            
            //data not downloaded
            AuxillaryFunctions.WriteLog(@"Data Colletion Timeout.", logFilename, -1, false, parameters);
            throw new TimeoutException();
        }

        private void UploadAndExecute()
        {
            using (Session session = new Session())
            {
                TransferOperationResult transferResult;
                TransferOptions transferOptions = new TransferOptions();
                transferOptions.TransferMode = TransferMode.Binary;
                StringBuilder bashCommand = new StringBuilder();
                string RemoteSampleFolder;

                //set up logging
                session.SessionLogPath = winscpLogPath;

                try //4a
                {
                    //Connect to iridis4a
                    AuxillaryFunctions.WriteLog(@"Connecting To Iridis4a...", logFilename, 0, false, parameters);
                    session.Open(iridis4a);
                }
                catch (Exception a)
                {
                    AuxillaryFunctions.WriteLog(@"Could not connect: " + a.ToString(), logFilename, -1, false, parameters);

                    try //4b
                    {
                        //Connect to iridis4b
                        AuxillaryFunctions.WriteLog(@"Connecting To Iridis4b...", logFilename, 0, false, parameters);
                        session.Open(iridis4b);
                    }
                    catch (Exception b)
                    {
                        AuxillaryFunctions.WriteLog(@"Could not connect: " + b.ToString(), logFilename, -1, false, parameters);

                        try //4c
                        {
                            //Connect to iridis4c
                            AuxillaryFunctions.WriteLog(@"Connecting To Iridis4c...", logFilename, 0, false, parameters);
                            session.Open(iridis4c);
                        }
                        catch (Exception c)
                        {
                            AuxillaryFunctions.WriteLog(@"Could not connect: " + c.ToString(), logFilename, -1, false, parameters);
                            throw;
                        }
                    }
                }

                //make remote project directory
                try
                {
                    AuxillaryFunctions.WriteLog(@"Creating remote directory " + scratchDir + runID, logFilename, 0, false, parameters);
                    session.CreateDirectory(scratchDir + runID);
                }
                catch (Exception ex)
                {
                    AuxillaryFunctions.WriteLog(@"Could not create remote directory: " + ex.ToString(), logFilename, -1, false, parameters);
                    throw;
                }

                //upload preferred transcripts file
                transferResult = session.PutFiles(parameters.getPreferredTranscriptsFile, scratchDir + runID + @"/", false, transferOptions);
                transferResult.Check(); // Throw on any error

                //loop over Sample_IDs and upload FASTQs
                foreach (SampleRecord record in sampleSheet.getSampleRecords)
                {
                    if (record.Analysis != @"P"){
                        continue;
                    }

                    //output to user
                    AuxillaryFunctions.WriteLog(@"Uploading data for " + record.Sample_ID, logFilename, 0, false, parameters);

                    RemoteSampleFolder = scratchDir + runID + @"/" + record.Sample_ID;

                    //make remote folder for Sample
                    session.CreateDirectory(RemoteSampleFolder);

                    //upload R1 FASTQ
                    transferResult = session.PutFiles(fastqFileNames[record.Sample_ID].Item1, RemoteSampleFolder + @"/", false, transferOptions);
                    transferResult.Check(); // Throw on any error

                    //upload R2 FASTQ
                    transferResult = session.PutFiles(fastqFileNames[record.Sample_ID].Item2, RemoteSampleFolder + @"/", false, transferOptions);
                    transferResult.Check(); // Throw on any error

                    //upload MD5 files
                    transferResult = session.PutFiles(fastqFileNames[record.Sample_ID].Item1 + @".md5", RemoteSampleFolder + @"/", false, transferOptions);
                    transferResult.Check(); // Throw on any error

                    //upload MD5 files
                    transferResult = session.PutFiles(fastqFileNames[record.Sample_ID].Item2 + @".md5", RemoteSampleFolder + @"/", false, transferOptions);
                    transferResult.Check(); // Throw on any error

                    //upload variables file
                    transferResult = session.PutFiles(localAnalysisDir + @"\" + record.Sample_ID + @".variables", RemoteSampleFolder + @"/", false, transferOptions);
                    transferResult.Check(); // Throw on any error

                    //copy WRGL pipeline scripts to RemoteSamplefolder
                    transferResult = session.PutFiles(parameters.getPanelScriptsDir + @"\*.sh", RemoteSampleFolder + @"/", false, transferOptions);
                    transferResult.Check(); // Throw on any error

                    //copy BEDfile to RemoteSamplefolder
                    transferResult = session.PutFiles(sampleSheet.getAnalyses[@"P"], RemoteSampleFolder + @"/", false, transferOptions);
                    transferResult.Check(); // Throw on any error

                    //build BASH command
                    bashCommand.Append(@"cd ");
                    bashCommand.Append(RemoteSampleFolder);
                    bashCommand.Append(@" ");
                    bashCommand.Append(@"&& ");
                    bashCommand.Append(@"qsub 1_Nextera_ReadAlignment.sh");

                    //execute pipeline
                    session.ExecuteCommand(bashCommand.ToString());
                    bashCommand.Clear();
                }
            }
        }

        private bool GetData()
        {
            using (Session session = new Session())
            {
                //set up logging
                session.SessionLogPath = winscpLogPath;

                try //4a
                {
                    //Connect to iridis4a
                    AuxillaryFunctions.WriteLog(@"Connecting to Iridis4a...", logFilename, 0, false, parameters);
                    session.Open(iridis4a);
                }
                catch (Exception a)
                {
                    AuxillaryFunctions.WriteLog(@"Could not connect " + a.ToString(), logFilename, -1, false, parameters);

                    try //4b
                    {
                        //Connect to iridis4b
                        AuxillaryFunctions.WriteLog(@"Connecting to Iridis4b...", logFilename, 0, false, parameters);
                        session.Open(iridis4b);
                    }
                    catch (Exception b)
                    {
                        AuxillaryFunctions.WriteLog(@"Could not connect " + b.ToString(), logFilename, -1, false, parameters);

                        try //4c
                        {
                            //Connect to iridis4c
                            AuxillaryFunctions.WriteLog(@"Connecting to Iridis4c...", logFilename, 0, false, parameters);
                            session.Open(iridis4c);
                        }
                        catch (Exception c)
                        {
                            AuxillaryFunctions.WriteLog(@"Could not connect " + c.ToString(), logFilename, -1, false, parameters);
                            throw;
                        }
                    }
                }
                
                //check if job is complete
                if (session.FileExists(scratchDir + runID + @"/complete"))
                {
                    AuxillaryFunctions.WriteLog(@"Analysis complete. Retrieving data...", logFilename, 0, false, parameters);

                    //loop over folders and get logs
                    foreach (SampleRecord record in sampleSheet.getSampleRecords)
                    {
                        if (record.Analysis == @"P")
                        {
                            session.GetFiles(scratchDir + runID + @"/" + record.Sample_ID + "/*.sh.o*", localAnalysisDir + @"\").Check();
                            session.GetFiles(scratchDir + runID + @"/" + record.Sample_ID + "/*.sh.e*", localAnalysisDir + @"\").Check();
                        }
                    }

                    // Download files and throw on any error
                    session.GetFiles(scratchDir + runID + @"/" + runID + "_Filtered_Annotated.vcf", localAnalysisDir + @"\").Check();
                    session.GetFiles(scratchDir + runID + @"/" + runID + "_recalibration_plots.pdf", localAnalysisDir + @"\").Check();
                    session.GetFiles(scratchDir + runID + @"/BAMsforDepthAnalysis.list", localAnalysisDir + @"\").Check();
                    session.GetFiles(scratchDir + runID + @"/" + runID + "_Coverage.txt", localAnalysisDir + @"\").Check();
                    session.GetFiles(scratchDir + runID + @"/PreferredTranscripts.txt", localAnalysisDir + @"\").Check();
                    session.GetFiles(scratchDir + runID + "/*.bed", localAnalysisDir + @"\").Check();
                    session.GetFiles(scratchDir + runID + "/*.sh.o*", localAnalysisDir + @"\").Check();
                    session.GetFiles(scratchDir + runID + "/*.sh.e*", localAnalysisDir + @"\").Check();

                    //copy to network
                    File.Copy(localAnalysisDir + @"\" + runID + "_Filtered_Annotated.vcf", networkAnalysisDir + @"\" + runID + "_Filtered_Annotated.vcf");
                    File.Copy(localAnalysisDir + @"\" + runID + "_recalibration_plots.pdf", networkAnalysisDir + @"\" + runID + "_recalibration_plots.pdf");
                    File.Copy(localAnalysisDir + @"\BAMsforDepthAnalysis.list", networkAnalysisDir + @"\BAMsforDepthAnalysis.list");
                    File.Copy(localAnalysisDir + @"\" + runID + "_Coverage.txt", networkAnalysisDir + @"\" + runID + "_Coverage.txt");
                    File.Copy(localAnalysisDir + @"\PreferredTranscripts.txt", networkAnalysisDir + @"\PreferredTranscripts.txt");

                    //copy files to the network
                    foreach (var f in Directory.GetFiles(localAnalysisDir).Where(path => Regex.Match(path, @".*.bed").Success)) { File.Copy(f, networkAnalysisDir + @"\" + Path.GetFileName(f)); }
                    foreach (var f in Directory.GetFiles(localAnalysisDir).Where(path => Regex.Match(path, @".*.sh.o..*").Success)){ File.Copy(f, networkAnalysisDir + @"\" + Path.GetFileName(f)); }
                    foreach (var f in Directory.GetFiles(localAnalysisDir).Where(path => Regex.Match(path, @".*.sh.e..*").Success)) { File.Copy(f, networkAnalysisDir + @"\" + Path.GetFileName(f)); }

                    return true;
                }
                else
                {
                    AuxillaryFunctions.WriteLog(@"Marker File Not Found. Run is Pending", logFilename, 0, false, parameters);
                    return false;
                }
            }
        }

        private void WritePanelReport()
        {
            AuxillaryFunctions.WriteLog(@"Writing panel report...", logFilename, 0, false, parameters);

            string[] hgvs;
            GenomicVariant tempGenomicVariant;
            StreamWriter panelReport = new StreamWriter(localReportFilename);
            ParseVCF VCFFile = new ParseVCF(localAnalysisDir + @"\" + runID + "_Filtered_Annotated.vcf", logFilename, parameters);
            ParseVCF BCIntFile = new ParseVCF(parameters.getInterpretationsFile, logFilename, parameters);
            Dictionary<GenomicVariant, string> interpretations = new Dictionary<GenomicVariant, string>();

            panelReport.WriteLine("SampleID\tSampleName\tGene\tHGVSc\tHGVSp\tExon\tGenotype\tTranscriptID\tFunction\tInterpretation\tHTSFAmplicon\tChromosome\tPosition\tReference\tAlternative");

            //get failed amplicons
            AnalyseCoverageData(localAnalysisDir + @"\" + runID + "_Coverage.txt", localAnalysisDir + @"\BAMsforDepthAnalysis.list");

            //load BC ints
            foreach (VCFRecordWithGenotype record in BCIntFile.getVCFRecords[""]) //loop over interpretations
            {
                tempGenomicVariant.CHROM = record.CHROM;
                tempGenomicVariant.POS = record.POS;
                tempGenomicVariant.REF = record.REF;
                tempGenomicVariant.ALT = record.ALT;

                interpretations.Add(tempGenomicVariant, record.infoSubFields["INT"]);
            }

            //loop over samples for panel analysis
            foreach (SampleRecord record in sampleSheet.getSampleRecords)
            {
                if (record.Analysis != @"P"){
                    continue;
                }

                if (!VCFFile.getVCFRecords.ContainsKey(record.Sample_ID))
                {
                    AuxillaryFunctions.WriteLog(@"Sample " + record.Sample_ID + " has no genotypes in panel VCF" , logFilename, -1, false, parameters);
                    throw new FileLoadException();
                }

                //loop over VCF records for this sample
                foreach (VCFRecordWithGenotype VCFrecord in VCFFile.getVCFRecords[record.Sample_ID])
                {
                    tempGenomicVariant.CHROM = VCFrecord.CHROM;
                    tempGenomicVariant.REF = VCFrecord.REF;
                    tempGenomicVariant.ALT = VCFrecord.ALT;
                    tempGenomicVariant.POS = VCFrecord.POS;

                    if (VCFrecord.formatSubFields["GT"] == @"0/0") //skip hom ref variants
                    {
                        continue;
                    }

                    if (VCFFile.getSnpEffAnnotations.ContainsKey(tempGenomicVariant)) //annotation is available
                    {
                        //loop over annotations
                        foreach (Annotation ann in VCFFile.getSnpEffAnnotations[tempGenomicVariant])
                        {
                            hgvs = ann.Amino_Acid_Change.Split('/');

                            panelReport.Write(record.Sample_ID + "\t");
                            panelReport.Write(record.Sample_Name + "\t");
                            panelReport.Write(ann.Gene_Name + "\t");

                            if (hgvs.Length == 2){ //c. and p.
                                panelReport.Write(hgvs[1] + "\t"); //c.
                                panelReport.Write(hgvs[0] + "\t"); //p.
                            }
                            else
                            {
                                panelReport.Write(hgvs[0] + "\t"); //c.
                                panelReport.Write("\t");
                            }
                            
                            //print exon
                            panelReport.Write("{0}\t", ann.Exon_Rank);

                            if (VCFrecord.formatSubFields["GT"] == @"0/1")
                            {
                                panelReport.Write("HET\t");
                            } else if (VCFrecord.formatSubFields["GT"] == @"1/1"){
                                panelReport.Write("HOM_ALT\t");
                            }
                            else if (VCFrecord.formatSubFields["GT"] == @"")
                            {
                                panelReport.Write("Unknown\t");
                            }
                            else
                            {
                                panelReport.Write("Complex\t");
                            }

                            panelReport.Write(ann.Transcript_ID + "\t");
                            panelReport.Write(ann.Effect + "\t");

                            if (interpretations.ContainsKey(tempGenomicVariant)){
                                panelReport.Write(interpretations[tempGenomicVariant] + "\t");
                            }
                            else
                            {
                                panelReport.Write("\t");
                            }

                            //print HTSF amplicon (if present)
                            panelReport.Write(AuxillaryFunctions.LookupAmpliconID(new Tuple<string, UInt32>(VCFrecord.CHROM, VCFrecord.POS), coreBEDRecords.getBEDRecords) + "\t");

                            panelReport.Write(VCFrecord.CHROM + "\t");
                            panelReport.Write(VCFrecord.POS + "\t");
                            panelReport.Write(VCFrecord.REF + "\t");
                            panelReport.Write(VCFrecord.ALT + "\t");

                            panelReport.WriteLine();
                        }

                    }
                    else //no annotation
                    {
                        panelReport.Write(record.Sample_ID + "\t");
                        panelReport.Write(record.Sample_Name + "\t");
                        panelReport.Write("\t\t\t\t");

                        if (VCFrecord.formatSubFields["GT"] == @"0/1")
                        {
                            panelReport.Write("HET\t");
                        }
                        else if (VCFrecord.formatSubFields["GT"] == @"1/1")
                        {
                            panelReport.Write("HOM_ALT\t");
                        }
                        else
                        {
                            panelReport.Write("OTHER\t");
                        }

                        panelReport.Write("\t\t");

                        if (interpretations.ContainsKey(tempGenomicVariant))
                        {
                            panelReport.Write(interpretations[tempGenomicVariant] + "\t");
                        }
                        else
                        {
                            panelReport.Write("\t");
                        }

                        //print HTSF amplicon (if present)
                        panelReport.Write(AuxillaryFunctions.LookupAmpliconID(new Tuple<string, UInt32>(VCFrecord.CHROM, VCFrecord.POS), coreBEDRecords.getBEDRecords) + "\t");

                        panelReport.Write(VCFrecord.CHROM + "\t");
                        panelReport.Write(VCFrecord.POS + "\t");
                        panelReport.Write(VCFrecord.REF + "\t");
                        panelReport.Write(VCFrecord.ALT + "\t");

                        panelReport.WriteLine();
                    }

                } //done looping over VCF Records

                //add gaps
                if (!failedAmplicons.ContainsKey(record.Sample_ID))
                {
                    AuxillaryFunctions.WriteLog(@"Sample " + record.Sample_ID + " coverage data not loaded successfully!", logFilename, -1, false, parameters);
                    throw new FileLoadException();
                }

                //print failed regions
                foreach (string failedCoreRegion in failedAmplicons[record.Sample_ID])
                {
                    panelReport.WriteLine(record.Sample_ID + "\t" + record.Sample_Name + "\t\t\t\t\t\t\t\tGAP\t" + failedCoreRegion);
                }

            } //done looping over SampleIDs

            panelReport.Close();
        }

        private void AnalyseCoverageData(string samtoolsDepthFilePath, string samtoolsDepthSampleIDFilePath)
        {
            AuxillaryFunctions.WriteLog(@"Analysing coverage data...", logFilename, 0, false, parameters);

            string line, failedAmpliconID;
            UInt32 pos;
            Dictionary<Tuple<string, UInt32>, bool> isBaseCovered = new Dictionary<Tuple<string, UInt32>, bool>(); //bool = observed
            List<string> sampleIDs = new List<string>();

            //read sampleID order
            using (StreamReader fileInput = new StreamReader(samtoolsDepthSampleIDFilePath))
            {
                while ((line = fileInput.ReadLine()) != null)
                {
                    string[] fields = line.Split('_', '.');
                    sampleIDs.Add(fields[fields.Length - 2]);
                    failedAmplicons.Add(fields[fields.Length - 2], new HashSet<string>());
                }
            }

            //loop over target ROI, hash bases
            foreach (BEDRecord record in targetBEDRecords.getBEDRecords)
            {
                //iterate over region
                for (pos = record.start + 2; pos < record.end + 1; ++pos)
                {
                    if (!isBaseCovered.ContainsKey(new Tuple<string, UInt32>(record.chromosome, pos)))
                    {
                        isBaseCovered.Add(new Tuple<string, UInt32>(record.chromosome, pos), false);
                    }
                }

            }

            // loop over output and assign failed to low coverage amplicons
            using (StreamReader reader = new StreamReader(samtoolsDepthFilePath))
            {
                // Loop over the lines in the string.
                while ((line = reader.ReadLine()) != null)
                {
                    string[] fields = line.Split('\t');

                    pos = Convert.ToUInt32(fields[1]);

                    //mark base as observed in the dataset
                    isBaseCovered[new Tuple<string, UInt32>(fields[0], pos)] = true;

                    for (int n = 2; n < fields.Length; ++n) //skip chrom & pos
                    {
                        if (Convert.ToUInt32(fields[n]) < 30) //base has failed
                        {
                            //mark amplicon as failed
                            failedAmpliconID = AuxillaryFunctions.LookupAmpliconID(new Tuple<string, UInt32>(fields[0], pos), coreBEDRecords.getBEDRecords);

                            if (failedAmpliconID != "") //skip off target
                            {
                                failedAmplicons[sampleIDs[n - 2]].Add(failedAmpliconID);
                            }
                        }
                    }
                }
            }

            //report missing bases as failed
            foreach (KeyValuePair<Tuple<string, UInt32>, bool> nucl in isBaseCovered)
            {
                if (nucl.Value == false) //base not present in dataset
                {
                    failedAmpliconID = AuxillaryFunctions.LookupAmpliconID(nucl.Key, coreBEDRecords.getBEDRecords);

                    if (failedAmpliconID != "") //skip off target
                    {
                        foreach (string sampleID in sampleIDs)
                        {
                            //mark amplicon as failed
                            failedAmplicons[sampleID].Add(failedAmpliconID);
                        }
                    }

                }
            }
        }

        private void WriteVariablesFile()
        {
            AuxillaryFunctions.WriteLog(@"Writing variable files...", logFilename, 0, false, parameters);

            //concatinate analysisdirs
            StringBuilder AnalysisDirs = new StringBuilder(@"AnalysisDirs=( ");

            foreach (SampleRecord record in sampleSheet.getSampleRecords)
            {
                if (record.Analysis == @"P"){
                    AnalysisDirs.Append('"');
                    AnalysisDirs.Append(scratchDir);
                    AnalysisDirs.Append(runID);
                    AnalysisDirs.Append('/');
                    AnalysisDirs.Append(record.Sample_ID);
                    AnalysisDirs.Append('"');
                    AnalysisDirs.Append(' ');
                }
            }

            AnalysisDirs.Append(')');

            foreach (SampleRecord record in sampleSheet.getSampleRecords)
            {
                if (record.Analysis != @"P"){
                    continue;
                }

                //open variables output file
                StreamWriter VarFile = new StreamWriter(localAnalysisDir + @"\" + record.Sample_ID + @".variables");
             
                VarFile.Write("#Description: Pipeline Variables File\n");

                VarFile.Write("\n#Sample_ID\n");
                VarFile.Write(@"Sample_ID=" + record.Sample_ID + "\n");

                VarFile.Write("\n#FASTQ MD5 checksum\n");
                VarFile.Write(@"R1MD5Filename=" + Path.GetFileName(fastqFileNames[record.Sample_ID].Item1) + ".md5\n");
                VarFile.Write(@"R2MD5Filename=" + Path.GetFileName(fastqFileNames[record.Sample_ID].Item2) + ".md5\n");

                VarFile.Write("\n#FASTQ filenames\n");
                VarFile.Write(@"R1Filename=" + Path.GetFileName(fastqFileNames[record.Sample_ID].Item1) + "\n");
                VarFile.Write(@"R2Filename=" + Path.GetFileName(fastqFileNames[record.Sample_ID].Item2) + "\n");

                VarFile.Write("\n#Capture ROI\n");
                VarFile.Write(@"BEDFilename=" + Path.GetFileName(sampleSheet.getAnalyses[@"P"]) + "\n");

                VarFile.Write("\n#RunDetails\n");
                VarFile.Write(@"RunID=" + runID + "\n");
                VarFile.Write(@"ExperimentName=" + sampleSheet.getExperimentName + "\n");
                VarFile.Write("Platform=ILLUMINA\n");

                VarFile.Write("\n#Annotation\n");
                VarFile.Write(@"PreferredTranscriptsFile=" + scratchDir + runID + @"/" + Path.GetFileName(parameters.getPreferredTranscriptsFile) + "\n");

                VarFile.Write("\n#AnalysisFolders\n");
                VarFile.Write(AnalysisDirs.ToString() + "\n");

                VarFile.Close();
            }

        }

    }
}
