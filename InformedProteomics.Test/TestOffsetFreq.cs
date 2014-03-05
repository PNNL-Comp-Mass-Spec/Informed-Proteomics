using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    public class IonProbability
    {
        public int Found { get; set; }
        public int Total { get; set; }
        public IonProbability(int f, int t) { Found = f; Total = t; }
    }

    [TestFixture]
    public class OffsetTableTemp
    {
        private string fileList;
        private string preRes;
        private string preRaw;
        private string outPre;
        private string outFileName;
        private string outSuff;
        private double PepQThreshold;
        private int PrecursorCharge;
        private string[] ionTypes;
        private IonTypeFactory ionTypeFactory;
        private ActivationMethod act;


        private const string PrecChargeHeader = "Charge";
        private const string ScanHeader = "ScanNum";
        const double RelativeIntensityThreshold = 1.0;
        Tolerance defaultTolerance = new Tolerance(15, ToleranceUnit.Ppm);

        private IEnumerable<Tuple<string, Spectrum>> CleanScans(string txtFileName, string rawFileName)
        {
            var tsvParser = new TsvFileParser(txtFileName);

            var scans = tsvParser.GetData(ScanHeader);
            var peptides = tsvParser.GetPeptides(PepQThreshold);
            var charges = tsvParser.GetData(PrecChargeHeader);
            var lcms = LcMsRun.GetLcMsRun(rawFileName, MassSpecDataType.XCaliburRun, 1.4826, 1.4826);

            var clean = new List<Tuple<string, Spectrum>>();
            using (var chargei = charges.GetEnumerator())
            using (var scani = scans.GetEnumerator())
            using (var peptidei = peptides.GetEnumerator())
            {
                while (scani.MoveNext() && peptidei.MoveNext() && chargei.MoveNext())
                {
                    var spec = lcms.GetSpectrum(Convert.ToInt32(scani.Current));
                    int precCharge = Convert.ToInt32(chargei.Current);
                    if (precCharge != PrecursorCharge) continue;
                    clean.Add(new Tuple<string, Spectrum>(peptidei.Current, spec));
                }
            }
            return clean;
        }

        private Dictionary<string, IonProbability> GetOffsetCounts(IEnumerable<Tuple<string, Spectrum>> cleanScans,
                                        string[] ionTypes, ActivationMethod act, IonTypeFactory ionTypeFactory)
        {
            var probabilities = new Dictionary<string, IonProbability>();
            foreach (string ionTypeStr in ionTypes)
            {
                if (!probabilities.ContainsKey(ionTypeStr))
                    probabilities.Add(ionTypeStr, new IonProbability(0, 0));
            }

            var pepDict = new Dictionary<string, int>();

            foreach (var node in cleanScans)
            {
                var protein = node.Item1;
                var spectrum = node.Item2;
                if (pepDict.ContainsKey(protein)) continue;
                else
                    pepDict.Add(protein, 0);

                var sequence = Sequence.GetSequenceFromMsGfPlusPeptideStr(protein);
                var spec = spectrum as ProductSpectrum;
                if (spec == null) continue;
                for (int i = 0; i < protein.Length - 1; i++)
                {
                    if (spec.ActivationMethod == act)
                    {
                        foreach (var ionTypeStr in ionTypes)
                        {
                            Composition sequenceComposition = null;
                            if (ionTypeStr[0] == 'a' || ionTypeStr[0] == 'b' || ionTypeStr[0] == 'c')
                                sequenceComposition = sequence.GetComposition(0, i);
                            else if (ionTypeStr[0] == 'x' || ionTypeStr[0] == 'y' || ionTypeStr[0] == 'z')
                                sequenceComposition = sequence.GetComposition(protein.Length - i, protein.Length);
                            else
                                throw new FormatException();

                            var ionType = ionTypeFactory.GetIonType(ionTypeStr);
                            var ion = ionType.GetIon(sequenceComposition);
                            ion.Composition.ComputeApproximateIsotopomerEnvelop();

//                            var mz = ion.GetMonoIsotopicMz();
//                            var peak = spec.FindPeak(mz, defaultTolerance);

                            probabilities[ionTypeStr].Total++;
                            if (spec.ContainsIon(ion, defaultTolerance, RelativeIntensityThreshold))
//                            if (peak != null)
                                probabilities[ionTypeStr].Found++;
                        }
                    }
                }
            }
            return probabilities;
        }

        private void WriteProbFile(string outFile, Dictionary<string, IonProbability> offsetCounts, string[] ionTypes)
        {
            using (var file = new StreamWriter(outFile))
            {
                for (int i = 0; i < ionTypes.Length; i++)
                {
                    file.WriteLine("{0}\t{1}", ionTypes[i], Math.Round((double)(offsetCounts[ionTypes[i]].Found) / (offsetCounts[ionTypes[i]].Total), 5));
                }
            }
        }

        private void WriteOffsetCountsFile(string outFile, Dictionary<string, IonProbability> offsetCounts, string[] ionTypes)
        {
            using (var file = new StreamWriter(outFile))
            {
                for (int i = 0; i < ionTypes.Length; i++)
                {
                    file.Write("{0}Found\t{0}Total", ionTypes[i]);
                    if (i != ionTypes.Length - 1)
                        file.Write("\t");
                }
                file.WriteLine();
                for (int i = 0; i < ionTypes.Length; i++)
                {
                    file.Write("{0}\t{1}", offsetCounts[ionTypes[i]].Found, offsetCounts[ionTypes[i]].Total);
                    if (i != ionTypes.Length - 1)
                        file.Write("\t");
                }
            }
        }

        private void InitTest(INIReader reader)
        {
            // Read program variables
            var config = reader.getNodes("vars").First();
            PrecursorCharge = Convert.ToInt32(config.Contents["precursorcharge"]);
            PepQThreshold = Convert.ToDouble(config.Contents["pepqvalue"]);
            var actStr = config.Contents["activationmethod"].ToLower();
            switch (actStr)
            {
                case "hcd":
                    act = ActivationMethod.HCD;
                    break;
                case "cid":
                    act = ActivationMethod.CID;
                    break;
                case "etd":
                    act = ActivationMethod.ETD;
                    break;
            }

            // Read ion data
            var ionInfo = reader.getNodes("ion").First();
            int totalCharges = Convert.ToInt32(ionInfo.Contents["totalcharges"]);
            var ionNames = ionInfo.Contents["ions"].Split(',');
            ionTypes = new string[ionNames.Length];
            Array.Copy(ionNames, ionTypes, ionTypes.Length);
            var ionTypeStr = ionInfo.Contents["iontype"].Split(',');
            var ions = new BaseIonType[ionTypeStr.Length];
            for (int i = 0; i < ionTypeStr.Length; i++)
            {
                switch (ionTypeStr[i].ToLower())
                {
                    case "a":
                        ions[i] = BaseIonType.A;
                        break;
                    case "b":
                        ions[i] = BaseIonType.B;
                        break;
                    case "c":
                        ions[i] = BaseIonType.C;
                        break;
                    case "x":
                        ions[i] = BaseIonType.X;
                        break;
                    case "y":
                        ions[i] = BaseIonType.Y;
                        break;
                    case "z":
                        ions[i] = BaseIonType.Z;
                        break;
                }
            }
            var ionLossStr = ionInfo.Contents["losses"].Split(',');
            var ionLosses = new NeutralLoss[ionLossStr.Length];
            for (int i = 0; i < ionLossStr.Length; i++)
            {
                switch (ionLossStr[i].ToLower())
                {
                    case "noloss":
                        ionLosses[i] = NeutralLoss.NoLoss;
                        break;
                    case "nh3":
                        ionLosses[i] = NeutralLoss.NH3;
                        break;
                    case "h2o":
                        ionLosses[i] = NeutralLoss.H2O;
                        break;
                }
            }
            ionTypeFactory = new IonTypeFactory(ions, ionLosses, totalCharges);

            // Read input and output file names
            var fileInfo = reader.getNodes("fileinfo").First();
            var name = fileInfo.Contents["name"];
            var fileListtemp = fileInfo.Contents["filelist"];
            fileList = fileListtemp.Replace("@", name);

            var tsvtemp = fileInfo.Contents["tsvpath"];
            preRes = tsvtemp.Replace("@", name);

            var rawtemp = fileInfo.Contents["rawpath"];
            preRaw = rawtemp.Replace("@", name);

            var outPathtemp = fileInfo.Contents["outpath"];
            outPathtemp = outPathtemp.Replace("@", name);
            outPre = outPathtemp.Replace("*", PrecursorCharge.ToString());

            var outFiletemp = fileInfo.Contents["outfile"];
            outFiletemp = outFiletemp.Replace("@", name);
            outFileName = outPre + outFiletemp.Replace("*", PrecursorCharge.ToString());

            var outSufftemp = fileInfo.Contents["outsuff"];
            outSufftemp = outSufftemp.Replace("@", name);
            outSuff = outSufftemp.Replace("*", PrecursorCharge.ToString());
        }

        [Test]
        public void OffsetFreq()
        {
            InitTest(new INIReader(@"C:\Users\wilk011\Documents\DataFiles\OffsetFreqConfig.ini"));

            var fileNameParser = new TsvFileParser(fileList);

            var txtFiles = fileNameParser.GetData("text");
            var rawFiles = fileNameParser.GetData("raw");

            using (var txtFileIt = txtFiles.GetEnumerator())
            using (var rawFileIt = rawFiles.GetEnumerator())
            {
                while (txtFileIt.MoveNext() && rawFileIt.MoveNext())
                {
                    string textFile = preRes + txtFileIt.Current;
                    string rawFile = preRaw + rawFileIt.Current;
                    Console.WriteLine(rawFile);
                    var scans = CleanScans(textFile, rawFile);
                    var cleanScans = scans as Tuple<string, Spectrum>[] ?? scans.ToArray();

                    var offsetCounts = GetOffsetCounts(cleanScans, ionTypes, act, ionTypeFactory);
                    var outFile = outPre + rawFileIt.Current + ".Charge" + PrecursorCharge + outSuff;
                    WriteOffsetCountsFile(outFile, offsetCounts, ionTypes);
                    WriteProbFile(outFile + ".prob", offsetCounts, ionTypes);
                }
            }

            // Consolidate files
            var found = new int[ionTypes.Length];
            var total = new int[ionTypes.Length];
            foreach (var rawFile in rawFiles)
            {
                string fileName = outPre + rawFile + ".Charge" + PrecursorCharge + outSuff;
                var outReader = new TsvFileParser(fileName);
                for (int i = 0; i < ionTypes.Length; i++)
                {
                    var foundStr = outReader.GetData(ionTypes[i] + "Found")[0];
                    var totalStr = outReader.GetData(ionTypes[i] + "Total")[0];
                    found[i] = Convert.ToInt32(foundStr);
                    total[i] = Convert.ToInt32(totalStr);
                }
            }
            using (var finalOutputFile = new StreamWriter(outFileName))
            {
                for (int i = 0; i < ionTypes.Length; i++)
                {
                    finalOutputFile.WriteLine("{0}\t{1}", ionTypes[i], Math.Round((double)(found[i]) / (total[i]), 5));
                }
            }
        }
    }

    class Node
    {
        public string Header { get; private set; }      // The header tag

        // All of the key/value pairs:
        public Dictionary<String, String> Contents { get; private set; }

        // constructor
        public Node(string header, Dictionary<String, String> contents)
        {
            Header = header;
            Contents = contents;
        }
    }
    class INIReader
    {
        public List<Node> Nodes { get; protected set; }

        // Exception
        public class InvalidHeader : Exception { }

        private bool ValidHeader(string header)
        {
            return (header[0] == '[' && header[header.Length - 1] == ']');
        }

        /* read() Read the file and store the result
         * in Nodes.
         */
        private void read(string fileName)
        {
            Node currentNode = null;
            Dictionary<String, String> keyvalue = new Dictionary<String, String>();
            string[] lines = System.IO.File.ReadAllLines(fileName);
            char[] headerbrackets = { '[', ']' };
            string header = "";
            foreach (var line in lines)
            {
                string commentsStripped = line.Split('#')[0];      // remove comments
                string[] parts = commentsStripped.Split('=');       // split key/value
                if (parts.Length < 2)
                {
                    // The line is either a header, empty line, or invalid
                    parts[0] = parts[0].Trim().ToLower();
                    if (parts[0] == "")
                        // empty line
                        continue;
                    else if (currentNode == null)
                    {
                        // first node in the file
                        currentNode = new Node(null, null);
                        header = parts[0].Trim(headerbrackets);
                    }
                    else if (currentNode != null)
                    {
                        // this isn't the first node in the file
                        currentNode = new Node(header, keyvalue);
                        keyvalue = new Dictionary<String, String>();
                        Nodes.Add(currentNode);
                        header = parts[0].Trim(headerbrackets);
                    }
                    if (!ValidHeader(parts[0]))
                        // invalid header
                        throw new InvalidHeader();
                }
                else
                {
                    // key value pair
                    string key = parts[0].Trim().ToLower();
                    string value = parts[1].Trim();
                    keyvalue.Add(key, value);
                }
            }
            currentNode = new Node(header, keyvalue);
            Nodes.Add(currentNode);
        }

        /*
         *  Constructor
         */
        public INIReader(string fileName)
        {
            Nodes = new List<Node>();
            read(fileName);
        }


        /*
         * getNodes() return a list of all the nodes with a particular
         * header tag.
         */
        public List<Node> getNodes(string headerTag)
        {
            return (from i in Nodes
                    where i.Header == headerTag
                    select i).ToList();
        }
    }
}
