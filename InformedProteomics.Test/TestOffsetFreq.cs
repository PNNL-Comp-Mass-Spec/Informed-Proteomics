using System;
using System.Globalization;
using System.IO;
using System.Collections.Generic;
using System.Linq;
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
        private string _preTsv;
        private string _preRaw;
        private string _outPre;
        private string _outFileName;
        private double _pepQThreshold;
        private int _precursorCharge;
        private List<IonType> _ionTypes;
        private IonTypeFactory _ionTypeFactory;
        private ActivationMethod _act;
        double _relativeIntensityThreshold = 1.0;
        private bool _precursorOff; 
        
        private const string PrecChargeHeader = "Charge";
        private const string ScanHeader = "ScanNum";
        private const string PeptideHeader = "Peptide";
        private const string PepQValueHeader = "PepQValue";
        private const string FormulaHeader = "Formula";
        private readonly char[] _prefixes = {'a', 'b', 'c'};
        private readonly char[] _suffixes = {'x', 'y', 'z'};
        readonly Tolerance _defaultTolerance = new Tolerance(15, ToleranceUnit.Ppm);

        private IEnumerable<Tuple<string, Spectrum>> CleanScans(string txtFileName, string rawFileName)
        {
            var tsvParser = new TsvFileParser(txtFileName);

            var scans = tsvParser.GetData(ScanHeader);
            var peptides = tsvParser.GetData(PeptideHeader);
            var charges = tsvParser.GetData(PrecChargeHeader);
            var pepQValues = tsvParser.GetData(PepQValueHeader);
            var compositions = tsvParser.GetData(FormulaHeader);

//            var lcms = LcMsRun.GetLcMsRun(rawFileName, MassSpecDataType.XCaliburRun, 1.4826, 1.4826);
            var lcms = LcMsRun.GetLcMsRun(rawFileName, MassSpecDataType.XCaliburRun, 0, 0);

            var clean = new List<Tuple<string, Spectrum>>();
            var numRows = scans.Count;
            var peptideSet = new HashSet<string>();

            for (var i = 0; i < numRows; i++)
            {
                if (Convert.ToDouble(pepQValues[i]) > _pepQThreshold) continue;

                var precCharge = Convert.ToInt32(charges[i]);
                if (precCharge != _precursorCharge) continue;

                var peptide = peptides[i];
                if (!peptideSet.Add(peptide)) continue;

                var sequence = Sequence.GetSequenceFromMsGfPlusPeptideStr(peptide);
                var composition = Composition.Parse(compositions[i]);
                Assert.True(composition.Equals(sequence.Composition+Composition.H2O));
                var spec = lcms.GetSpectrum(Convert.ToInt32(scans[i]));
                clean.Add(new Tuple<string, Spectrum>(peptide, spec));
            }

            return clean;
        }

        private Dictionary<IonType, IonProbability> PrecursorOff(IEnumerable<Tuple<string, Spectrum>> cleanScans)
        {
            var probabilities = new Dictionary<IonType, IonProbability>();

            foreach (var ionType in _ionTypes)
            {
                if (!probabilities.ContainsKey(ionType))
                    probabilities.Add(ionType, new IonProbability(0, 0));
            }

            foreach (var node in cleanScans)
            {
                var peptide = node.Item1;
                var spectrum = node.Item2;

                var spec = spectrum as ProductSpectrum;
                if (spec == null) continue;
                if (spec.ActivationMethod != _act) continue;

                var sequence = Sequence.GetSequenceFromMsGfPlusPeptideStr(peptide);
                Composition comp = sequence.GetComposition();

                foreach (var ionType in _ionTypes)
                {
                    var ion = ionType.GetIon(comp);
                    ion.Composition.ComputeApproximateIsotopomerEnvelop();

                    probabilities[ionType].Total++;
                    if (spec.ContainsIon(ion, _defaultTolerance, _relativeIntensityThreshold))
                        probabilities[ionType].Found++;
                }
            }

            return probabilities;
        }

        private Dictionary<IonType, IonProbability> GetOffsetCounts(IEnumerable<Tuple<string, Spectrum>> cleanScans)
        {
            var probabilities = new Dictionary<IonType, IonProbability>();
            foreach (var ionType in _ionTypes)
            {
                if (!probabilities.ContainsKey(ionType))
                    probabilities.Add(ionType, new IonProbability(0, 0));
            }

            foreach (var node in cleanScans)
            {
                var protein = node.Item1;
                var spectrum = node.Item2;

                var sequence = Sequence.GetSequenceFromMsGfPlusPeptideStr(protein);
                var spec = spectrum as ProductSpectrum;
                if (spec == null) continue;
                for (int i = 1; i < protein.Length; i++)
                {
                    if (spec.ActivationMethod == _act)
                    {
                        foreach (var ionType in _ionTypes)
                        {
                            Composition sequenceComposition;
                            if (_prefixes.Contains(ionType.Name[0]))
                                sequenceComposition = sequence.GetComposition(0, i);
                            else if (_suffixes.Contains(ionType.Name[0]))
                                sequenceComposition = sequence.GetComposition(protein.Length - i, protein.Length);
                            else
                                throw new FormatException();

                            var ion = ionType.GetIon(sequenceComposition);
                            ion.Composition.ComputeApproximateIsotopomerEnvelop();

                            probabilities[ionType].Total++;
                            if (spec.ContainsIon(ion, _defaultTolerance, _relativeIntensityThreshold))
                                probabilities[ionType].Found++;

                            // Added by Sangtae for debugging
                            //Console.WriteLine("{0}{1} {2} {3}", ionTypeStr, i, ion.GetMonoIsotopicMz(), spec.ContainsIon(ion, _defaultTolerance, RelativeIntensityThreshold));
                        }
                    }
                }
            }
            return probabilities;
        }

        private void InitTest(INIReader reader)
        {
            // Read program variables
            var config = reader.getNodes("vars").First();
            _precursorCharge = Convert.ToInt32(config.Contents["precursorcharge"]);
            _pepQThreshold = Convert.ToDouble(config.Contents["pepqvalue"]);
            var actStr = config.Contents["activationmethod"].ToLower();
            switch (actStr)
            {
                case "hcd":
                    _act = ActivationMethod.HCD;
                    break;
                case "cid":
                    _act = ActivationMethod.CID;
                    break;
                case "etd":
                    _act = ActivationMethod.ETD;
                    break;
            }

            _precursorOff = (config.Contents.ContainsKey("precursorfrequencies") &&
                             config.Contents["precursorfrequencies"].ToLower() == "true");

            _relativeIntensityThreshold = Convert.ToDouble(config.Contents["relativeintensitythreshold"]);

            // Read ion data
            var ionInfo = reader.getNodes("ion").First();
            int totalCharges = Convert.ToInt32(ionInfo.Contents["totalcharges"]);
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
            _ionTypeFactory = new IonTypeFactory(ions, ionLosses, totalCharges);
            _ionTypes = _ionTypeFactory.GetAllKnownIonTypes().ToList();
            if (ionInfo.Contents.ContainsKey("exclusions"))
            {
                var ionExclusions = ionInfo.Contents["exclusions"].Split(',');
                foreach (var ionType in _ionTypes)
                {
                    if (ionExclusions.Contains(ionType.Name))
                        _ionTypes.Remove(ionType);
                }
            }

            // Read input and output file names
            var fileInfo = reader.getNodes("fileinfo").First();
            var name = fileInfo.Contents["name"];

            var tsvtemp = fileInfo.Contents["tsvpath"];
            _preTsv = tsvtemp.Replace("@", name);

            var rawtemp = fileInfo.Contents["rawpath"];
            _preRaw = rawtemp.Replace("@", name);

            var outPathtemp = fileInfo.Contents["outpath"];
            outPathtemp = outPathtemp.Replace("@", name);
            _outPre = outPathtemp.Replace("*", _precursorCharge.ToString(CultureInfo.InvariantCulture));

            var outFiletemp = fileInfo.Contents["outfile"];
            outFiletemp = outFiletemp.Replace("@", name);
            _outFileName = _outPre + outFiletemp.Replace("*", _precursorCharge.ToString(CultureInfo.InvariantCulture));
        }

        [Test]
        public void OffsetFreq()
        {
            InitTest(new INIReader(@"\\protoapps\UserData\Wilkins\ForChris\PrecursorOffConfig.ini"));

            var txtFiles = Directory.GetFiles(_preTsv).ToList();
            var rawFilesTemp = Directory.GetFiles(_preRaw).ToList();
            var rawFiles = rawFilesTemp.Where(rawFile => Path.GetExtension(rawFile) == ".raw").ToList();

            Assert.True(rawFiles.Count == txtFiles.Count);

            var found = new int[_ionTypes.Count];
            var total = new int[_ionTypes.Count];

            for (int i = 0; i < _ionTypes.Count; i++)
            {
                found[i] = 0;
                total[i] = 0;
            }

            for (int i=0;i<txtFiles.Count;i++)
            {
                string textFile = txtFiles[i];
                string rawFile = rawFiles[i];

                Console.WriteLine("{0}\t{1}", textFile, rawFile);
                var scans = CleanScans(textFile, rawFile);
                var cleanScans = scans as Tuple<string, Spectrum>[] ?? scans.ToArray();

                Dictionary<IonType, IonProbability> offsetCounts;

                if (_precursorOff)
                    offsetCounts = PrecursorOff(cleanScans);
                else
                    offsetCounts = GetOffsetCounts(cleanScans);

                for (int j = 0; j < _ionTypes.Count; j++)
                {
                    found[j] += offsetCounts[_ionTypes[j]].Found;
                    total[j] += offsetCounts[_ionTypes[j]].Total;
                }
            }

            using (var finalOutputFile = new StreamWriter(_outFileName))
            {
                for (int i = 0; i < _ionTypes.Count; i++)
                {
                    finalOutputFile.WriteLine("{0}\t{1}", _ionTypes[i].Name, Math.Round((double)(found[i]) / (total[i]), 5));
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
