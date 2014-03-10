using System;
using System.Globalization;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
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
        public IonType Ion { get; private set; }

        public double Probability
        {
            get { return Math.Round((double) (Found)/(Total), 5); }
        }
        public IonProbability(int f, int t, IonType ion)
        {
            Found = f;
            Total = t;
            Ion = ion;
        }

        public static IonProbability operator +(IonProbability l, IonProbability r)
        {
            var added = new IonProbability(0, 0, l.Ion);
            added.Found = l.Found + r.Found;
            added.Total = l.Total + r.Total;
            return added;
        }
    }

    [TestFixture]
    public class OffsetTableTemp
    {
        private string[] _names;
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
        private bool _combineCharges;
        
        private const string PrecChargeHeader = "Charge";
        private const string ScanHeader = "ScanNum";
        private const string PeptideHeader = "Peptide";
        private const string PepQValueHeader = "PepQValue";
        private const string FormulaHeader = "Formula";
        private readonly char[] _prefixes = {'a', 'b', 'c'};
        private readonly char[] _suffixes = {'x', 'y', 'z'};
        readonly Tolerance _defaultTolerance = new Tolerance(15, ToleranceUnit.Ppm);

        private IEnumerable<Tuple<string, Spectrum>> CleanScans(LcMsRun lcms, string txtFileName, int charge)
        {
            var tsvParser = new TsvFileParser(txtFileName);

            var scans = tsvParser.GetData(ScanHeader);
            var peptides = tsvParser.GetData(PeptideHeader);
            var charges = tsvParser.GetData(PrecChargeHeader);
            var pepQValues = tsvParser.GetData(PepQValueHeader);
            var compositions = tsvParser.GetData(FormulaHeader);

            var clean = new List<Tuple<string, Spectrum>>();
            var numRows = scans.Count;
            var peptideSet = new HashSet<string>();

            for (var i = 0; i < numRows; i++)
            {
                if (Convert.ToDouble(pepQValues[i]) > _pepQThreshold) continue;

                var precCharge = Convert.ToInt32(charges[i]);
                if (precCharge != charge) continue;

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
                    probabilities.Add(ionType, new IonProbability(0, 0, ionType));
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
                    probabilities.Add(ionType, new IonProbability(0, 0, ionType));
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
                            //Console.WriteLine("{0}{1} {2} {3}", ionType.Name, i, ion.GetMonoIsotopicMz(), spec.ContainsIon(ion, _defaultTolerance, RelativeIntensityThreshold));
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

            _combineCharges = (config.Contents.ContainsKey("combinecharges") &&
                 config.Contents["combinecharges"].ToLower() == "true");

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
            var tempIonList = new List<IonType>();
            if (ionInfo.Contents.ContainsKey("exclusions"))
            {
                var ionExclusions = ionInfo.Contents["exclusions"].Split(',');
                foreach (var ionType in _ionTypes)
                {
                    if (!ionExclusions.Contains(ionType.Name))
                        tempIonList.Add(ionType);
                }
            }
            _ionTypes = tempIonList;


            // Read input and output file names
            var fileInfo = reader.getNodes("fileinfo").First();
            _names = fileInfo.Contents["name"].Split(',');

            _preTsv = fileInfo.Contents["tsvpath"];

            _preRaw = fileInfo.Contents["rawpath"];

            var outPathtemp = fileInfo.Contents["outpath"];
            _outPre = outPathtemp;
//            _outPre = outPathtemp.Replace("*", _precursorCharge.ToString(CultureInfo.InvariantCulture));

            var outFiletemp = fileInfo.Contents["outfile"];

            _outFileName = _outPre + outFiletemp;
//            _outFileName = _outPre + outFiletemp.Replace("*", _precursorCharge.ToString(CultureInfo.InvariantCulture));
        }

        [Test]
        public void OffsetFreq()
        {
            InitTest(new INIReader(@"\\protoapps\UserData\Wilkins\ForChris\OffsetFreqConfig.ini"));

            foreach (var name in _names)
            {
                var tsvName = _preTsv.Replace("@", name);
                var rawName = _preRaw.Replace("@", name);
                var txtFiles = Directory.GetFiles(tsvName).ToList();
                var rawFilesTemp = Directory.GetFiles(rawName).ToList();
                var rawFiles = rawFilesTemp.Where(rawFile => Path.GetExtension(rawFile) == ".raw").ToList();

                Assert.True(rawFiles.Count == txtFiles.Count);

                var probabilities = new Dictionary<string, IonProbability>[_precursorCharge];

                for (int i = 0; i < _precursorCharge; i++)
                {
                    probabilities[i] = new Dictionary<string, IonProbability>();
                    foreach (var ion in _ionTypes)
                    {
                        probabilities[i].Add(ion.Name, new IonProbability(0, 0, ion));
                    }
                }

                for (int i = 0; i < txtFiles.Count; i++)
                {
                    string textFile = txtFiles[i];
                    string rawFile = rawFiles[i];

                    Console.WriteLine("{0}\t{1}", textFile, rawFile);

//                  var lcms = LcMsRun.GetLcMsRun(rawFile, MassSpecDataType.XCaliburRun, 1.4826, 1.4826);
                    var lcms = LcMsRun.GetLcMsRun(rawFile, MassSpecDataType.XCaliburRun, 0, 0);

                    for (int j = 1; j <= _precursorCharge; j++)
                    {
                        var scans = CleanScans(lcms, textFile, j);
                        var cleanScans = scans as Tuple<string, Spectrum>[] ?? scans.ToArray();

                        Dictionary<IonType, IonProbability> offsetCounts;

                        if (_precursorOff)
                            offsetCounts = PrecursorOff(cleanScans);
                        else
                            offsetCounts = GetOffsetCounts(cleanScans);

                        foreach (var ion in _ionTypes)
                        {
                            string ionName = ion.Name;
                            probabilities[j - 1][ionName] += offsetCounts[ion];
                        }
                    }
                }

                var prob = new Dictionary<string, double>[_precursorCharge];
                for (int i = 0; i < _precursorCharge; i++)
                {
                    prob[i] = new Dictionary<string, double>();
                    foreach (var ionType in _ionTypes)
                    {
                        string ionName = ionType.Name;
                        prob[i].Add(ionType.Name, probabilities[i][ionType.Name].Probability);
                        if (_combineCharges && ionType.Charge > 1)
                        {
                            var chargeRemover = new StringBuilder(ionType.Name);
                            chargeRemover.Remove(1, 1);
                            ionName = chargeRemover.ToString();
                            prob[i][ionName] += probabilities[i][ionType.Name].Probability;
                        }
                    }
                }


                var outFile = _outFileName.Replace("@", name);

                for (int i = 0; i < _precursorCharge; i++)
                {
                    var outFileCharge = outFile.Replace("*", (i+1).ToString());
                    using (var finalOutputFile = new StreamWriter(outFileCharge))
                    {
                        foreach (var ion in _ionTypes)
                        {
                            if (_combineCharges && ion.Charge > 1) continue;
                            finalOutputFile.WriteLine("{0}\t{1}", ion.Name, prob[i][ion.Name]);
                        }
                    }
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
            var keyvalue = new Dictionary<String, String>();
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
