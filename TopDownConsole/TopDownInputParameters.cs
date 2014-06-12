using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Sequence;

namespace TopDownConsole
{
    public class TopDownInputParameters
    {
        public const string ParameterFileExtension = ".param";

        public string SpecFilePath { get; set; }
        public string DatabaseFilePath { get; set; }
        public string OutputDir { get; set; }
        public AminoAcidSet AminoAcidSet { get; set; }
        public int SearchMode { get; set; }
        public bool Tda { get; set; }
        public double PrecursorIonTolerancePpm { get; set; }
        public double ProductIonTolerancePpm { get; set; }
        public int MinSequenceLength { get; set; }
        public int MaxSequenceLength { get; set; }
        public int MinPrecursorIonCharge { get; set; }
        public int MaxPrecursorIonCharge { get; set; }
        public int MinProductIonCharge { get; set; }
        public int MaxProductIonCharge { get; set; }
        public double MinSequenceMass { get; set; }
        public double MaxSequenceMass { get; set; }

        private IEnumerable<SearchModification> _searchModifications;
        private int _maxNumDynModsPerSequence;

        public void Display()
        {
            Console.WriteLine("SpecFilePath: " + SpecFilePath);
            Console.WriteLine("DatabaseFilePath: " + DatabaseFilePath);
            Console.WriteLine("OutputDir: " + OutputDir);
            Console.WriteLine("SearchMode: " + SearchMode);
            Console.WriteLine("Tda: " + Tda);
            Console.WriteLine("PrecursorIonTolerancePpm: " + PrecursorIonTolerancePpm);
            Console.WriteLine("ProductIonTolerancePpm: " + ProductIonTolerancePpm);
            Console.WriteLine("MinSequenceLength: " + MinSequenceLength);
            Console.WriteLine("MaxSequenceLength: " + MaxSequenceLength);
            Console.WriteLine("MinPrecursorIonCharge: " + MinPrecursorIonCharge);
            Console.WriteLine("MaxPrecursorIonCharge: " + MaxPrecursorIonCharge);
            Console.WriteLine("MinProductIonCharge: " + MinProductIonCharge);
            Console.WriteLine("MaxProductIonCharge: " + MaxProductIonCharge);
            Console.WriteLine("MinSequenceMass: " + MinSequenceMass);
            Console.WriteLine("MaxSequenceMass: " + MaxSequenceMass);
            Console.WriteLine("MaxDynamicModificationsPerSequence: " + _maxNumDynModsPerSequence);
            Console.WriteLine("Modifications: ");
            foreach (var searchMod in _searchModifications)
            {
                Console.WriteLine(searchMod);
            }
        }

        public void Write()
        {
            var outputFilePath = OutputDir + Path.DirectorySeparatorChar +
                                       Path.GetFileNameWithoutExtension(SpecFilePath) + ParameterFileExtension;

            using (var writer = new StreamWriter(outputFilePath))
            {
                writer.WriteLine("SpecFile\t" + Path.GetFileName(SpecFilePath));
                writer.WriteLine("DatabaseFile\t" + Path.GetFileName(DatabaseFilePath));
                writer.WriteLine("SearchMode\t" + SearchMode);
                writer.WriteLine("Tda\t" + Tda);
                writer.WriteLine("PrecursorIonTolerancePpm\t" + PrecursorIonTolerancePpm);
                writer.WriteLine("ProductIonTolerancePpm\t" + ProductIonTolerancePpm);
                writer.WriteLine("MinSequenceLength\t" + MinSequenceLength);
                writer.WriteLine("MaxSequenceLength\t" + MaxSequenceLength);
                writer.WriteLine("MinPrecursorIonCharge\t" + MinPrecursorIonCharge);
                writer.WriteLine("MaxPrecursorIonCharge\t" + MaxPrecursorIonCharge);
                writer.WriteLine("MinProductIonCharge\t" + MinProductIonCharge);
                writer.WriteLine("MaxProductIonCharge\t" + MaxProductIonCharge);
                writer.WriteLine("MinSequenceMass\t" + MinSequenceMass);
                writer.WriteLine("MaxSequenceMass\t" + MaxSequenceMass);
                writer.WriteLine("MaxDynamicModificationsPerSequence\t" + _maxNumDynModsPerSequence);
                foreach (var searchMod in _searchModifications)
                {
                    writer.WriteLine("Modification\t"+searchMod);
                }
            }
        }

        public string Parse(Dictionary<string, string> parameters)
        {
            var message = CheckIsValid(parameters);
            if (message != null) return message;

            SpecFilePath = parameters["-s"];
            DatabaseFilePath = parameters["-d"];
            OutputDir = parameters["-o"] ?? Environment.CurrentDirectory;

            var modFilePath = parameters["-mod"];
            if (modFilePath != null)
            {
                var parser = new ModFileParser(modFilePath);
                _searchModifications = parser.SearchModifications;
                _maxNumDynModsPerSequence = parser.MaxNumDynModsPerSequence;

                if (_searchModifications == null) return "Error while parsing " + modFilePath + "!";

                AminoAcidSet = new AminoAcidSet(_searchModifications, _maxNumDynModsPerSequence);
            }
            else
            {
                AminoAcidSet = new AminoAcidSet();
            }

            SearchMode = Convert.ToInt32(parameters["-m"]);
            if (SearchMode < 0 || SearchMode > 2)
            {
                return "Invalid value (" + SearchMode + ") for parameter -m";
            }

            PrecursorIonTolerancePpm = Convert.ToDouble(parameters["-t"]);
            ProductIonTolerancePpm = Convert.ToDouble(parameters["-f"]);

            var tdaVal = Convert.ToInt32(parameters["-tda"]);
            if (tdaVal != 0 && tdaVal != 1)
            {
                return "Invalid value (" + tdaVal + ") for parameter -tda";
            }
            Tda = (tdaVal == 1);

            MinSequenceLength = Convert.ToInt32(parameters["-minLength"]);
            MaxSequenceLength = Convert.ToInt32(parameters["-maxLength"]);
            if (MinSequenceLength > MaxSequenceLength)
            {
                return "MinSequenceLength (" + MinSequenceLength + ") is larger than MaxSequenceLength (" + MaxSequenceLength + ")!";
            }

            MinPrecursorIonCharge = Convert.ToInt32(parameters["-minCharge"]);
            MaxPrecursorIonCharge = Convert.ToInt32(parameters["-maxCharge"]);
            if (MinSequenceLength > MaxSequenceLength)
            {
                return "MinPrecursorCharge (" + MinPrecursorIonCharge + ") is larger than MaxPrecursorCharge (" + MaxPrecursorIonCharge + ")!";
            }

            MinProductIonCharge = Convert.ToInt32(parameters["-minFragCharge"]);
            MaxProductIonCharge = Convert.ToInt32(parameters["-maxFragCharge"]);
            if (MinSequenceLength > MaxSequenceLength)
            {
                return "MinFragmentCharge (" + MinProductIonCharge + ") is larger than MaxFragmentCharge (" + MaxProductIonCharge + ")!";
            }

            MinSequenceMass = Convert.ToDouble(parameters["-minMass"]);
            MaxSequenceMass = Convert.ToDouble(parameters["-maxMass"]);
            if (MinSequenceMass > MaxSequenceMass)
            {
                return "MinSequenceMassInDa (" + MinSequenceMass + ") is larger than MaxSequenceMassInDa (" + MaxSequenceMass + ")!";
            }

            return null;
        }

        static string CheckIsValid(Dictionary<string, string> parameters)
        {
            foreach (var keyValuePair in parameters)
            {
                var key = keyValuePair.Key;
                var value = keyValuePair.Value;
                if (keyValuePair.Value == null && keyValuePair.Key != "-mod")
                {
                    return "Missing required parameter " + key + "!";
                }

                if (key.Equals("-s"))
                {
                    if (value == null)
                    {
                        return "Missing parameter " + key + "!";
                    }
                    if (!File.Exists(value))
                    {
                        return "File not found: " + value + "!";
                    }
                    var extension = Path.GetExtension(value);
                    if (!Path.GetExtension(value).ToLower().Equals(".raw"))
                    {
                        return "Invalid extension for the parameter " + key + " (" + extension + ")!";
                    }
                }
                else if(key.Equals("-d"))
                {
                    if (value == null)
                    {
                        return "Missing parameter " + key + "!";
                    }
                    if (!File.Exists(value))
                    {
                        return "File not found." + value + "!"; ;
                    }
                    var extension = Path.GetExtension(value).ToLower();
                    if (!extension.Equals(".fa") && !extension.Equals(".fasta"))
                    {
                        return "Invalid extension for the parameter " + key + " (" + extension + ")!";
                    }
                }
                else if (key.Equals("-o"))
                {
                    
                }
                else if (key.Equals("-mod"))
                {
                    if (value != null && !File.Exists(value))
                    {
                        return "File not found." + value + "!";
                    }
                }
                else
                {
                    double num;
                    if (!double.TryParse(value, out num))
                    {
                        return "Invalid value (" + value + ") for the parameter " + key + "!";
                    }
                }
            }

            return null;
        }
    }
}
