using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace InformedProteomics.TopDown.Execution
{
    public class Ms1FeatureFinderInputParameter
    {
        public double MinSearchMass;
        public double MaxSearchMass;
        public int MinSearchCharge;
        public int MaxSearchCharge;
        public string InputPath;
        public bool ScoreReport;
        public bool CsvOutput;
        public bool TmpOutput;
        public int MaxThreads;
        public string OutputPath;
        //public bool Quant;

        public Ms1FeatureFinderInputParameter()
        {
            MinSearchMass = 3000;
            MaxSearchMass = 50000;
            MinSearchCharge = 2;
            MaxSearchCharge = 60;
            ScoreReport = false;
            CsvOutput = true;
            TmpOutput = false;
            MaxThreads = -1;
        }

        public Ms1FeatureFinderInputParameter(Dictionary<string, string> paramDic)
        {
            Parse(paramDic);
        }

        public void Parse(Dictionary<string, string> paramDic)
        {

            MinSearchMass = Math.Max(double.Parse(paramDic["-minMass"]), 100);
            MaxSearchMass = Math.Min(double.Parse(paramDic["-maxMass"]), 500000);


            MinSearchCharge = (int)Math.Max(double.Parse(paramDic["-minCharge"]), 2);
            MaxSearchCharge = (int)Math.Min(double.Parse(paramDic["-maxCharge"]), 60);
            InputPath = paramDic["-i"];
            OutputPath = paramDic["-o"];

            if (OutputPath == null) OutputPath = Path.GetDirectoryName(InputPath);

            MaxThreads = Int32.Parse(paramDic["-maxThreads"]);

            ScoreReport = Str2Bool(paramDic["-score"]);
            CsvOutput = Str2Bool(paramDic["-csv"]);
            TmpOutput = Str2Bool(paramDic["-tmp"]);
            //_quant = Str2Bool(_paramDic["-quant"]);

        }

        public void Display()
        {
            Console.WriteLine("InputPath\t{0}", InputPath);
            Console.WriteLine("OutputPath\t{0}", OutputPath);

            Console.WriteLine("MinMass\t{0}", MinSearchMass);
            Console.WriteLine("MaxMass\t{0}", MaxSearchMass);
            Console.WriteLine("MinCharge\t{0}", MinSearchCharge);
            Console.WriteLine("MaxCharge\t{0}", MaxSearchCharge);
            Console.WriteLine("ScoreReport\t{0}", ScoreReport ? "Y" : "N");
            Console.WriteLine("MaxThreads\t{0}", MaxThreads);
        }

        private static bool Str2Bool(string value)
        {
            return (value.Equals("y") || value.Equals("Y"));
        }

    }
}
