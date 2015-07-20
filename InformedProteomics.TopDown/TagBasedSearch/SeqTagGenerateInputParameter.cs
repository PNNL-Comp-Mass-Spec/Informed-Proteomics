using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace InformedProteomics.TopDown.TagBasedSearch
{
    public class SeqTagGenerateInputParameter
    {
        public string InputPath;
        public string OutputPath;
        public double TolerancePpm;
        public int MinSeqTagLength;
        public int MaxSeqTagCount;

        public SeqTagGenerateInputParameter()
        {
            TolerancePpm = 5;
            MinSeqTagLength = 4;

        }

        public SeqTagGenerateInputParameter(Dictionary<string, string> paramDic)
        {
            Parse(paramDic);
        }

        public void Parse(Dictionary<string, string> paramDic)
        {
        
            InputPath = paramDic["-i"];
            OutputPath = paramDic["-o"];

        }
        /*
        public void Display()
        {
            Console.WriteLine("InputPath\t{0}", InputPath);
            Console.WriteLine("OutputPath\t{0}", OutputPath);

            Console.WriteLine("Tolerance\t{0}", MinSearchMass);
            Console.WriteLine("MaxMass\t{0}", MaxSearchMass);
            Console.WriteLine("MinCharge\t{0}", MinSearchCharge);
            Console.WriteLine("MaxCharge\t{0}", MaxSearchCharge);
            
            Console.WriteLine("ScoreReport\t{0}", ScoreReport ? "Y" : "N");

            Console.WriteLine("LikelihoodRatioThreshold\t{0}", LikelihoodScoreThreshold);

            Console.WriteLine("MaxThreads\t{0}", MaxThreads);


 "Usage: " + Name + ".exe\n" +
                "\t[-i InputFolder or InputFile]\n" +
                "\t[-o OutFolder (default: InputFolder)]\n" +
                "\t[-t Tolerance (default: 5 ppm)]\n"+
                "\t[-minLen MinSequenceTagLength] (minimum length of sequence tag, default: 5)\n" +
                "\t[-maxTags MaxNumberOfSequenceTags] (maximum number of sequence tags per spectrum, default: -1, unlimited: -1)\n" +
                "\t[-maxThreads 0 (default: 0 (no limit))]\n"

        }

        private static bool Str2Bool(string value)
        {
            return (value.Equals("y") || value.Equals("Y"));
        }*/
    }
}
