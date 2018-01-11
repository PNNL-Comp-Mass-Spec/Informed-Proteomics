using System;
using System.Collections.Generic;
using System.IO;

namespace InformedProteomics.TopDown.Quantification
{
    public class MSDeconvReader
    {
        public MSDeconvReader(double minMass = double.MinValue, double maxMass = double.MaxValue, int minCharge = Int32.MinValue, int maxCharge = int.MaxValue, int minScan = Int32.MinValue, int maxScan = Int32.MaxValue)
        {
            MinMass = minMass;
            MaxMass = maxMass;
            MinCharge = minCharge;
            MaxCharge = maxCharge;
            MinScanNum = minScan;
            MaxScanNum = maxScan;
        }

        public List<MSDeconvNode> GetDeconvNodesForMsDeconv(String filePath)
        {
            var txtFile = String.Format(@"{0}", filePath);
            var lineInfoList = new List<MSDeconvNode>();

            if (!File.Exists(txtFile))
            {
                Console.WriteLine(@"Error: {0} not found", txtFile);
                return lineInfoList;
            }

            var file = new StreamReader(txtFile);
            var line = "";
            var isInfoLine = false;
            while ((line = file.ReadLine()) != null)
            {
                //checks for ending of info
                if (String.IsNullOrWhiteSpace(line))
                {
                    isInfoLine = false;
                    continue;
                }

                //checks for start of important information
                var splitString = line.Split(' ');
                if (splitString[0] == "Ms")
                {
                    var msLevel = Int32.Parse(splitString[splitString.Length - 1]);
                    if (MsLevel == 0 || msLevel == MsLevel)
                    {
                        isInfoLine = true;
                        file.ReadLine();
                        if (msLevel == 2) file.ReadLine();
                        continue;
                    }
                }
                //Process line
                if (isInfoLine)
                {
                    splitString = line.Split('\t');
                    var infoLine = GetDeconvLineMsDeconv(splitString);
                    if(infoLine != null) lineInfoList.Add(infoLine);
                }
            }

            file.Close();
            return lineInfoList;
        }

        public List<MSDeconvNode> GetDeconvNodesForDecon2Ls(String filePath)
        {
            var nodeList = new List<MSDeconvNode>();

            if (!File.Exists(filePath))
            {
                Console.WriteLine(@"Error: {0} not found", filePath);
                return nodeList;
            }

            var lineList = File.ReadAllText(filePath).Split('\n');
            Console.WriteLine(lineList.Length);
            for (var i = 1; i < lineList.Length-1; i++)
            {
                var line = lineList[i].Split(',');
                var node = GetDeconvLineDecon2Ls(line);
                if (node != null) nodeList.Add(node);
            }

            return nodeList;
        }

        private MSDeconvNode GetDeconvLineMsDeconv(String[] line)
        {
            var scanNum = Int32.Parse(line[0]);
            var charge = Int32.Parse(line[4]);
            var mass =  Double.Parse(line[11]);
            var intensity = Double.Parse(line[13]);

            if (Double.IsNaN(intensity)) intensity = 0;

            if (scanNum < MinScanNum || scanNum > MaxScanNum) return null;
            if (charge < MinCharge || charge > MaxCharge) return null;
            if (mass < MinMass || mass > MaxMass) return null;

            return new MSDeconvNode(scanNum,mass,intensity,charge);
        }

        private MSDeconvNode GetDeconvLineDecon2Ls(String[] line)
        {
            var scanNum = Int32.Parse(line[0]);
            var charge = Int32.Parse(line[1]);
            var mass = Double.Parse(line[6]);
            var intensity = Double.Parse(line[10]);

            if (scanNum < MinScanNum || scanNum > MaxScanNum) return null;
            if (charge < MinCharge || charge > MaxCharge) return null;
            if (mass < MinMass || mass > MaxMass) return null;

            return new MSDeconvNode(scanNum,mass,intensity,charge);
        }

        public double MinMass { get; set; }
        public double MaxMass { get; set; }
        public int MinCharge { get; set; }
        public double MaxCharge { get; set; }
        public double MinScanNum { get; set; }
        public double MaxScanNum { get; set; }
        public int MsLevel { get; set; }
    }
}
