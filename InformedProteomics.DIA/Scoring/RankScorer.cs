using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.DIA.Scoring
{
    public class RankScorer
    {
        public RankScorer(string paramFilePath) : this(paramFilePath, false)
        {
            
        }

        public RankScorer(string paramFilePath, bool verbose)
        {
            ReadFromFile(paramFilePath, verbose);
        }

        private Dictionary<int, IonType[]> _ionTypeTable;

        // spectral data type
        private int _numSegments = 1;
        private Dictionary<int, int> _chargeHistogram;  // charge, #spectra
        private SortedSet<Partition> _partitionSet;
        private Dictionary<int, IList<PrecursorOffsetFrequency>> _precursorOffMap; // charge -> precursor offsets to be removed
        private Dictionary<Partition, IEnumerable<FragmentOffsetFrequency>> _fragmentOffMap; // partition -> ion types
        private Dictionary<Partition, IEnumerable<FragmentOffsetFrequency>> _noiseOffMap; // partition -> ion types (noise)
        private Dictionary<Partition, Dictionary<IonType,float[]>> _rankDistTable; // rank distributions

        //private Tolerance _mme; // mass measurement error

        // Deconvolution
        private bool _applyDeconvolution;   // MS/MS scans are deconvoluted if true
        private float _deconvolutionErrorToleranceDa;   // Tolerance in Da for deconvolution

        private int _numPrecursorOffsets;
        private int _maxRank;

        private void ReadFromFile(string paramFile, Boolean verbose)
        {
            if (!File.Exists(paramFile))
                return;

            using (var reader = new BinaryReader(File.Open(paramFile, FileMode.Open)))
            {
                // Version
                var version = reader.ReadInt32();
                if(verbose) Console.WriteLine("Version: {0}", version);

                // Activation method
                var actMethod = reader.ReadChars(reader.ReadByte());

                // Instrument type
                var instType = reader.ReadChars(reader.ReadByte());

                // Enzyme
                var enzyme = reader.ReadChars(reader.ReadByte());

                // Protocol
                var protocol = reader.ReadChars(reader.ReadByte());

                if (verbose) Console.WriteLine("Spectral type: {0}_{1}_{2}_{3}", actMethod, instType, enzyme, protocol);

                // TODO: set up spectral type

                // Tolerance
                var isTolerancePpm = reader.ReadBoolean();
                var mmeVal = reader.ReadSingle();
                // Ignore mme
                //_mme = new Tolerance(mmeVal, isTolerancePpm ? ToleranceUnit.Ppm : ToleranceUnit.Dalton);

                // Deconvolution information
                _applyDeconvolution = reader.ReadBoolean();
                _deconvolutionErrorToleranceDa = reader.ReadSingle();
                if (verbose) Console.WriteLine("Apply deconvolution: {0}, Tolerance: {1}", _applyDeconvolution, _deconvolutionErrorToleranceDa);

                // Charge histogram
                if (verbose) Console.WriteLine("Charge Histogram");
                var chargeHistSize = reader.ReadInt32();
                _chargeHistogram = new Dictionary<int, int>();
                for (var i = 0; i < chargeHistSize; i++)
                {
                    var charge = reader.ReadInt32();
                    var numSpecs = reader.ReadInt32();
                    _chargeHistogram[charge] = numSpecs;
                    if(verbose) Console.WriteLine(charge+"\t"+numSpecs);
                }

                // Partition information
                if (verbose) Console.WriteLine("Partition Information"); 
                var numPartitions = reader.ReadInt32();
                _numSegments = reader.ReadInt32();
                _partitionSet = new SortedSet<Partition>();

                for (var i = 0; i < numPartitions; i++)
                {
                    var charge = reader.ReadInt32();
                    var neutralPeptideMass = reader.ReadInt32();
                    var segIndex = reader.ReadInt32();
                    _partitionSet.Add(new Partition(charge, neutralPeptideMass, segIndex));
                    if(verbose) Console.WriteLine("{0}\t{1}\t{2}", charge, neutralPeptideMass, segIndex);
                }

                // Precursor offset frequency function
                if (verbose) Console.WriteLine("Precursor Offset Frequency Function");
                _precursorOffMap = new Dictionary<int, IList<PrecursorOffsetFrequency>>();
                _numPrecursorOffsets = reader.ReadInt32();
                for (var i = 0; i < _numPrecursorOffsets; i++)
                {
                    var charge = reader.ReadInt32();
                    var reducedCharge = reader.ReadInt32();
                    var offset = reader.ReadSingle();
                    var isTolPpm = reader.ReadBoolean();
                    var tolVal = reader.ReadSingle();
                    var tolerance = new Tolerance(tolVal, isTolPpm ? ToleranceUnit.Ppm : ToleranceUnit.Dalton);
                    var frequency = reader.ReadSingle();

                    IList<PrecursorOffsetFrequency> offList;
                    if (!_precursorOffMap.TryGetValue(charge, out offList))
                    {
                        offList = new List<PrecursorOffsetFrequency>();
                        _precursorOffMap[charge] = offList;
                    }
                    offList.Add(new PrecursorOffsetFrequency(reducedCharge, offset, frequency, tolerance));
                    if (verbose) Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", 
                        charge, reducedCharge, offset, tolerance, frequency);
                }

                // Fragment ion offset frequency function
                if (verbose) Console.WriteLine("Fragment Offset Frequency Function");
                _fragmentOffMap = new Dictionary<Partition, IEnumerable<FragmentOffsetFrequency>>();
                foreach(var partition in _partitionSet)
                {
                    if(verbose) Console.WriteLine("{0}\t{1}\t{2}", partition.Charge, partition.SegmentIndex, partition.NeutralPeptideMass);
                    var numIons = reader.ReadInt32();
                    var fragmentOff = new FragmentOffsetFrequency[numIons];
                    for (var ionIndex = 0; ionIndex < numIons; ionIndex++)
                    {
                        var isPrefix = reader.ReadBoolean();
                        var charge = reader.ReadInt32();
                        var offset = reader.ReadSingle();
                        var frequency = reader.ReadSingle();
                    }
                }

                // TODO: determine fragment ions

                // Rank distributions
                var maxRank = reader.ReadInt32();
                for (var i = 0; i < numPartitions; i++)
                {
                    // repeat this for #Ions
                    for (var rank = 0; rank < maxRank + 1; rank++)
                    {
                        var frequency = reader.ReadSingle();
                    }
                }

                // Error distributions
                var errorScalingFactor = reader.ReadInt32();
                if (errorScalingFactor > 0)
                {
                    var numErrorBins = errorScalingFactor*2 + 1;

                    // Ion error distribution
                    for (var i = 0; i < numErrorBins; i++)
                    {
                        var ionErrorDist = reader.ReadSingle();
                    }

                    // Noise error distribution
                    for (var i = 0; i < numErrorBins; i++)
                    {
                        var noiseErrorDist = reader.ReadSingle();
                    }

                    // Ion existence table
                    for (var i = 0; i < 4; i++)
                    {
                        var ionExFreq = reader.ReadSingle();
                    }
                }

                // Validation
                var validationCode = reader.ReadInt32();
                if (validationCode != Int32.MaxValue)
                {
                    Console.WriteLine("Error while reading parameter file {0}", paramFile);
                    System.Environment.Exit(-1); // Error
                }
            }
        }

        private void DetermineIonTypes()
        {
            //_ionTypeTable = new Dictionary<int, IonType[]>();

            //for (Partition partition : partitionSet)
            //{
            //    ArrayList<FragmentOffsetFrequency> offList = fragOFFTable.get(partition);
            //    IonType[] ionTypes = new IonType[offList.size()];
            //    for (int i = 0; i < offList.size(); i++)
            //        ionTypes[i] = offList.get(i).getIonType();
            //    _ionTypeTable.put(partition, ionTypes);
            //}

            //mainIonTable = new HashMap<Partition, IonType>();
            //for (Partition partition :
            //partitionSet)
            //{
            //    if (partition.getSegNum() != 0)
            //        continue;
            //    HashMap<IonType, Float> ionProb = new HashMap<IonType, Float>();
            //    for (int seg = 0; seg < numSegments; seg++)
            //    {
            //        Partition part = new Partition(partition.getCharge(), partition.getParentMass(), seg);
            //        ArrayList<FragmentOffsetFrequency> offList = fragOFFTable.get(part);
            //        for (FragmentOffsetFrequency off :
            //        offList)
            //        {
            //            Float prob = ionProb.get(off.getIonType());
            //            if (prob == null)
            //                ionProb.put(off.getIonType(), off.getFrequency());
            //            else
            //                ionProb.put(off.getIonType(), prob + off.getFrequency());
            //        }
            //    }
            //    IonType mainIon = null;
            //    float prob = -1;
            //    for (IonType ion :
            //    ionProb.keySet())
            //    {
            //        if (ionProb.get(ion) > prob)
            //        {
            //            mainIon = ion;
            //            prob = ionProb.get(ion);
            //        }
            //    }
            //    assert(mainIon != null);
            //    for (int seg = 0; seg < numSegments; seg++)
            //    {
            //        Partition part = new Partition(partition.getCharge(), partition.getParentMass(), seg);
            //        mainIonTable.put(part, mainIon);
            //    }
            //}
        }
    }
}
