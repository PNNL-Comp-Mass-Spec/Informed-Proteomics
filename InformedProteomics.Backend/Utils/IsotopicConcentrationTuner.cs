using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MathAndStats;

namespace InformedProteomics.Backend.Utils
{
    /// <summary>
    /// Class attempts to find ideal concentration of particular isotope by comparing it to an
    /// observed isotopic profile
    /// </summary>
    public class IsotopicConcentrationTuner
    {
        // Ignore Spelling: averagine, monoisotope

        /// <summary>
        /// Initializes new instance of the <see cref="IsotopicConcentrationTuner" /> class
        /// </summary>
        public IsotopicConcentrationTuner()
        {
            // Set default parameters
            Element = Atom.Get("C");
            IsotopeIndex = 1;
            Tolerance = new Tolerance(10, ToleranceUnit.Ppm);
            ObservedPeaks = new List<Peak>();
            Mass = 0.0;
            Charge = 1;
            StepSize = 0.1;
            MaxConcentration = 20.0;
            RelativeIntensityThreshold = 0.01;
        }

        /// <summary>
        /// Gets or sets the element for which isotope proportions will be manipulated
        /// </summary>
        public Atom Element { get; set; }

        /// <summary>
        /// Gets or sets the index of the isotope to manipulate, relative to the monoisotope
        /// </summary>
        /// <remarks>Monoisotope is index 0</remarks>
        public int IsotopeIndex { get; set; }

        /// <summary>
        /// Gets the peak tolerance for matching observed peaks to theoretical peaks
        /// </summary>
        public Tolerance Tolerance { get; set; }

        /// <summary>
        /// Gets or sets the list of observed peaks to compare to the theoretical isotope profile
        /// </summary>
        public List<Peak> ObservedPeaks { get; set; }

        /// <summary>
        /// Gets or sets the monoisotopic mass of the ion for which isotope peaks will be calculated
        /// </summary>
        public double Mass { get; set; }

        /// <summary>
        /// Gets or sets the charge of the ion for which isotope peaks will be calculated
        /// </summary>
        public int Charge { get; set; }

        /// <summary>
        /// Gets or sets the amount to increase the concentration of the selected isotope index for each iteration
        /// </summary>
        public double StepSize { get; set; }

        /// <summary>
        /// Gets or sets the maximum concentration of the selected isotope to consider
        /// </summary>
        public double MaxConcentration { get; set; }

        /// <summary>
        /// Gets or sets the least abundant theoretical isotope peak to consider, relative to the highest theoretical isotope peak
        /// </summary>
        public double RelativeIntensityThreshold { get; set; }

        /// <summary>
        /// Try to find the best concentration of the selected isotope
        /// by stepping through the concentrations and fitting a theoretical
        /// isotopic profile to the provided observed peaks
        /// </summary>
        /// <param name="progress">The progress reporter</param>
        public IsotopeConcentrationCorrelationCurve Tune(IProgress<PRISM.ProgressData> progress = null)
        {
            // Set up progress reporter
            var progressData = new PRISM.ProgressData(progress);

            ValidateParameters();

            // Get default proportions for the selected element.
            // Copy it to a new array so we can manipulate it.
            var proportions = GetDefaultProportions(Element).ToArray();

            // Make sure this is an isotope we know about and that it isn't the monoisotope
            if (IsotopeIndex < 1 || IsotopeIndex >= proportions.Length)
            {
                throw new ArgumentOutOfRangeException(nameof(IsotopeIndex));
            }

            var defaultProportion = proportions[IsotopeIndex];

            // Set the default best point (the first one).
            var results = new IsotopeConcentrationCorrelationCurve
            {
                BestConcentration = new IsotopeConcentrationCorrelationCurve.ConcentrationCorrelationPoint
                {
                    IsotopeConcentration = defaultProportion,
                    MonoisotopeConcentration = proportions[0],
                    PearsonCorrelation = 0.0
                }
            };

            // Iterate over concentration values
            var numberOfSteps = (int)(MaxConcentration - defaultProportion / StepSize);
            for (var i = 0; i <= numberOfSteps; i++)
            {
                // Update percent complete
                progressData.Report(i, numberOfSteps);

                // Calculate concentration
                var concentrationStep = i * StepSize;
                proportions[IsotopeIndex] += concentrationStep; // Increase isotope of interest
                proportions[0] -= concentrationStep;                 // Decrease monoisotope

                // Get theoretical isotope profile and align the observed peaks to it
                var theoreticalIsotopeProfile = GetTheoreticalIsotopeProfile(proportions);
                var alignedObservedPeaks = AlignObservedPeaks(ObservedPeaks, theoreticalIsotopeProfile, Tolerance);

                // Break out the intensities of the isotope profiles
                var theoreticalIntensities = theoreticalIsotopeProfile.Select(peak => peak.Intensity).ToArray();
                var obsIntensities = alignedObservedPeaks.Select(peak => peak.Intensity).ToArray();

                // Compute Pearson correlation
                var pearsonCorrelation = FitScoreCalculator.GetPearsonCorrelation(obsIntensities, theoreticalIntensities);

                // Add data point for this concentration to result curve
                var dataPoint = new IsotopeConcentrationCorrelationCurve.ConcentrationCorrelationPoint
                {
                    IsotopeConcentration = proportions[IsotopeIndex],
                    MonoisotopeConcentration = proportions[0],
                    PearsonCorrelation = pearsonCorrelation
                };

                results.DataPoints.Add(dataPoint);

                // If this concentration has a better fit, update the stored results
                if (pearsonCorrelation >= results.BestConcentration.PearsonCorrelation)
                {
                    results.BestConcentration = dataPoint;
                }
            }

            return results;
        }

        /// <summary>
        /// Gets the theoretical isotope profile calculated using Averagine with the provided
        /// isotope proportions
        /// </summary>
        /// <param name="proportions">The proportions of each isotope</param>
        /// <returns>The theoretical isotope profile peaks</returns>
        public List<Peak> GetTheoreticalIsotopeProfile(double[] proportions)
        {
            // Get IsoProfilePredictor with updated proportions
            var isoProfilePredictor = new IsoProfilePredictor(
                Element.Code == "C" ? proportions : IsoProfilePredictor.DefaultProbC,
                Element.Code == "H" ? proportions : IsoProfilePredictor.DefaultProbH,
                Element.Code == "N" ? proportions : IsoProfilePredictor.DefaultProbN,
                Element.Code == "O" ? proportions : IsoProfilePredictor.DefaultProbO,
                Element.Code == "S" ? proportions : IsoProfilePredictor.DefaultProbS
            );

            var averagine = new Averagine();

            return averagine.GetTheoreticalIsotopeProfileInst(
                                Mass,
                                Charge,
                                RelativeIntensityThreshold,
                                isoProfilePredictor);
        }

        /// <summary>
        /// Aligns observed peak list to theoretical peak list
        /// </summary>
        /// <param name="observedPeaks"></param>
        /// <param name="theoreticalPeaks"></param>
        /// <param name="tolerance"></param>
        /// <returns>List of peaks</returns>
        public List<Peak> AlignObservedPeaks(IList<Peak> observedPeaks, IList<Peak> theoreticalPeaks, Tolerance tolerance = null)
        {
            tolerance ??= new Tolerance(10, ToleranceUnit.Ppm);

            // Remove empty peaks
            observedPeaks = observedPeaks.Where(peak => peak.Mz > 0.0).ToList();

            var alignedPeaks = new List<Peak> { Capacity = theoreticalPeaks.Count };

            var j = 0;
            foreach (var theoreticalPeak in theoreticalPeaks)
            {
                var tolDa = tolerance.GetToleranceAsMz(theoreticalPeak.Mz);
                var maxMz = theoreticalPeak.Mz + tolDa;
                var obsPeak = observedPeaks[j];

                var selectedPeak = new Peak(theoreticalPeak.Mz, 0);
                while (obsPeak.Mz <= maxMz)
                {
                    var diff = Math.Abs(obsPeak.Mz - theoreticalPeak.Mz);
                    if (diff < tolDa && obsPeak.Intensity > selectedPeak.Intensity)
                    {
                        selectedPeak = obsPeak;
                    }

                    j = Math.Min(observedPeaks.Count - 1, j + 1); // Increment, but do not go out of bounds
                    obsPeak = observedPeaks[j];
                }

                alignedPeaks.Add(selectedPeak);
            }

            return alignedPeaks;
        }

        /// <summary>
        /// Get the default isotope proportions for the given element
        /// </summary>
        /// <param name="element">The element to get the default proportions for</param>
        /// <returns>An array where each index is the </returns>
        private double[] GetDefaultProportions(Atom element)
        {
            double[] proportions = null;

            if (element.Code == "C")
            {
                proportions = IsoProfilePredictor.DefaultProbC;
            }
            else if (element.Code == "H")
            {
                proportions = IsoProfilePredictor.DefaultProbH;
            }
            else if (element.Code == "N")
            {
                proportions = IsoProfilePredictor.DefaultProbN;
            }
            else if (element.Code == "O")
            {
                proportions = IsoProfilePredictor.DefaultProbO;
            }
            else if (element.Code == "S")
            {
                proportions = IsoProfilePredictor.DefaultProbS;
            }

            return proportions;
        }

        /// <summary>
        /// Checks to make sure that the selected element is one that is possible to manipulate
        /// </summary>
        private void ValidateParameters()
        {
            // Make sure that the selected element is an element that we can manipulate
            if (Element.Code != "C" && Element.Code != "H" && Element.Code != "N" && Element.Code != "O" && Element.Code != "S")
            {
                throw new ArgumentException(string.Format("Cannot manipulate isotope proportions for {0}.", Element.Name));
            }
        }

        /// <summary>
        /// Class representing a the results of the isotope concentration tuning
        /// </summary>
        public class IsotopeConcentrationCorrelationCurve
        {
            /// <summary>
            /// Initializes a new instance of the <see cref="IsotopeConcentrationCorrelationCurve" />
            /// </summary>
            public IsotopeConcentrationCorrelationCurve()
            {
                DataPoints = new List<ConcentrationCorrelationPoint>();
                BestConcentration = new ConcentrationCorrelationPoint();
            }

            /// <summary>
            /// Gets or sets the curve showing isotope concentration vs Pearson correlation with fit to observed peaks
            /// </summary>
            public List<ConcentrationCorrelationPoint> DataPoints { get; set; }

            /// <summary>
            /// The concentration with the best fit with observed peaks
            /// </summary>
            public ConcentrationCorrelationPoint BestConcentration { get; set; }

            /// <summary>
            /// Class representing a single point in the curve showing concentration vs correlation
            /// </summary>
            public class ConcentrationCorrelationPoint
            {
                /// <summary>
                /// Gets the concentration of selected isotope
                /// </summary>
                public double IsotopeConcentration { get; set; }

                /// <summary>
                /// Gets the concentration of the monoisotope
                /// </summary>
                public double MonoisotopeConcentration { get; set; }

                /// <summary>
                /// Gets Pearson correlation of the isotope concentration fit to the observed peaks
                /// </summary>
                public double PearsonCorrelation { get; set; }
            }
        }
    }
}
