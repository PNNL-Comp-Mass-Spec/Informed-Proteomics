using System;
using System.Collections.Generic;
using System.Linq;

using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;

using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.Utils
{
    /// <summary>
    /// Class attempts to find ideal concentration of particular isotope by comparing it to an
    /// observed isotopic profile.
    /// </summary>
    public class IsotopicConcentrationTuner
    {
        /// <summary>
        /// Initializes new instance of the <see cref="IsotopicConcentrationTuner" /> class.
        /// </summary>
        public IsotopicConcentrationTuner()
        {
            // Set default parameters
            this.Element = Atom.Get("C");
            this.IsotopeIndex = 1;
            this.Tolerance = new Tolerance(10, ToleranceUnit.Ppm);
            this.ObservedPeaks = new List<Peak>();
            this.Mass = 0.0;
            this.Charge = 1;
            this.StepSize = 0.1;
            this.MaxConcentration = 20.0;
            this.RelativeIntensityThreshold = 0.01;
        }

        /// <summary>
        /// Gets or sets the element to manipulte isotope proportions for.
        /// </summary>
        public Atom Element { get; set; }

        /// <summary>
        /// Gets or sets the index of the isotope to manipulate, relative to the monoisotope.
        /// </summary>
        /// <remarks>Monoisotope is index 0.</remarks>
        public int IsotopeIndex { get; set; }

        /// <summary>
        /// Gets the peak tolerance for matching observed peaks to theoretical peaks.
        /// </summary>
        public Tolerance Tolerance { get; set; }

        /// <summary>
        /// Gets or sets the list of observed peaks to compare to the theoretical isotope profile.
        /// </summary>
        public List<Peak> ObservedPeaks { get; set; }

        /// <summary>
        /// Gets or sets the monoisotopic mass of the ion to calculate isotope peaks for.
        /// </summary>
        public double Mass { get; set; }

        /// <summary>
        /// Gets or sets the charge of the ion to calculate isotope peaks for.
        /// </summary>
        public int Charge { get; set; }

        /// <summary>
        /// Gets or sets the amount to increase the concentration of the selected isotope index for each iteration.
        /// </summary>
        public double StepSize { get; set; }

        /// <summary>
        /// Gets or sets the maximum concentration of the selected isotope to consider.
        /// </summary>
        public double MaxConcentration { get; set; }

        /// <summary>
        /// Gets or sets the least abundant theoretical isotope peak to consider, relative to the highest theoretical isotope peak.
        /// </summary>
        public double RelativeIntensityThreshold { get; set; }

        /// <summary>
        /// Try to find the best concentration of the selected isotope
        /// by stepping through the concentrations and fitting a theoretical 
        /// isotopic profile to the provided observed peaks.
        /// </summary>
        /// <param name="progress">The progress reporter.</param>
        public IsotopeConcentrationCorrelationCurve Tune(IProgress<ProgressData> progress = null)
        {
            // Set up progress reporter
            progress = progress ?? new Progress<ProgressData>();
            var progressData = new ProgressData();

            this.ValidateParameters();

            // Get default proportions for the selected element.
            // Copy it to a new array so we can manipulate it.
            var proportions = this.GetDefaultProportions(this.Element).ToArray();

            // Make sure this is an isotope we know about and that it isn't the monoisotope
            if (this.IsotopeIndex < 1 || this.IsotopeIndex >= proportions.Length)
            {
                throw new ArgumentOutOfRangeException("isotopeIndex");
            }

            var defaultProportion = proportions[this.IsotopeIndex];

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
            int numberOfSteps = (int)(this.MaxConcentration - defaultProportion / this.StepSize);
            for (int i = 0; i <= numberOfSteps; i++)
            {
                // Update percent complete
                progress.Report(progressData.UpdatePercent(100.0 * i / numberOfSteps));

                // Calculate concentration
                var concentrationStep = i * this.StepSize;
                proportions[this.IsotopeIndex] += concentrationStep; // Increase isotope of interest
                proportions[0] -= concentrationStep;                 // Decrease monoisotope

                // Get theoretical isotope profile and align the observed peaks to it
                var theoreticalIsotopeProfile = this.GetTheoreticalIsotopeProfile(proportions);
                var alignedObservedPeaks = this.AlignObservedPeaks(this.ObservedPeaks, theoreticalIsotopeProfile, this.Tolerance);

                // Break out the intensities of the isotope profiles
                var theoIntensities = theoreticalIsotopeProfile.Select(peak => peak.Intensity).ToArray();
                var obsIntensities = alignedObservedPeaks.Select(peak => peak.Intensity).ToArray();

                // Compute pearson correlation
                var pearsonCorrelation = FitScoreCalculator.GetPearsonCorrelation(obsIntensities, theoIntensities);

                // Add data point for this concentration to result curve
                var dataPoint = new IsotopeConcentrationCorrelationCurve.ConcentrationCorrelationPoint
                {
                    IsotopeConcentration = proportions[this.IsotopeIndex],
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
        /// isotope proportions.
        /// </summary>
        /// <param name="proportions">The proportions of each isotope.</param>
        /// <returns>The theoretical isotope profile peaks.</returns>
        public List<Peak> GetTheoreticalIsotopeProfile(double[] proportions)
        {
            // Get IsoProfilePredictor with updated proportions
            var isoProfilePredictor = new IsoProfilePredictor(
                this.Element.Code == "C" ? proportions : IsoProfilePredictor.DefaultProbC,
                this.Element.Code == "H" ? proportions : IsoProfilePredictor.DefaultProbH,
                this.Element.Code == "N" ? proportions : IsoProfilePredictor.DefaultProbN,
                this.Element.Code == "O" ? proportions : IsoProfilePredictor.DefaultProbO,
                this.Element.Code == "S" ? proportions : IsoProfilePredictor.DefaultProbS
            );

            return Averagine.GetTheoreticalIsotopeProfile(
                                this.Mass, 
                                this.Charge,
                                this.RelativeIntensityThreshold,
                                isoProfilePredictor);
        }

        /// <summary>
        /// Aligns observed peak list to theoretical peak list.
        /// </summary>
        /// <param name="observedPeaks"></param>
        /// <param name="theoreticalPeaks"></param>
        /// <param name="tolerance"></param>
        /// <returns></returns>
        public List<Peak> AlignObservedPeaks(IList<Peak> observedPeaks, IList<Peak> theoreticalPeaks, Tolerance tolerance = null)
        {
            tolerance = tolerance ?? new Tolerance(10, ToleranceUnit.Ppm);

            // Remove empty peaks
            observedPeaks = observedPeaks.Where(peak => peak.Mz > 0.0).ToList();

            List<Peak> alignedPeaks = new List<Peak> { Capacity = theoreticalPeaks.Count };

            int j = 0;
            foreach (var theoPeak in theoreticalPeaks)
            {
                double tolDa = tolerance.GetToleranceAsTh(theoPeak.Mz);
                double maxMz = theoPeak.Mz + tolDa;
                Peak obsPeak = observedPeaks[j];

                Peak selectedPeak = new Peak(theoPeak.Mz, 0);
                while (obsPeak.Mz <= maxMz)
                {
                    double diff = Math.Abs(obsPeak.Mz - theoPeak.Mz);
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
        /// Get the default isotope proportions for the given element.
        /// </summary>
        /// <param name="element">The element to get the default proportions for.</param>
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
        /// Checks to make sure that the selected element is one that is possible to manipulate.
        /// </summary>
        private void ValidateParameters()
        {
            // Make sure that the selected element is an element that we can manipulate
            if (this.Element.Code != "C" && this.Element.Code != "H" && this.Element.Code != "N" && this.Element.Code != "O" && this.Element.Code != "S")
            {
                throw new ArgumentException(string.Format("Cannot manipulate isotope proportions for {0}.", this.Element.Name));
            }
        }

        /// <summary>
        /// Class representing a the results of the isotope concentration tuning.
        /// </summary>
        public class IsotopeConcentrationCorrelationCurve
        {
            /// <summary>
            /// Gets or sets the curve showing isotope concentration vs pearson correlation with fit to observed peaks.
            /// </summary>
            public List<ConcentrationCorrelationPoint> DataPoints { get; set; }
            
            /// <summary>
            /// The concentration with the best fit with observed peaks.
            /// </summary>
            public ConcentrationCorrelationPoint BestConcentration { get; set; }

            /// <summary>
            /// Class representing a single point in the curve showing concentration vs correlation.
            /// </summary>
            public class ConcentrationCorrelationPoint
            {
                /// <summary>
                /// Gets the concentration of selected isotope.
                /// </summary>
                public double IsotopeConcentration { get; set; }

                /// <summary>
                /// Gets the concentration of monoisotope..
                /// </summary>
                public double MonoisotopeConcentration { get; set; }

                /// <summary>
                /// Gets pearson correlation of the isotope concentration fit to the observed peaks.
                /// </summary>
                public double PearsonCorrelation { get; set; }
            }
        }
    }
}
