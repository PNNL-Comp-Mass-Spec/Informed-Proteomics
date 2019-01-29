using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.TopDown.Scoring.FlipScoring
{
    /// <summary>
    /// The scoring parameters for scoring ions in a single mass bin range.
    /// </summary>
    public class ScoringParameters
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="ScoringParameters" /> class.
        /// </summary>
        /// <param name="mass">The mass</param>
        /// <param name="tolerance"></param>
        /// <param name="precursorPeakFilter"></param>
        public ScoringParameters(double mass, Tolerance tolerance, PrecursorPeakFilter precursorPeakFilter)
        {
            this.Mass = mass;
            this.PrecursorPeakFilter = new PrecursorPeakFilter();
        }

        /// <summary>
        /// Gets the peak tolerance value used during training.
        /// </summary>
        public Tolerance TrainerTolerance { get; private set; }

        /// <summary>
        /// Gets the maximum bin mass for the scoring parameters.
        /// </summary>
        public double Mass { get; }

        /// <summary>
        /// Gets the trained precursor peak filter for removing precursor ions selected for this bin.
        /// </summary>
        public PrecursorPeakFilter PrecursorPeakFilter { get; }

        /// <summary>
        /// Gets the weights for the selected ion features.
        /// </summary>
        public FeatureWeights FeatureWeights { get; private set; }

        /// <summary>
        /// Gets the the ions selected for scoring for this bin.
        /// </summary>
        public BaseIonType[] SelectedIonTypes => this.FeatureWeights.IonWeights.Keys.ToArray();

        /// <summary>
        /// Parse a single scoring param file.
        /// </summary>
        /// <param name="filePath">The path to the scoring param file.</param>
        /// <returns>The scoring parameters sorted by max bin mass.</returns>
        public static ScoringParameters ParseFile(string filePath)
        {
            if (!File.Exists(filePath))
            {
                throw new FileNotFoundException(string.Format("Cannot find scoring param file: {0}", filePath));
            }

            throw new NotImplementedException();
        }

        /// <summary>
        /// Parse either a single scoring param file or an entire directory of scoring param files.
        /// </summary>
        /// <param name="path">The path to the scoring param directory or file.</param>
        /// <returns>The scoring parameters sorted by max bin mass.</returns>
        public static ScoringParameters[] Parse(string path)
        {
            if (File.Exists(path))
            {   // It's actually a single file path.
                return new[] { ParseFile(path) };
            }

            if (!Directory.Exists(path))
            {
                throw new DirectoryNotFoundException(string.Format("Cannot find scoring param source directory: {0}", path));
            }

            var files = Directory.GetFiles(path);
            var scoringParams = new List<ScoringParameters>();
            foreach (var file in files)
            {
                // See if this is a valid scoring param file
                var fileName = Path.GetFileNameWithoutExtension(path);
                var parts = fileName.Split('_');
                if (parts.Length < 2)
                {   // Not a valid file name for a scoring param file. Just skip it
                    continue;
                }

                // Parse the file
                scoringParams.Add(ParseFile(file));
            }

            // Sort by max bin mass and convert to an array.
            return scoringParams.OrderBy(x => x.Mass).ToArray();
        }
    }
}
