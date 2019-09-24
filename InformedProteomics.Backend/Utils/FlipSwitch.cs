namespace InformedProteomics.Backend.Utils
{
    /// <summary>
    /// Singleton class for conditionalizing use of FLIP scoring code.
    /// </summary>
    public class FlipSwitch
    {
        /// <summary>
        /// Singleton instance
        /// </summary>
        public static FlipSwitch Instance { get; }

        /// <summary>
        /// If true, FLIP scoring code should be used.
        /// </summary>
        public static bool UseFlipScoring => Instance.useFlipScoring;

        static FlipSwitch()
        {
            Instance = new FlipSwitch();
        }

        private FlipSwitch()
        {
            useFlipScoring = false;
        }

        private bool useFlipScoring;

        /// <summary>
        /// Set whether FLIP Scoring should be used or not.
        /// </summary>
        /// <param name="useIt">True if FLIP Scoring should be used</param>
        public void SetUseFlipScoring(bool useIt = false)
        {
            useFlipScoring = useIt;
        }
    }
}
