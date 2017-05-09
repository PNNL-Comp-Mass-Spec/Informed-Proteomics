using System;
using System.IO;
using System.Linq;

namespace InformedProteomics.Tests.Base
{
    public static class Utils
    {
        public const string DEFAULT_TEST_FILE_FOLDER = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles";

        public const string DEFAULT_SPEC_FILES_FOLDER = DEFAULT_TEST_FILE_FOLDER + @"\SpecFiles";

        private const string MZML_TEST_FILE = @"QC_Shew_Intact_26Sep14_Bane_C2Column3_Excerpt.mzML";

        private const string PBF_TEST_FILE = @"QC_Shew_Intact_26Sep14_Bane_C2Column3_Excerpt.pbf";

        /// <summary>
        /// Look for the default .mzML testing file in the default test file folder,
        /// plus also in the UnitTest_Files directory in this project or solution
        /// </summary>
        /// <returns>Path to the file if found, otherwise the default path on Proto-2</returns>
        public static string GetMzmlTestFilePath()
        {
            var mzmlFilePath = Path.Combine(DEFAULT_SPEC_FILES_FOLDER, MZML_TEST_FILE);

            var mzmlFileInfo = GetTestFile("MzmlTestFilePath", mzmlFilePath, false);

            if (mzmlFileInfo != null)
                return mzmlFileInfo.FullName;

            // The file is missing; return the default path, despite the fact that the file does not exist
            // The calling method will discover that the file is missing and act accordingly
            return mzmlFilePath;
        }

        /// <summary>
        /// Look for the default .pbf testing file in the default test file folder,
        /// plus also in the UnitTest_Files directory in this project or solution
        /// </summary>
        /// <param name="createIfMissing">If true and the .mzML file is ofund, create the .pbf file</param>
        /// <returns>Path to the file if found, otherwise the default path on Proto-2</returns>
        public static string GetPbfTestFilePath(bool createIfMissing)
        {
            var pbfFilePath = Path.Combine(DEFAULT_SPEC_FILES_FOLDER, PBF_TEST_FILE);

            var pbfFileInfo = GetTestFile("PbfTestFilePath", pbfFilePath, false);
            if (pbfFileInfo != null)
            {
                    return pbfFileInfo.FullName;
            }

                // Return the default path, despite the fact that the file does not exist
                return pbfFilePath;
        }
        /// Checks for the existence of filePath
        /// If not found, but the file exists in the UnitTest_Files directory in this project or solution, returns the local file path
        /// If not found in either place, raises an IgnoreException by invoking Assert.Ignore
        /// </summary>
        /// <param name="callingMethod">Calling method name</param>
        /// <param name="filePath">File path to check for, typically at \\proto-2\UnitTest_Files\InformedProteomics_TestFiles</param>
        /// <param name="assertIgnoreIfNotFound">When true, invoke Assert.Ignore if the file is not found</param>
        /// <returns>File path if found, otherwise an empty string</returns>
        /// <remarks>If the file is not found, Assert.Ignore is called, which will throw and IgnoreException</remarks>
        public static FileInfo GetTestFile(string callingMethod, string filePath, bool assertIgnoreIfNotFound = true)
        {
            var testFile = new FileInfo(filePath);
            if (testFile.Exists)
            {
                return testFile;
            }

            var altFile = new FileInfo(Path.Combine("UnitTest_Files", Path.GetFileName(filePath)));
            if (altFile.Exists)
            {
                return altFile;
            }

            // File still not found; try parent directories, up to 4 levels up
            var iteration = 0;
            if (altFile.Directory != null)
            {
                var parentFolder = altFile.Directory.Parent;

                while (iteration < 4 && parentFolder != null)
                {
                    var altFile2 = new FileInfo(Path.Combine(parentFolder.FullName, altFile.Name));
                    if (altFile2.Exists)
                    {
                        return altFile2;
                    }

                    if (altFile2.Directory == null)
                        break;

                    parentFolder = altFile2.Directory.Parent;
                    iteration++;
                }
            }

            if (assertIgnoreIfNotFound)
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", callingMethod, filePath);
            }

            return null;

        }

        public static void ShowEnding(string methodName)
        {
            ShowMessage(methodName, "Ending");
        }

        public static void ShowStarting(string methodName)
        {
            ShowMessage(methodName, "Starting");
        }

        public static void ShowStarting(string methodName, string filePathForTestCase)
        {
            if (string.IsNullOrWhiteSpace(filePathForTestCase))
                ShowStarting(methodName + " ( ?? UnknownFilePath ?? )");
            else if (filePathForTestCase.Any(x => Path.GetInvalidPathChars().Contains(x)))
                ShowStarting(methodName + " (" + filePathForTestCase + ")");
            else
                ShowStarting(methodName + " (" + Path.GetFileName(filePathForTestCase) + ")");
        }

        private static void ShowMessage(string methodName, string message)
        {
            Console.WriteLine(DateTime.Now.ToString("yyyy-MM-dd hh:mm:ss") + @", {0}, {1}", message, methodName);
        }
    }
}