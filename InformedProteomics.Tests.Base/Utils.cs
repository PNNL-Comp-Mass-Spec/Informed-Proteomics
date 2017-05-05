using System;
using System.IO;
using System.Linq;

namespace InformedProteomics.Tests.Base
{
    public static class Utils
    {
        /// <summary>
        /// Checks for the existence of filePath
        /// If not found, but the file exists in the UnitTest_Files directory in this project or solution, returns the local file path
        /// If not found in either place, raises and Ignore
        /// </summary>
        /// <param name="callingMethod">Calling method name</param>
        /// <param name="filePath">File path to check for, typically at \\proto-2\UnitTest_Files\InformedProteomics_TestFiles</param>
        /// <returns>File path if found, otherwise an empty string</returns>
        /// <remarks>If the file is not found, Assert.Ignore is called, which will throw and IgnoreException</remarks>
        public static FileInfo GetTestFile(string callingMethod, string filePath)
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

            NUnit.Framework.Assert.Ignore(@"Skipping test {0} since file not found: {1}", callingMethod, filePath);
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