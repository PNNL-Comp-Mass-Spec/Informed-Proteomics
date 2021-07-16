using System;
using System.IO;
using System.Linq;
using System.Reflection;
using InformedProteomics.Backend.MassSpecData;
using NUnit.Framework;
using PRISM;

namespace InformedProteomics.Tests.Base
{
    public static class Utils
    {
        // Ignore Spelling: yyyy-MM-dd, hh:mm:ss tt, pbf

        public const string DEFAULT_TEST_FILE_FOLDER = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles";

        public const string DEFAULT_SPEC_FILES_FOLDER = DEFAULT_TEST_FILE_FOLDER + @"\SpecFiles";

        private const string MZML_TEST_FILE = "QC_Shew_Intact_26Sep14_Bane_C2Column3_Excerpt.mzML";

        private const string PBF_TEST_FILE = "QC_Shew_Intact_26Sep14_Bane_C2Column3_Excerpt.pbf";

        private const string DATE_TIME_FORMAT = "yyyy-MM-dd hh:mm:ss tt";

        private static string mLastStatus = string.Empty;

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
            {
                return mzmlFileInfo.FullName;
            }

            // The file is missing; return the default path, despite the fact that the file does not exist
            // The calling method will discover that the file is missing and act accordingly
            return mzmlFilePath;
        }

        /// <summary>
        /// Look for the default .pbf testing file in the default test file folder,
        /// plus also in the UnitTest_Files directory in this project or solution
        /// </summary>
        /// <param name="createIfMissing">If true and the .mzML file is found, create the .pbf file</param>
        /// <returns>Path to the file if found, otherwise the default path on Proto-2</returns>
        public static string GetPbfTestFilePath(bool createIfMissing)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;

            var pbfFilePath = Path.Combine(DEFAULT_SPEC_FILES_FOLDER, PBF_TEST_FILE);

            var pbfFileInfo = GetTestFile("PbfTestFilePath", pbfFilePath, false);

            if (pbfFileInfo != null)
            {
                // Check for a lock file, which would indicate another process is creating the .pbf file
                var existingLockFile = new FileInfo(Path.Combine(pbfFileInfo.FullName + ".lock"));
                WaitForLockFile(existingLockFile, false);

                pbfFileInfo.Refresh();
                if (pbfFileInfo.Exists)
                {
                    DeleteLockFile(existingLockFile);
                    return pbfFileInfo.FullName;
                }
            }

            if (!createIfMissing)
            {
                // Not creating the file if missing; return the default path, despite the fact that the file does not exist
                // The calling method will discover that the file is missing and act accordingly
                return pbfFilePath;
            }

            // Create the missing file

            var mzmlFilePath = Path.Combine(DEFAULT_SPEC_FILES_FOLDER, MZML_TEST_FILE);

            var mzmlFileInfo = GetTestFile("MzmlTestFilePath", mzmlFilePath, false);
            if (mzmlFileInfo?.DirectoryName == null)
            {
                // Unable to create the pbf file; return the default path, despite the fact that the file does not exist
                return pbfFilePath;
            }

            ShowMessage(methodName, string.Format("Creating {0} using {1}", PBF_TEST_FILE, mzmlFileInfo.FullName));

            var lockFile = new FileInfo(Path.Combine(mzmlFileInfo.DirectoryName, Path.GetFileNameWithoutExtension(mzmlFileInfo.Name) + ".pbf.lock"));
            var newPbfFilePath = Path.Combine(mzmlFileInfo.DirectoryName, PBF_TEST_FILE);

            try
            {
                // Create a new lock file so other processes know this thread is creating the .pbf file
                WaitForLockFile(lockFile, true);

                mLastStatus = string.Empty;
                var startTime = DateTime.UtcNow;

                var reader = MassSpecDataReaderFactory.GetMassSpecDataReader(mzmlFileInfo.FullName);
                var progress = new Progress<ProgressData>(p =>
                {
                    p.UpdateFrequencySeconds = 2;
                    if (p.Percent < 100 && (p.Percent % 25).Equals(0) || p.ShouldUpdate())
                    {
                        var statusMessage = string.Format("{0}, {1:00.0}% complete                        ", p.Status, p.Percent);

                        if (string.Equals(mLastStatus, statusMessage))
                        {
                            return;
                        }

                        mLastStatus = statusMessage;
                        Console.Write("\r{0}, {1:00.0}% complete                        ", p.Status, p.Percent);
                    }
                });

                var run = new PbfLcMsRun(mzmlFileInfo.FullName, reader, newPbfFilePath, 0, 0, progress);
                Console.WriteLine();
                ShowMessage(methodName, string.Format("Created {0} in {1:F0} seconds", run.PbfFilePath, DateTime.UtcNow.Subtract(startTime).TotalSeconds));

                DeleteLockFile(lockFile);

                return run.PbfFilePath;
            }
            catch (Exception ex)
            {
                ShowMessage(methodName, string.Format("Exception creating {0} using {1}: {2}", PBF_TEST_FILE, mzmlFileInfo.FullName, ex.Message));

                try
                {
                    var incompletePbfFile = new FileInfo(newPbfFilePath);
                    if (incompletePbfFile.Exists)
                    {
                        incompletePbfFile.Delete();
                    }

                    DeleteLockFile(lockFile);
                }
                catch
                {
                    // Ignore errors here
                }

                // Return the default path, despite the fact that the file does not exist
                return pbfFilePath;
            }
        }

        private static void CreateLockFile(FileSystemInfo lockFile)
        {
            using var writer = new StreamWriter(new FileStream(lockFile.FullName, FileMode.Create, FileAccess.Write, FileShare.ReadWrite));
            writer.WriteLine("Creating PBF file: " + DateTime.Now.ToString(DATE_TIME_FORMAT));
        }

        private static void DeleteLockFile(FileSystemInfo lockFile)
        {
            try
            {
                if (lockFile == null)
                {
                    return;
                }

                lockFile.Refresh();
                if (lockFile.Exists)
                {
                    lockFile.Delete();
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine("Error in DeleteLockFile: " + ex.Message);
            }
        }

        private static void WaitForLockFile(FileSystemInfo lockFile, bool createNewLockFile, int waitTimeMinutes = 5)
        {
            try
            {
                if (lockFile == null)
                {
                    return;
                }

                if (!lockFile.Exists)
                {
                    if (createNewLockFile)
                    {
                        CreateLockFile(lockFile);
                    }

                    return;
                }

                var timestamp = lockFile.LastWriteTimeUtc;
                while (DateTime.UtcNow.Subtract(timestamp).TotalMinutes < waitTimeMinutes)
                {
                    lockFile.Refresh();
                    if (!lockFile.Exists)
                    {
                        // The other process has deleted the lock file
                        break;
                    }

                    System.Threading.Thread.Sleep(5000);
                }

                if (createNewLockFile)
                {
                    CreateLockFile(lockFile);
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine("Error in WaitForLockFile: " + ex.Message);
            }
        }

        /// <summary>
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
                    {
                        break;
                    }

                    parentFolder = altFile2.Directory.Parent;
                    iteration++;
                }
            }

            if (assertIgnoreIfNotFound)
            {
                Assert.Ignore("Skipping test {0} since file not found: {1}", callingMethod, filePath);
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
            {
                ShowStarting(methodName + " ( ?? UnknownFilePath ?? )");
            }
            else if (filePathForTestCase.Any(x => Path.GetInvalidPathChars().Contains(x)))
            {
                ShowStarting(methodName + " (" + filePathForTestCase + ")");
            }
            else
            {
                ShowStarting(methodName + " (" + Path.GetFileName(filePathForTestCase) + ")");
            }
        }

        private static void ShowMessage(string methodName, string message)
        {
            Console.WriteLine(DateTime.Now.ToString(DATE_TIME_FORMAT) + ", {0}, {1}", message, methodName);
        }
    }
}