using System;
using System.IO;
using System.Linq;

namespace InformedProteomics.Tests.Base
{
    public static class Utils
    {
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
