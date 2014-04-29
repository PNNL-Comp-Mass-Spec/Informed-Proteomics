using System;
using System.Collections.Generic;

namespace InformedProteomics.Backend.Utils
{
    /// <summary>
    /// Parses a list of strings for command line processing.
    /// Originally developed by Brian LaMarche and modified by Sangtae Kim
    /// </summary>
    public class CommandLineParser
    {
        /// <summary>
        /// Processes the command line args into a dictionary.  Skips any non-option items at beggnining of args list.
        /// </summary>
        /// <param name="args">Args to process</param>
        /// <param name="offset">Starting index of argument array.</param>
        /// <returns>Dictionary of arguments to process.</returns>
        public static IDictionary<string, List<string>> ProcessArgs(string[] args, int offset)
        {
            IDictionary<string, List<string>> options = new Dictionary<string, List<string>>();

            var i = offset;
            while (i < args.Length)
            {
                var key = args[i].ToLower();
                if (key.StartsWith("-"))
                {
                    if (options.ContainsKey(key))
                    {
                        throw new Exception("The argument key already exists " + key);
                    }
                    options.Add(key, new List<string>());
                    while (++i < args.Length)
                    {
                        if (!args[i].StartsWith("-"))
                        {
                            options[key].Add(args[i].ToLower());
                        }
                        else
                        {
                            break;
                        }
                    }
                }
                else
                {
                    i++;
                }
            }
            return options;
        }
    }
}