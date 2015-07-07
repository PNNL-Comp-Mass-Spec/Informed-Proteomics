﻿using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using InformedProteomics.Backend.Data.Spectrometry;
using pwiz.CLI.analysis;
using pwiz.CLI.cv;
using pwiz.CLI.msdata;
using Spectrum = InformedProteomics.Backend.Data.Spectrometry.Spectrum;

namespace InformedProteomics.Backend.MassSpecData
{
    /// <summary>
    /// This class currently doesn't work because of an unresolved dll dependency.
    /// </summary>
    public sealed class ProteoWizardReader: IMassSpecDataReader, IDisposable
    {
        public static Assembly ProteoWizardAssemblyResolver(object sender, ResolveEventArgs args)
        {
            // https://support.microsoft.com/en-us/kb/837908
            //This handler is called only when the common language runtime tries to bind to the assembly and fails.
            string pwizPath = Environment.GetEnvironmentVariable("ProteoWizard");
            if (string.IsNullOrWhiteSpace(pwizPath) && Directory.Exists(@"C:\DMS_Programs\ProteoWizard"))
            {
                pwizPath = @"C:\DMS_Programs\ProteoWizard";
            }
            if (string.IsNullOrWhiteSpace(pwizPath))
            {
                var progFiles = Environment.GetEnvironmentVariable("ProgramFiles");
                if (string.IsNullOrWhiteSpace(progFiles))
                {
                    return null;
                }
                var progPwiz = Path.Combine(progFiles, "ProteoWizard");
                if (!Directory.Exists(progPwiz))
                {
                    return null;
                }
                var posPaths = Directory.GetDirectories(progPwiz, "ProteoWizard *");
                pwizPath = posPaths.Max(); // Try to get the "newest" folder
            }

            //Retrieve the list of referenced assemblies in an array of AssemblyName.
            string strTempAssmbPath = "";

            AssemblyName[] arrReferencedAssmbNames = Assembly.GetExecutingAssembly().GetReferencedAssemblies();

            //Loop through the array of referenced assembly names.
            foreach (AssemblyName strAssmbName in arrReferencedAssmbNames)
            {
                //Check for the assembly names that have raised the "AssemblyResolve" event.
                if (strAssmbName.FullName.Substring(0, strAssmbName.FullName.IndexOf(",")) == args.Name.Substring(0, args.Name.IndexOf(",")))
                {
                    //Build the path of the assembly from where it has to be loaded.                
                    strTempAssmbPath = Path.Combine(pwizPath, args.Name.Substring(0, args.Name.IndexOf(",")) + ".dll");
                    break;
                }
            }

            //Load the assembly from the specified path.                    
            var myAssembly = Assembly.LoadFrom(strTempAssmbPath);

            //Return the loaded assembly.
            return myAssembly;
        }

        public ProteoWizardReader(string filePath)
        {
            AppDomain.CurrentDomain.AssemblyResolve += ProteoWizardAssemblyResolver;

            var readers = ReaderList.FullReaderList;
            readers.read(filePath, _dataFile);
            if ((new string[] {".mzml", ".mzml.gz", ".mzxml", ".mzxml.gz", ".mgf", ".mgf.gz", ".txt"})
                .Any(ext => filePath.ToLower().EndsWith(ext)))
            {
                Console.WriteLine("Using cwt Centroiding");
                _filters.Add(_cwtCentroiding);
            }
            else
            {
                Console.WriteLine("Using vendor Centroiding");
                _filters.Add(_vendorCentroiding);
            }
            SpectrumListFactory.wrap(_dataFile, _filters);
        }

        public IEnumerable<Spectrum> ReadAllSpectra()
        {
            for (int i = 1; i <= NumSpectra; i++)
            {
                yield return ReadMassSpectrum(i);
            }
        }

        public bool TryMakeRandomAccessCapable()
        {
            return true;
        }

        public Spectrum ReadMassSpectrum(int scanIndex)
        {
            var pwizSpec = _dataFile.run.spectrumList.spectrum(scanIndex - 1, true);

            var msLevel = (int)(pwizSpec.cvParam(CVID.MS_ms_level).value);
            double[] mzArray = new double[0];
            double[] intensityArray = new double[0];
            foreach (var bda in pwizSpec.binaryDataArrays)
            {
                if (bda.hasCVParam(CVID.MS_m_z_array))
                {
                    mzArray = bda.data.ToArray();
                }
                if (bda.hasCVParam(CVID.MS_intensity_array))
                {
                    intensityArray = bda.data.ToArray();
                }
            }
            double scanTime = 0;
            foreach (var s in pwizSpec.scanList.scans)
            {
                if (s.hasCVParam(CVID.MS_scan_start_time))
                {
                    scanTime = (double) (s.cvParam(CVID.MS_scan_start_time).value);
                }
            }
            if (msLevel > 1)
            {
                double? thermoMonoMass = null;
                foreach (var up in pwizSpec.userParams)
                {
                    if (up.name == "[Thermo Trailer Extra]Monoisotopic M/Z:")
                    {
                        thermoMonoMass = (double) (up.value);
                    }
                }
                ActivationMethod am = ActivationMethod.Unknown;
                Data.Spectrometry.IsolationWindow iw = null;
                foreach (var precursor in pwizSpec.precursors)
                {
                    var act = precursor.activation;
                    foreach (var param in act.cvParams)
                    {
                        switch (param.cvid)
                        {
                            case CVID.MS_collision_induced_dissociation:
                                am = ActivationMethod.CID;
                                break;
                            case CVID.MS_electron_transfer_dissociation:
                                am = ActivationMethod.ETD;
                                break;
                            case CVID.MS_beam_type_collision_induced_dissociation:
                                am = ActivationMethod.HCD;
                                break;
                            case CVID.MS_electron_capture_dissociation:
                                am = ActivationMethod.ECD;
                                break;
                            case CVID.MS_pulsed_q_dissociation:
                                am = ActivationMethod.PQD;
                                break;
                        }
                    }
                    var piw = precursor.isolationWindow;
                    var target = (double)(piw.cvParam(CVID.MS_isolation_window_target_m_z).value);
                    var lowOff = (double)(piw.cvParam(CVID.MS_isolation_window_lower_offset).value);
                    var uppOff = (double)(piw.cvParam(CVID.MS_isolation_window_upper_offset).value);
                    int? charge = null;
                    double selectedIonMz = 0;
                    foreach (var si in precursor.selectedIons)
                    {
                        if (si.hasCVParam(CVID.MS_charge_state))
                        {
                            charge = (int)(si.cvParam(CVID.MS_charge_state).value);
                        }
                        selectedIonMz = (double) (si.cvParam(CVID.MS_selected_ion_m_z).value);
                    }
                    iw = new Data.Spectrometry.IsolationWindow(target, lowOff, uppOff, thermoMonoMass, charge);
                }
                return new ProductSpectrum(mzArray, intensityArray, scanIndex)
                {
                    NativeId = pwizSpec.id,
                    ActivationMethod = am,
                    IsolationWindow = iw,
                    MsLevel = msLevel,
                    ElutionTime = scanTime,
                };
            }
            return new Spectrum(mzArray, intensityArray, scanIndex)
            {
                NativeId = pwizSpec.id,
                ElutionTime = scanTime,
            };
        }

        public void Close()
        {
            _dataFile.Dispose();
        }

        public int NumSpectra
        {
            get { return _dataFile.run.spectrumList.size(); }
        }

        private readonly MSData _dataFile = new MSData();
        private readonly string _vendorCentroiding = "peakPicking true 1-";
        private readonly string _cwtCentroiding = "peakPicking cwt snr=1.0 peakSpace=0.1 msLevel=1-";
        private readonly List<string> _filters = new List<string>();

        public void Dispose()
        {
            _dataFile.Dispose();
        }
    }
}
