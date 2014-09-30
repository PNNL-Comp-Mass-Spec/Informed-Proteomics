﻿using System;
using System.Collections.Generic;
using System.IO;
using System.IO.Compression;
using System.Xml;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Backend.MassSpecData
{
    public sealed class MzMLReader: IMassSpecDataReader
    {
        private XmlReader _reader;

        /// <summary>
        /// Enumeration of common mzML versions
        /// </summary>
        private enum MzML_Version
        {
            mzML1_0_0,
            mzML1_1_0
        }

        /// <summary>
        /// Store the mzML version, so that we can use it to adjust how some things are processed.
        /// </summary>
        private MzML_Version _version;

        private class ScanData
        {
            public double MonoisotopicMz;
            public double StartTime;

            public ScanData()
            {
                MonoisotopicMz = 0.0;
                StartTime = 0.0;
            }
        }

        private enum ParamType
        {
            cvParam,
            userParam
        }

        private abstract class Param
        {
            public ParamType ParamType { get; protected set; }
            public string Name;          // Required
            public string Value;         // Optional
            public string UnitCVRef;     // Optional
            public string UnitAccession; // Optional
            public string UnitName;      // Optional

            protected string _cvRef;
            protected string _accession;
            protected string _type;

            public virtual string CVRef
            {
                get { return string.Empty; }
                set { _cvRef = string.Empty; }
            }

            public virtual string Accession
            {
                get { return string.Empty; }
                set { _accession = string.Empty; }
            }

            public virtual string Type
            {
                get { return string.Empty; }
                set { _type = string.Empty; }
            }

            protected Param()
            {
                Name = string.Empty;
                Value = string.Empty;
                UnitAccession = string.Empty;
                UnitCVRef = string.Empty;
                UnitName = string.Empty;
                _cvRef = string.Empty;
                _accession = string.Empty;
                _type = string.Empty;
            }
        }
        private class CVParam : Param
        {
            public override string CVRef      // Required
            {
                get { return _cvRef; }
                set { _cvRef = value; }
            }

            public override string Accession  // Required
            {
                get { return _accession; }
                set { _accession = value; }
            }

            public CVParam()
            {
                ParamType = ParamType.cvParam;
            }
        }

        private class UserParam : Param
        {
            public override string Type       // Optional
            {
                get { return _type; }
                set { _type = value; }
            }

            public UserParam()
            {
                ParamType = ParamType.userParam;
            }
        }

        private readonly Dictionary<string, List<Param>> _referenceableParamGroups = new Dictionary<string, List<Param>>();

        private class SelectedIon
        {
            public double SelectedIonMz;
            public int Charge;

            public SelectedIon()
            {
                SelectedIonMz = 0.0;
                Charge = 0;
            }
        }

        private class Precursor
        {
            public List<SelectedIon> Ions;
            public double IsolationWindowTargetMz;
            public double IsolationWindowLowerOffset;
            public double IsolationWindowUpperOffset;
            public ActivationMethod Activation;

            public Precursor()
            {
                Ions = new List<SelectedIon>();
                IsolationWindowTargetMz = 0.0;
                IsolationWindowLowerOffset = 0.0;
                IsolationWindowUpperOffset = 0.0;
                Activation = ActivationMethod.Unknown;
            }
        }

        private enum Precision
        {
            Precision32,
            Precision64
        };

        private enum ArrayType
        {
            m_z_array,
            intensity_array,
            charge_array,
            signal_to_noise_array,
            time_array,
            wavelength_array,
            non_standard_data_array,
            flow_rate_array,
            pressure_array,
            temperature_array
        };

        private enum Instrument
        {
            ABI_WIFF_format, //MS_ABI_WIFF_format, "MS:1000562", "ABI WIFF format", "Applied Biosystems WIFF file format."
            Thermo_RAW_format, //MS_Thermo_RAW_format, "MS:1000563", "Thermo RAW format", "Thermo Scientific RAW file format."
            Waters_raw_format, //MS_Waters_raw_format, "MS:1000526", "Waters raw format", "Waters data file format found in a Waters RAW directory, generated from an MS acquisition."
            Unknown
        };

        private Instrument _instrument;
        private int _artificialScanNum;

        private class BinaryDataArray
        {
            public int ArrayLength;
            public Precision Precision;
            public ArrayType ArrayType;
            public double[] Data;

            public BinaryDataArray()
            {
                Data = null;
                Precision = Precision.Precision32;
                ArrayType = ArrayType.m_z_array;
                ArrayLength = 0;
            }
        }

        private readonly List<Spectrum> _spectra = new List<Spectrum>();

        /// <summary>
        /// Initialize a MzMlReader object
        /// </summary>
        /// <param name="filePath">Path to mzML file</param>
        public MzMLReader(string filePath)
        {
            _instrument = Instrument.Unknown;
            _artificialScanNum = 1;
            _version = MzML_Version.mzML1_1_0;

            // Set a very large read buffer, it does decrease the read times for uncompressed files.
            Stream file = new FileStream(filePath, FileMode.Open, FileAccess.Read, FileShare.Read, 65536);

            if (filePath.EndsWith(".mzML.gz"))
            {
                file = new GZipStream(file, CompressionMode.Decompress);
            }

            var xSettings = new XmlReaderSettings { IgnoreWhitespace = true };
            _reader = XmlReader.Create(new StreamReader(file, System.Text.Encoding.UTF8, true, 65536), xSettings);
        }

        public IEnumerable<Spectrum> ReadAllSpectra()
        {
            ReadMzMl();
            return _spectra;
        }

        public void Close()
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Read and parse a .mzML file
        /// Files are commonly larger than 100 MB, so use a streaming reader instead of a DOM reader
        /// </summary>
        private void ReadMzMl()
        {
            // Handle disposal of allocated object correctly
            using (_reader)
            {
                // Guarantee a move to the root node
                _reader.MoveToContent();
                if (_reader.Name == "indexedmzML")
                {
                    // Read to the mzML root tag, and ignore the extra indexedmzML data
                    _reader.ReadToDescendant("mzML");
                    _reader = _reader.ReadSubtree();
                    _reader.MoveToContent();
                }
                string schemaName = _reader.GetAttribute("xsi:schemaLocation");
                // We automatically assume it uses the mzML_1.1.0 schema. Check for the old version.
                //if (!schemaName.Contains("mzML1.1.0.xsd"))
                if (schemaName.Contains("mzML1.0.0.xsd"))
                {
                    _version = MzML_Version.mzML1_0_0;
                }
                // Consume the mzML root tag
                // Throws exception if we are not at the "mzML" tag.
                // This is a critical error; we want to stop processing for this file if we encounter this error
                _reader.ReadStartElement("mzML");
                // Read the next node - should be the first child node
                while (_reader.ReadState == ReadState.Interactive)
                {
                    // Handle exiting out properly at EndElement tags
                    if (_reader.NodeType != XmlNodeType.Element)
                    {
                        _reader.Read();
                        continue;
                    }
                    // Handle each 1st level as a chunk
                    switch (_reader.Name)
                    {
                        case "cvList":
                            // Schema requirements: one instance of this element
                            _reader.Skip();
                            break;
                        case "fileDescription":
                            // Schema requirements: one instance of this element
                            ReadFileDescription(_reader.ReadSubtree());
                            _reader.ReadEndElement(); // "fileDescription" must have child nodes
                            break;
                        case "referenceableParamGroupList":
                            // Schema requirements: zero to one instances of this element
                            ReadReferenceableParamGroupList(_reader.ReadSubtree());
                            _reader.ReadEndElement(); // "referenceableParamGroupList" must have child nodes
                            break;
                        case "sampleList":
                            // Schema requirements: zero to one instances of this element
                            _reader.Skip();
                            break;
                        case "softwareList":
                            // Schema requirements: one instance of this element
                            _reader.Skip();
                            break;
                        case "scanSettingsList":
                            // Schema requirements: zero to one instances of this element
                            _reader.Skip();
                            break;
                        case "instrumentConfigurationList":
                            // Schema requirements: one instance of this element
                            _reader.Skip();
                            break;
                        case "dataProcessingList":
                            // Schema requirements: one instance of this element
                            _reader.Skip();
                            break;
                        case "acquisitionSettingsList": // mzML 1.0.0 compatibility
                            // Schema requirements: zero to one instances of this element
                            _reader.Skip();
                            break;
                        case "run":
                            // Schema requirements: one instance of this element
                            // Use reader.ReadSubtree() to provide an XmlReader that is only valid for the element and child nodes
                            ReadRunData(_reader.ReadSubtree());
                            // "run" might not have any child nodes
                            // We will either consume the EndElement, or the same element that was passed to ReadRunData (in case of no child nodes)
                            _reader.Read();
                            break;
                        default:
                            // We are not reading anything out of the tag, so bypass it
                            _reader.Skip();
                            break;
                    }
                }
            }
        }

        /// <summary>
        /// Handle the child nodes of the fileDescription element
        /// Called by ReadMzML (xml hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single "fileDescription" element</param>
        private void ReadFileDescription(XmlReader reader)
        {
            reader.MoveToContent();
            reader.ReadStartElement("fileDescription"); // Throws exception if we are not at the "fileDescription" tag.
            while (reader.ReadState == ReadState.Interactive)
            {
                // Handle exiting out properly at EndElement tags
                if (reader.NodeType != XmlNodeType.Element)
                {
                    reader.Read();
                    continue;
                }
                switch (reader.Name)
                {
                    case "fileContent":
                        // Schema requirements: one instance of this element
                        reader.Skip();
                        break;
                    case "sourceFileList":
                        // Schema requirements: zero to one instances of this element
                        ReadSourceFileList(reader.ReadSubtree());
                        reader.ReadEndElement(); // "sourceFileList" must have child nodes
                        break;
                    case "contact":
                        // Schema requirements: zero to many instances of this element
                        reader.Skip();
                        break;
                    default:
                        reader.Skip();
                        break;
                }
            }
            reader.Close();
        }

        /// <summary>
        /// Handle a single sourceFileList element and child nodes
        /// Called by ReadMzML (xml hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single sourceFileList element</param>
        private void ReadSourceFileList(XmlReader reader)
        {
            reader.MoveToContent();
            //int count = Convert.ToInt32(reader.GetAttribute("count"));
            reader.ReadStartElement("sourceFileList"); // Throws exception if we are not at the "sourceFileList" tag.
            while (reader.ReadState == ReadState.Interactive && _instrument == Instrument.Unknown)
            {
                // Handle exiting out properly at EndElement tags
                if (reader.NodeType != XmlNodeType.Element)
                {
                    reader.Read();
                    continue;
                }
                if (reader.Name == "sourceFile")
                {
                    // Schema requirements: one to many instances of this element
                    var innerReader = reader.ReadSubtree();
                    innerReader.MoveToContent();
                    innerReader.ReadStartElement("sourceFile"); // Throws exception if we are not at the "sourceFile" tag.
                    while (innerReader.ReadState == ReadState.Interactive)
                    {
                        // Handle exiting out properly at EndElement tags
                        if (innerReader.NodeType != XmlNodeType.Element)
                        {
                            innerReader.Read();
                            continue;
                        }
                        switch (innerReader.Name)
                        {
                            case "referenceableParamGroupRef":
                                // Schema requirements: zero to many instances of this element
                                innerReader.Skip();
                                break;
                            case "cvParam":
                                // Schema requirements: zero to many instances of this element
                                /* MUST supply a *child* term of MS:1000767 (native spectrum identifier format) only once
                                 *   e.g.: MS:1000768 (Thermo nativeID format)
                                 *   e.g.: MS:1000769 (Waters nativeID format)
                                 *   e.g.: MS:1000770 (WIFF nativeID format)
                                 *   e.g.: MS:1000771 (Bruker/Agilent YEP nativeID format)
                                 *   e.g.: MS:1000772 (Bruker BAF nativeID format)
                                 *   e.g.: MS:1000773 (Bruker FID nativeID format)
                                 *   e.g.: MS:1000774 (multiple peak list nativeID format)
                                 *   e.g.: MS:1000775 (single peak list nativeID format)
                                 *   e.g.: MS:1000776 (scan number only nativeID format)
                                 *   e.g.: MS:1000777 (spectrum identifier nativeID format)
                                 *   e.g.: MS:1000823 "Bruker U2 nativeID format"
                                 *   e.g.: MS:1000824 "no nativeID format"
                                 *   e.g.: MS:1000929 "Shimadzu Biotech nativeID format"
                                 *   e.g.: MS:1001480 "AB SCIEX TOF/TOF nativeID format"
                                 *   e.g.: MS:1001508 "Agilent MassHunter nativeID format"
                                 *   e.g.: MS:1001526 "spectrum from database integer nativeID format"
                                 *   e.g.: MS:1001528 "Mascot query number"
                                 *   e.g.: MS:1001531 "spectrum from ProteinScape database nativeID format"
                                 *   e.g.: MS:1001532 "spectrum from database string nativeID format"
                                 *   e.g.: MS:1001559 "AB SCIEX TOF/TOF T2D nativeID format"
                                 *   e.g.: MS:1001562 "Scaffold nativeID format"
                                 *   e.g.: MS:1002303 "Bruker Container nativeID format"
                                 *   et al.
                                 * MUST supply a *child* term of MS:1000561 (data file checksum type) one or more times
                                 *   e.g.: MS:1000568 (MD5)
                                 *   e.g.: MS:1000569 (SHA-1)
                                 * MUST supply a *child* term of MS:1000560 (mass spectrometer file format) only once
                                 *   e.g.: MS:1000526 (Waters raw file)             1.0.0: (MassLynx raw format)
                                 *   e.g.: MS:1000562 (ABI WIFF file)               1.0.0: (wiff file)
                                 *   e.g.: MS:1000563 (Thermo RAW file)             1.0.0: (Xcalibur RAW file)
                                 *   e.g.: MS:1000564 (PSI mzData file)             1.0.0: (mzData file)
                                 *   e.g.: MS:1000565 (Micromass PKL file)          1.0.0: (pkl file)
                                 *   e.g.: MS:1000566 (ISB mzXML file)              1.0.0: (mzXML file)
                                 *   e.g.: MS:1000567 (Bruker/Agilent YEP file)     1.0.0: (yep file)
                                 *   e.g.: MS:1000584 (mzML file)
                                 *   e.g.: MS:1000613 (DTA file)                    1.0.0: (dta file)
                                 *   e.g.: MS:1000614 (ProteinLynx Global Server mass spectrum XML file)
                                 *   e.g.: MS:1000740 "parameter file"
                                 *   e.g.: MS:1000742 "Bioworks SRF format"
                                 *   e.g.: MS:1000815 "Bruker BAF format"
                                 *   e.g.: MS:1000816 "Bruker U2 format"
                                 *   e.g.: MS:1000825 "Bruker FID format"
                                 *   e.g.: MS:1000930 "Shimadzu Biotech database entity"
                                 *   e.g.: MS:1001062 "Mascot MGF format"
                                 *   e.g.: MS:1001245 "PerSeptive PKS format"
                                 *   e.g.: MS:1001246 "Sciex API III format"
                                 *   e.g.: MS:1001247 "Bruker XML format"
                                 *   e.g.: MS:1001369 "text format"
                                 *   e.g.: MS:1001463 "Phenyx XML format"
                                 *   e.g.: MS:1001481 "AB SCIEX TOF/TOF database"
                                 *   e.g.: MS:1001509 "Agilent MassHunter format"
                                 *   e.g.: MS:1001527 "Proteinscape spectra"
                                 *   e.g.: MS:1001560 "AB SCIEX TOF/TOF T2D format"
                                 *   e.g.: MS:1001881 "mz5 format"
                                 *   e.g.: MS:1002302 "Bruker Container format"
                                 *   e.g.: MS:1002385 "SCiLS Lab format"
                                 *   e.g.: MS:1002441 "Andi-MS format"
                                 *   et al.
                                 */
                                switch (innerReader.GetAttribute("accession"))
                                {

                                    case "MS:1000768":
                                        // name="Thermo nativeID format"
                                        _instrument = Instrument.Thermo_RAW_format;
                                        break;
                                    case "MS:1000769":
                                        // name="Waters nativeID format"
                                        _instrument = Instrument.Waters_raw_format;
                                        break;
                                    case "MS:1000770":
                                        // name="WIFF nativeID format"
                                        _instrument = Instrument.ABI_WIFF_format;
                                        break;
                                    case "MS:1000771":
                                        // name="Bruker/Agilent YEP nativeID format"
                                        break;
                                    case "MS:1000772":
                                        // name="Bruker BAF nativeID format"
                                        break;
                                    case "MS:1000773":
                                        // name="Bruker FID nativeID format"
                                        break;
                                    case "MS:1000774":
                                        // name="multiple peak list nativeID format"
                                        break;
                                    case "MS:1000775":
                                        // name="single peak list nativeID format"
                                        break;
                                    case "MS:1000776":
                                        // name="scan number only nativeID format"
                                        break;
                                    case "MS:1000777":
                                        // name="spectrum identifier nativeID format"
                                        break;
                                    case "MS:1000823":
                                        // name="Bruker U2 nativeID format"
                                        break;
                                    case "MS:1000824":
                                        // name="no nativeID format"
                                        break;
                                    case "MS:1000929":
                                        // name="Shimadzu Biotech nativeID format"
                                        break;
                                    case "MS:1001480":
                                        // name="AB SCIEX TOF/TOF nativeID format"
                                        break;
                                    case "MS:1001508":
                                        // name="Agilent MassHunter nativeID format"
                                        break;
                                    case "MS:1001526":
                                        // name="spectrum from database integer nativeID format"
                                        break;
                                    case "MS:1001528":
                                        // name="Mascot query number"
                                        break;
                                    case "MS:1001531":
                                        // name="spectrum from ProteinScape database nativeID format"
                                        break;
                                    case "MS:1001532":
                                        // name="spectrum from database string nativeID format"
                                        break;
                                    case "MS:1001559":
                                        // name="AB SCIEX TOF/TOF T2D nativeID format"
                                        break;
                                    case "MS:1001562":
                                        // name="Scaffold nativeID format"
                                        break;
                                    case "MS:1002303":
                                        // name="Bruker Container nativeID format"
                                        break;
                                }
                                innerReader.Read(); // Consume the cvParam element (no child nodes)
                                break;
                            case "userParam":
                                // Schema requirements: zero to many instances of this element
                                innerReader.Skip();
                                break;
                            default:
                                innerReader.Skip();
                                break;
                        }
                    }
                    innerReader.Close();
                    reader.Read();
                }
                else
                {
                    reader.Read();
                }
            }
            reader.Close();
        }

        /// <summary>
        /// Handle the child nodes of the referenceableParamGroupList element
        /// Called by ReadMzML (xml hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single "referenceableParamGroupList" element</param>
        private void ReadReferenceableParamGroupList(XmlReader reader)
        {
            reader.MoveToContent();
            int count = Convert.ToInt32(reader.GetAttribute("count"));
            reader.ReadStartElement("referenceableParamGroupList"); // Throws exception if we are not at the "referenceableParamGroupList" tag.
            while (reader.ReadState == ReadState.Interactive)
            {
                // Handle exiting out properly at EndElement tags
                if (reader.NodeType != XmlNodeType.Element)
                {
                    reader.Read();
                    continue;
                }
                if (reader.Name == "referenceableParamGroup")
                {
                    // Schema requirements: one to many instances of this element
                    string id = reader.GetAttribute("id");
                    List<Param> paramList = new List<Param>();
                    var innerReader = reader.ReadSubtree();
                    innerReader.MoveToContent();
                    innerReader.ReadStartElement("referenceableParamGroup"); // Throws exception if we are not at the "sourceFile" tag.
                    while (innerReader.ReadState == ReadState.Interactive)
                    {
                        // Handle exiting out properly at EndElement tags
                        if (innerReader.NodeType != XmlNodeType.Element)
                        {
                            innerReader.Read();
                            continue;
                        }
                        switch (innerReader.Name)
                        {
                            case "cvParam":
                                // Schema requirements: zero to many instances of this element
                                paramList.Add(ReadCvParam(innerReader.ReadSubtree()));
                                innerReader.Read(); // Consume the cvParam element (no child nodes)
                                break;
                            case "userParam":
                                // Schema requirements: zero to many instances of this element
                                paramList.Add(ReadUserParam(innerReader.ReadSubtree()));
                                innerReader.Read(); // Consume the userParam element (no child nodes)
                                break;
                            default:
                                innerReader.Skip();
                                break;
                        }
                    }
                    innerReader.Close();
                    reader.Read();
                    _referenceableParamGroups.Add(id, paramList);
                }
                else
                {
                    reader.Read();
                }
            }
            reader.Close();
        }
        
        /// <summary>
        /// Handle the cvParam element
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single "cvParam" element</param>
        private CVParam ReadCvParam(XmlReader reader)
        {
            reader.MoveToContent();
            CVParam cvParam = new CVParam();
            cvParam.Accession = reader.GetAttribute("accession");
            cvParam.CVRef = reader.GetAttribute("cvRef");
            cvParam.Name = reader.GetAttribute("name");
            cvParam.Value = reader.GetAttribute("value");
            cvParam.UnitAccession = reader.GetAttribute("unitAccession");
            cvParam.UnitCVRef = reader.GetAttribute("unitCVRef");
            cvParam.UnitName = reader.GetAttribute("unitName");

            reader.Close();
            return cvParam;
        }

        /// <summary>
        /// Handle the userParam element
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single "userParam" element</param>
        private UserParam ReadUserParam(XmlReader reader)
        {
            reader.MoveToContent();
            UserParam userParam = new UserParam();
            userParam.Name = reader.GetAttribute("name");
            userParam.Type = reader.GetAttribute("type");
            userParam.Value = reader.GetAttribute("value");
            userParam.UnitAccession = reader.GetAttribute("unitAccession");
            userParam.UnitCVRef = reader.GetAttribute("unitCVRef");
            userParam.UnitName = reader.GetAttribute("unitName");

            reader.Close();
            return userParam;
        }

        /// <summary>
        /// Handle the child nodes of the run element
        /// Called by ReadMzML (xml hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single "run" element</param>
        private void ReadRunData(XmlReader reader)
        {
            reader.MoveToContent();
            reader.ReadStartElement("run"); // Throws exception if we are not at the "run" tag.
            while (reader.ReadState == ReadState.Interactive)
            {
                // Handle exiting out properly at EndElement tags
                if (reader.NodeType != XmlNodeType.Element)
                {
                    reader.Read();
                    continue;
                }
                switch (reader.Name)
                {
                    case "referenceableParamGroupRef":
                        // Schema requirements: zero to many instances of this element
                        reader.Skip();
                        break;
                    case "cvParam":
                        // Schema requirements: zero to many instances of this element
                        reader.Skip();
                        break;
                    case "userParam":
                        // Schema requirements: zero to many instances of this element
                        reader.Skip();
                        break;
                    case "sourceFileRefList": // mzML_1.0.0 compatibility
                        // Schema requirements: zero to many instances of this element
                        reader.Skip();
                        break;
                    case "spectrumList":
                        // Schema requirements: zero to one instances of this element
                        // Use reader.ReadSubtree() to provide an XmlReader that is only valid for the element and child nodes
                        ReadSpectrumList(reader.ReadSubtree());
                        // "spectrumList" might not have any child nodes
                        // We will either consume the EndElement, or the same element that was passed to ReadSpectrumList (in case of no child nodes)
                        reader.Read();
                        break;
                    case "chromatogramList":
                        // Schema requirements: zero to one instances of this element
                        reader.Skip();
                        break;
                    default:
                        reader.Skip();
                        break;
                }
            }
            reader.Close();
        }

        /// <summary>
        /// Handle the child nodes of a spectrumList element
        /// Called by ReadRunData (xml hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single spectrumList element</param>
        private void ReadSpectrumList(XmlReader reader)
        {
            reader.MoveToContent();
            reader.ReadStartElement("spectrumList"); // Throws exception if we are not at the "SpectrumIdentificationList" tag.
            while (reader.ReadState == ReadState.Interactive)
            {
                // Handle exiting out properly at EndElement tags
                if (reader.NodeType != XmlNodeType.Element)
                {
                    reader.Read();
                    continue;
                }
                if (reader.Name == "spectrum")
                {
                    // Schema requirements: zero to many instances of this element
                    // Use reader.ReadSubtree() to provide an XmlReader that is only valid for the element and child nodes
                    ReadSpectrum(reader.ReadSubtree());
                    // "spectrum" might not have any child nodes
                    // We will either consume the EndElement, or the same element that was passed to ReadSpectrum (in case of no child nodes)
                    reader.Read();
                }
                else
                {
                    reader.Skip();
                }
            }
            reader.Close();
        }

        /// <summary>
        /// Handle a single spectrum element and child nodes
        /// Called by ReadSpectrumList (xml hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single spectrum element</param>
        private void ReadSpectrum(XmlReader reader)
        {
            reader.MoveToContent();
            string index = reader.GetAttribute("index");
            //Console.WriteLine("Reading spectrum indexed by " + index);
            // This is correct for Thermo files converted by msConvert, but need to implement for others as well
            string spectrumId = reader.GetAttribute("id"); // Native ID in mzML_1.1.0; unique identifier in mzML_1.0.0, often same as nativeID
            string nativeId = spectrumId;
            if (_version == MzML_Version.mzML1_0_0)
            {
                nativeId = reader.GetAttribute("nativeID"); // Native ID in mzML_1.0.0
            }
            //int scanNum = Convert.ToInt32(spectrumId.Substring(spectrumId.LastIndexOf("scan=") + 5));
            // TODO: Get rid of this hack, use something with nativeID. May involve special checks for mzML version
            int scanNum = _artificialScanNum++;
            int defaultArraySize = Convert.ToInt32(reader.GetAttribute("defaultArrayLength"));
            reader.ReadStartElement("spectrum"); // Throws exception if we are not at the "spectrum" tag.
            bool is_ms_ms = false;
            int msLevel = 0;
            bool centroided = false;
            List<Precursor> precursors = new List<Precursor>();
            List<ScanData> scans = new List<ScanData>();
            List<BinaryDataArray> bdas = new List<BinaryDataArray>();
            while (reader.ReadState == ReadState.Interactive)
            {
                // Handle exiting out properly at EndElement tags
                if (reader.NodeType != XmlNodeType.Element)
                {
                    reader.Read();
                    continue;
                }
                //////////////////////////////////////////////////////////////////////////////////////
                /// 
                /// MS1 Spectra: only need Spectrum data: scanNum, MSLevel, ElutionTime, mzArray, IntensityArray
                /// 
                /// MS2 Spectra: use ProductSpectrum; adds ActivationMethod and IsolationWindow
                /// 
                //////////////////////////////////////////////////////////////////////////////////////
                switch (reader.Name)
                {
                    case "referenceableParamGroupRef":
                        // Schema requirements: zero to many instances of this element
                        reader.Skip();
                        break;
                    case "cvParam":
                        // Schema requirements: zero to many instances of this element
                        /* MAY supply a *child* term of MS:1000465 (scan polarity) only once
                         *   e.g.: MS:1000129 (negative scan)
                         *   e.g.: MS:1000130 (positive scan)
                         * MUST supply a *child* term of MS:1000559 (spectrum type) only once
                         *   e.g.: MS:1000322 (charge inversion mass spectrum)
                         *   e.g.: MS:1000325 (constant neutral gain spectrum)
                         *   e.g.: MS:1000326 (constant neutral loss spectrum)
                         *   e.g.: MS:1000328 (e/2 mass spectrum)
                         *   e.g.: MS:1000341 (precursor ion spectrum)
                         *   e.g.: MS:1000579 (MS1 spectrum)
                         *   e.g.: MS:1000580 (MSn spectrum)
                         *   e.g.: MS:1000581 (CRM spectrum)
                         *   e.g.: MS:1000582 (SIM spectrum)
                         *   e.g.: MS:1000583 (SRM spectrum)
                         *   e.g.: MS:1000620 (PDA spectrum)
                         *   e.g.: MS:1000627 (selected ion current chromatogram)
                         *   e.g.: MS:1000789 (enhanced multiply charged spectrum)
                         *   e.g.: MS:1000790 (time-delayed fragmentation spectrum)
                         *   et al.
                         * MUST supply term MS:1000525 (spectrum representation) or any of its children only once
                         *   e.g.: MS:1000127 (centroid spectrum)
                         *   e.g.: MS:1000128 (profile spectrum)
                         * MAY supply a *child* term of MS:1000499 (spectrum attribute) one or more times
                         *   e.g.: MS:1000285 (total ion current)
                         *   e.g.: MS:1000497 (zoom scan)
                         *   e.g.: MS:1000504 (base peak m/z)
                         *   e.g.: MS:1000505 (base peak intensity)
                         *   e.g.: MS:1000511 (ms level)
                         *   e.g.: MS:1000527 (highest observed m/z)
                         *   e.g.: MS:1000528 (lowest observed m/z)
                         *   e.g.: MS:1000618 (highest observed wavelength)
                         *   e.g.: MS:1000619 (lowest observed wavelength)
                         *   e.g.: MS:1000796 (spectrum title)
                         *   et al.
                         */
                        switch (reader.GetAttribute("accession"))
                        {
                            case "MS:1000127":
                                // name="centroid spectrum"
                                centroided = true;
                                break;
                            case "MS:1000128":
                                // name="profile spectrum"
                                is_ms_ms = false;
                                break;
                            case "MS:1000511":
                                // name="ms level"
                                msLevel = Convert.ToInt32(reader.GetAttribute("value"));
                                break;
                            case "MS:1000579":
                                // name="MS1 spectrum"
                                is_ms_ms = false;
                                break;
                            case "MS:1000580":
                                // name="MSn spectrum"
                                is_ms_ms = true;
                                break;
                        }
                        reader.Read(); // Consume the cvParam element (no child nodes)
                        break;
                    case "userParam":
                        // Schema requirements: zero to many instances of this element
                        reader.Skip();
                        break;
                    case "spectrumDescription": // mzML_1.0.0 compatibility
                        // Schema requirements: one instance of this element
                        ReadSpectrumDescription(reader.ReadSubtree(), ref scans, ref precursors, out centroided);
                        reader.ReadEndElement(); // "spectrumDescription" must have child nodes
                        break;
                    case "scanList":
                        // Schema requirements: zero to one instances of this element
                        scans.AddRange(ReadScanList(reader.ReadSubtree()));
                        reader.ReadEndElement(); // "scanList" must have child nodes
                        break;
                    case "precursorList":
                        // Schema requirements: zero to one instances of this element
                        precursors.AddRange(ReadPrecursorList(reader.ReadSubtree()));
                        reader.ReadEndElement(); // "precursorList" must have child nodes
                        break;
                    case "productList":
                        // Schema requirements: zero to one instances of this element
                        reader.Skip();
                        break;
                    case "binaryDataArrayList":
                        // Schema requirements: zero to one instances of this element
                        bdas.AddRange(ReadBinaryDataArrayList(reader.ReadSubtree(), defaultArraySize));
                        reader.ReadEndElement(); // "binaryDataArrayList" must have child nodes
                        break;
                    default:
                        reader.Skip();
                        break;
                }
            }
            reader.Close();
            // Process the spectrum data
            BinaryDataArray mzs = new BinaryDataArray();
            BinaryDataArray intensities = new BinaryDataArray();
            ScanData scan = new ScanData();
            Spectrum spectrum;
            foreach (var bda in bdas)
            {
                if (bda.ArrayType == ArrayType.m_z_array)
                {
                    mzs = bda;
                }
                else if (bda.ArrayType == ArrayType.intensity_array)
                {
                    intensities = bda;
                }
            }
            if (!centroided)
            {
                // Centroid spectrum
                // ProteoWizard
                var centroider = new Centroider(mzs.Data, intensities.Data);
                double[] centroidedMzs, centroidedIntensities;
                centroider.GetCentroidedData(out centroidedMzs, out centroidedIntensities);
                mzs.Data = centroidedMzs;
                intensities.Data = centroidedIntensities;
            }
            if (scans.Count == 1)
            {
                scan = scans[0];
            }
            else if (scans.Count > 1)
            {
                // TODO: Should do something else to appropriately handle combinations...
                scan = scans[0];
            }
            if (is_ms_ms)
            {
                Precursor precursor = new Precursor();
                if (precursors.Count == 1)
                {
                    precursor = precursors[0];
                }
                else if (precursors.Count > 1)
                {
                    // TODO: Should do something else to appropriately handle multiple precursors...
                    precursor = precursors[0];
                }
                SelectedIon ion = new SelectedIon();
                if (precursor.Ions.Count == 1)
                {
                    ion = precursor.Ions[0];
                }
                else if (precursor.Ions.Count > 1)
                {
                    // TODO: Should do something else to appropriately handle multiple selected ions...
                    ion = precursor.Ions[0];
                }

                var pspectrum = new ProductSpectrum(mzs.Data, intensities.Data, scanNum);
                pspectrum.ActivationMethod = precursor.Activation;
                // Select mz value to use based on presence of a Thermo-specific user param.
                // The user param has a slightly higher precision, if that matters.
                double mz = scan.MonoisotopicMz == 0.0 ? ion.SelectedIonMz : scan.MonoisotopicMz;
                pspectrum.IsolationWindow = new IsolationWindow(precursor.IsolationWindowTargetMz, precursor.IsolationWindowLowerOffset, precursor.IsolationWindowUpperOffset, mz, ion.Charge);
                spectrum = pspectrum;
            }
            else
            {
                spectrum = new Spectrum(mzs.Data, intensities.Data, scanNum);
            }
            spectrum.MsLevel = msLevel;
            spectrum.ElutionTime = scan.StartTime;
            spectrum.NativeId = nativeId;
            _spectra.Add(spectrum);
        }

        /// <summary>
        /// mzML_1.0.0 compatibility
        /// Handle a single spectrumDescription element and child nodes
        /// Called by ReadSpectrumList (xml hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single spectrum element</param>
        /// <param name="scans"></param>
        /// <param name="precursors"></param>
        /// <param name="centroided"></param>
        private void ReadSpectrumDescription(XmlReader reader, ref List<ScanData> scans, ref List<Precursor> precursors, out bool centroided)
        {
            reader.MoveToContent();
            // This is correct for Thermo files converted by msConvert, but need to implement for others as well
            reader.ReadStartElement("spectrumDescription"); // Throws exception if we are not at the "spectrumDescription" tag.
            centroided = false;
            while (reader.ReadState == ReadState.Interactive)
            {
                // Handle exiting out properly at EndElement tags
                if (reader.NodeType != XmlNodeType.Element)
                {
                    reader.Read();
                    continue;
                }
                switch (reader.Name)
                {
                    case "referenceableParamGroupRef":
                        // Schema requirements: zero to many instances of this element
                        reader.Skip();
                        break;
                    case "cvParam":
                        // Schema requirements: zero to many instances of this element
                        /* MUST supply a *child* term of MS:1000525 (spectrum representation) only once
                         *   e.g.: MS:1000127 (centroid mass spectrum)
                         *   e.g.: MS:1000128 (profile mass spectrum)
                         * MUST supply a *child* term of MS:1000499 (spectrum attribute) one or more times
                         *   e.g.: MS:1000285 (total ion current)
                         *   e.g.: MS:1000504 (base peak m/z)
                         *   e.g.: MS:1000505 (base peak intensity)
                         *   e.g.: MS:1000527 (highest m/z value)
                         *   e.g.: MS:1000528 (lowest m/z value)
                         *   e.g.: MS:1000618 (highest wavelength value)
                         *   e.g.: MS:1000619 (lowest wavelength value)
                         */
                        switch (reader.GetAttribute("accession"))
                        {
                            case "MS:1000127":
                                // name="centroid spectrum"
                                centroided = true;
                                break;
                        }
                        reader.Read(); // Consume the cvParam element (no child nodes)
                        break;
                    case "userParam":
                        // Schema requirements: zero to many instances of this element
                        reader.Skip();
                        break;
                    case "acquisitionList":
                        // Schema requirements: zero to one instances of this element
                        // Very comparable to mzML_1.1.0's scanList. Use it.
                        scans.AddRange(ReadScanList(reader.ReadSubtree()));
                        reader.ReadEndElement(); // "acquisitionList" mustt have child nodes
                        break;
                    case "precursorList":
                        // Schema requirements: zero to one instances of this element
                        precursors.AddRange(ReadPrecursorList(reader.ReadSubtree()));
                        reader.ReadEndElement(); // "precursorList" must have child nodes
                        break;
                    case "scan":
                        // Schema requirements: zero to one instances of this element
                        scans.Add(ReadScan(reader.ReadSubtree()));
                        reader.Read(); // "scan" might not have child nodes
                        break;
                    default:
                        reader.Skip();
                        break;
                }
            }
            reader.Close();
        }

        /// <summary>
        /// Handle a single scanList element and child nodes
        /// Called by ReadSpectrum (xml hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single scanList element</param>
        /// <returns></returns>
        private List<ScanData> ReadScanList(XmlReader reader)
        {
            reader.MoveToContent();
            int count = Convert.ToInt32(reader.GetAttribute("count"));
            List<ScanData> scans = new List<ScanData>();
            if (_version == MzML_Version.mzML1_0_0)
            {
                reader.ReadStartElement("acquisitionList"); // Throws exception if we are not at the "scanList" tag.
            }
            else
            {
                reader.ReadStartElement("scanList"); // Throws exception if we are not at the "scanList" tag.
            }
            while (reader.ReadState == ReadState.Interactive)
            {
                // Handle exiting out properly at EndElement tags
                if (reader.NodeType != XmlNodeType.Element)
                {
                    reader.Read();
                    continue;
                }

                switch (reader.Name)
                {
                    case "referenceableParamGroupRef":
                        // Schema requirements: zero to many instances of this element
                        reader.Skip();
                        break;
                    case "cvParam":
                        // Schema requirements: zero to many instances of this element
                        /* MUST supply a *child* term of MS:1000570 (spectra combination) only once
                         *   e.g.: MS:1000571 (sum of spectra)
                         *   e.g.: MS:1000573 (median of spectra)
                         *   e.g.: MS:1000575 (mean of spectra)
                         *   e.g.: MS:1000795 (no combination)
                         */
                        switch (reader.GetAttribute("accession"))
                        {
                            case "MS:1000571":
                                // name="sum of spectra"
                                break;
                            case "MS:1000573":
                                // name="median of spectra"
                                break;
                            case "MS:1000575":
                                // name="mean of spectra"
                                break;
                            case "MS:1000795":
                                // name="no combination"
                                break;
                        }
                        reader.Read(); // Consume the cvParam element (no child nodes)
                        break;
                    case "userParam":
                        // Schema requirements: zero to many instances of this element
                        reader.Skip();
                        break;
                    case "scan":
                    case "acquisition": // mzML_1.0.0 compatibility
                        // Schema requirements: one to many instances of this element
                        scans.Add(ReadScan(reader.ReadSubtree()));
                        reader.Read(); // "scan" might not have child nodes
                        break;
                    default:
                        reader.Skip();
                        break;
                }
            }
            reader.Close();
            return scans;
        }

        /// <summary>
        /// Handle a single scan element and child nodes
        /// Called by ReadSpectrum (xml hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single scan element</param>
        /// <returns></returns>
        private ScanData ReadScan(XmlReader reader)
        {
            reader.MoveToContent();
            if (_version == MzML_Version.mzML1_0_0)
            {
                string name = reader.Name;
                if (!name.Equals("scan") && !name.Equals("acquisition"))
                {
                    throw new XmlException("Invalid schema");
                }
                reader.ReadStartElement(name);
            }
            else
            {
                reader.ReadStartElement("scan"); // Throws exception if we are not at the "scan" tag.
            }
            ScanData scan = new ScanData();
            while (reader.ReadState == ReadState.Interactive)
            {
                // Handle exiting out properly at EndElement tags
                if (reader.NodeType != XmlNodeType.Element)
                {
                    reader.Read();
                    continue;
                }

                switch (reader.Name)
                {
                    case "referenceableParamGroupRef":
                        // Schema requirements: zero to many instances of this element
                        reader.Skip();
                        break;
                    case "cvParam":
                        // Schema requirements: zero to many instances of this element
                        /* MAY supply a *child* term of MS:1000503 (scan attribute) one or more times
                         *   e.g.: MS:1000011 (mass resolution)
                         *   e.g.: MS:1000015 (scan rate)
                         *   e.g.: MS:1000016 (scan start time)
                         *   e.g.: MS:1000502 (dwell time)
                         *   e.g.: MS:1000512 (filter string)
                         *   e.g.: MS:1000616 (preset scan configuration)
                         *   e.g.: MS:1000800 (mass resolving power)
                         *   e.g.: MS:1000803 (analyzer scan offset)
                         *   e.g.: MS:1000826 (elution time)
                         *   e.g.: MS:1000880 (interchannel delay)
                         * MAY supply a *child* term of MS:1000018 (scan direction) only once
                         *   e.g.: MS:1000092 (decreasing m/z scan)
                         *   e.g.: MS:1000093 (increasing m/z scan)
                         * MAY supply a *child* term of MS:1000019 (scan law) only once
                         *   e.g.: MS:1000094 (exponential)
                         *   e.g.: MS:1000095 (linear)
                         *   e.g.: MS:1000096 (quadratic)
                         */
                        switch (reader.GetAttribute("accession"))
                        {
                            case "MS:1000016":
                                // name="scan start time"
                                double time = Convert.ToDouble(reader.GetAttribute("value"));
                                bool isSeconds = reader.GetAttribute("unitName") == "second";
                                // Should only see "second" and "minute"
                                scan.StartTime = isSeconds ? time / 60.0 : time;
                                //scan.StartTime = Convert.ToDouble(reader.GetAttribute("value"));
                                break;
                            case "MS:1000512":
                                // name="filter string"
                                break;
                            case "MS:1000616":
                                // name="preset scan configuration"
                                break;
                            case "MS:1000927":
                                // name="ion injection time"
                                break;
                            case "MS:1000826":
                                // name="elution time"
                                //startTime = Convert.ToDouble(reader.GetAttribute("value"));
                                break;
                        }
                        reader.Read(); // Consume the cvParam element (no child nodes)
                        break;
                    case "userParam":
                        // Schema requirements: zero to many instances of this element
                        if (reader.GetAttribute("name") == "[Thermo Trailer Extra]Monoisotopic M/Z:")
                        {
                            scan.MonoisotopicMz = Convert.ToDouble(reader.GetAttribute("value"));
                        }
                        reader.Read(); // Consume the userParam element (no child nodes)
                        break;
                    case "scanWindowList":
                        // Schema requirements: zero to one instances of this element
                        //ReadScanList(reader.ReadSubtree());
                        //reader.ReadEndElement(); // "scanWindowList" must have child nodes
                        reader.Skip();
                        break;
                    default:
                        reader.Skip();
                        break;
                }
            }
            reader.Close();
            return scan;
        }

        /// <summary>
        /// Handle a single precursorList element and child nodes
        /// Called by ReadSpectrum (xml hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single precursorList element</param>
        /// <returns></returns>
        private List<Precursor> ReadPrecursorList(XmlReader reader)
        {
            reader.MoveToContent();
            int count = Convert.ToInt32(reader.GetAttribute("count"));
            List<Precursor> precursors = new List<Precursor>();
            reader.ReadStartElement("precursorList"); // Throws exception if we are not at the "precursorList" tag.
            while (reader.ReadState == ReadState.Interactive)
            {
                // Handle exiting out properly at EndElement tags
                if (reader.NodeType != XmlNodeType.Element)
                {
                    reader.Read();
                    continue;
                }

                switch (reader.Name)
                {
                    case "precursor":
                        // Schema requirements: one to many instances of this element
                        precursors.Add(ReadPrecursor(reader.ReadSubtree()));
                        reader.ReadEndElement(); // "SpectrumIdentificationItem" must have child nodes
                        break;
                    default:
                        reader.Skip();
                        break;
                }
            }
            reader.Close();
            return precursors;
        }

        /// <summary>
        /// Handle a single precursor element and child nodes
        /// Called by ReadPrecursorList (xml hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single precursor element</param>
        /// <returns></returns>
        private Precursor ReadPrecursor(XmlReader reader)
        {
            reader.MoveToContent();
            reader.ReadStartElement("precursor"); // Throws exception if we are not at the "precursor" tag.
            XmlReader innerReader;
            Precursor precursor = new Precursor();
            while (reader.ReadState == ReadState.Interactive)
            {
                // Handle exiting out properly at EndElement tags
                if (reader.NodeType != XmlNodeType.Element)
                {
                    reader.Read();
                    continue;
                }

                switch (reader.Name)
                {
                    case "isolationWindow":
                        // Schema requirements: zero to one instances of this element
                        innerReader = reader.ReadSubtree();
                        innerReader.MoveToContent();
                        innerReader.ReadStartElement("isolationWindow"); // Throws exception if we are not at the "selectedIon" tag.
                        while (innerReader.ReadState == ReadState.Interactive)
                        {
                            // Handle exiting out properly at EndElement tags
                            if (innerReader.NodeType != XmlNodeType.Element)
                            {
                                innerReader.Read();
                                continue;
                            }
                            switch (innerReader.Name)
                            {
                                case "referenceableParamGroupRef":
                                    // Schema requirements: zero to many instances of this element
                                    innerReader.Skip();
                                    break;
                                case "cvParam":
                                    // Schema requirements: zero to many instances of this element
                                    /* MUST supply a *child* term of MS:1000792 (isolation window attribute) one or more times
                                     *   e.g.: MS:1000827 (isolation window target m/z)
                                     *   e.g.: MS:1000828 (isolation window lower offset)
                                     *   e.g.: MS:1000829 (isolation window upper offset)
                                     */
                                    switch (innerReader.GetAttribute("accession"))
                                    {
                                        case "MS:1000827":
                                            // name="isolation window target m/z"
                                            precursor.IsolationWindowTargetMz = Convert.ToDouble(innerReader.GetAttribute("value"));
                                            break;
                                        case "MS:1000828":
                                            // name="isolation window lower offset"
                                            precursor.IsolationWindowLowerOffset = Convert.ToDouble(innerReader.GetAttribute("value"));
                                            break;
                                        case "MS:1000829":
                                            // name="isolation window upper offset"
                                            precursor.IsolationWindowUpperOffset = Convert.ToDouble(innerReader.GetAttribute("value"));
                                            break;
                                    }
                                    innerReader.Read(); // Consume the cvParam element (no child nodes)
                                    break;
                                case "userParam":
                                    // Schema requirements: zero to many instances of this element
                                    innerReader.Skip();
                                    break;
                                default:
                                    innerReader.Skip();
                                    break;
                            }
                        }
                        innerReader.Close();
                        reader.Read(); // "selectedIon" might not have child nodes
                        // We will either consume the EndElement, or the same element that was passed to ReadSpectrum (in case of no child nodes)
                        break;
                    case "selectedIonList":
                        // Schema requirements: zero to one instances of this element
                        // mzML_1.0.0: one instance of this element
                        precursor.Ions = ReadSelectedIonList(reader.ReadSubtree());
                        reader.ReadEndElement(); // "selectedIonList" must have child nodes
                        break;
                    case "activation":
                        // Schema requirements: one instance of this element
                        innerReader = reader.ReadSubtree();
                        innerReader.MoveToContent();
                        innerReader.ReadStartElement("activation"); // Throws exception if we are not at the "activation" tag.
                        while (innerReader.ReadState == ReadState.Interactive)
                        {
                            // Handle exiting out properly at EndElement tags
                            if (innerReader.NodeType != XmlNodeType.Element)
                            {
                                innerReader.Read();
                                continue;
                            }
                            switch (innerReader.Name)
                            {
                                case "referenceableParamGroupRef":
                                    // Schema requirements: zero to many instances of this element
                                    innerReader.Skip();
                                    break;
                                case "cvParam":
                                    // Schema requirements: zero to many instances of this element
                                    /* MAY supply a *child* term of MS:1000510 (precursor activation attribute) one or more times
                                     *   e.g.: MS:1000045 (collision energy)
                                     *   e.g.: MS:1000138 (percent collision energy)
                                     *   e.g.: MS:1000245 (charge stripping)
                                     *   e.g.: MS:1000412 (buffer gas)
                                     *   e.g.: MS:1000419 (collision gas)
                                     *   e.g.: MS:1000509 (activation energy)
                                     *   e.g.: MS:1000869 (collision gas pressure)
                                     * MUST supply term MS:1000044 (dissociation method) or any of its children one or more times
                                     *   e.g.: MS:1000133 (collision-induced dissociation)
                                     *   e.g.: MS:1000134 (plasma desorption)
                                     *   e.g.: MS:1000135 (post-source decay)
                                     *   e.g.: MS:1000136 (surface-induced dissociation)
                                     *   e.g.: MS:1000242 (blackbody infrared radiative dissociation)
                                     *   e.g.: MS:1000250 (electron capture dissociation)
                                     *   e.g.: MS:1000262 (infrared multiphoton dissociation)
                                     *   e.g.: MS:1000282 (sustained off-resonance irradiation)
                                     *   e.g.: MS:1000422 (high-energy collision-induced dissociation)
                                     *   e.g.: MS:1000433 (low-energy collision-induced dissociation)
                                     *   et al.
                                     *   
                                     *   e.g.: MS:1000133 "collision-induced dissociation"
                                     *   e.g.: MS:1000134 "plasma desorption"
                                     *   e.g.: MS:1000135 "post-source decay"
                                     *   e.g.: MS:1000136 "surface-induced dissociation"
                                     *   e.g.: MS:1000242 "blackbody infrared radiative dissociation"
                                     *   e.g.: MS:1000250 "electron capture dissociation"
                                     *   e.g.: MS:1000262 "infrared multiphoton dissociation"
                                     *   e.g.: MS:1000282 "sustained off-resonance irradiation"
                                     *   e.g.: MS:1000422 "beam-type collision-induced dissociation"
                                     *   e.g.: MS:1000433 "low-energy collision-induced dissociation"
                                     *   e.g.: MS:1000435 "photodissociation"
                                     *   e.g.: MS:1000598 "electron transfer dissociation"
                                     *   e.g.: MS:1000599 "pulsed q dissociation"
                                     *   e.g.: MS:1001880 "in-source collision-induced dissociation"
                                     *   e.g.: MS:1002000 "LIFT"
                                     *   e.g.: MS:1002472 "trap-type collision-induced dissociation"
                                     */
                                    switch (innerReader.GetAttribute("accession"))
                                    {
                                        case "MS:1000133":
                                            // name="collision-induced dissociation"
                                            precursor.Activation = ActivationMethod.CID;
                                            break;
                                        case "MS:1000598":
                                            // name="electron transfer dissociation"
                                            precursor.Activation = ActivationMethod.ETD;
                                            break;
                                        case "MS:1000422":
                                            // name="beam-type collision-induced dissociation", "high-energy collision-induced dissociation"
                                            precursor.Activation = ActivationMethod.HCD;
                                            break;
                                        case "MS:1000250":
                                            // name="electron capture dissociation"
                                            precursor.Activation = ActivationMethod.ECD;
                                            break;
                                        case "MS:1000599":
                                            // name="pulsed q dissociation"
                                            precursor.Activation = ActivationMethod.PQD;
                                            break;
                                        case "MS:1000045":
                                            // name="collision energy"
                                            //energy = Convert.ToDouble(innerReader.GetAttribute("value"));
                                            break;
                                    }
                                    innerReader.Read(); // Consume the cvParam element (no child nodes)
                                    break;
                                case "userParam":
                                    // Schema requirements: zero to many instances of this element
                                    innerReader.Skip();
                                    break;
                                default:
                                    innerReader.Skip();
                                    break;
                            }
                        }
                        innerReader.Close();
                        reader.Read(); // "selectedIon" might not have child nodes
                        // We will either consume the EndElement, or the same element that was passed to ReadSpectrum (in case of no child nodes)
                        break;
                    default:
                        reader.Skip();
                        break;
                }
            }
            reader.Close();
            return precursor;
        }

        /// <summary>
        /// Handle a single selectedIonList element and child nodes
        /// Called by ReadPrecursor (xml hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single selectedIonList element</param>
        /// <returns></returns>
        private List<SelectedIon> ReadSelectedIonList(XmlReader reader)
        {
            reader.MoveToContent();
            //int count = Convert.ToInt32(reader.GetAttribute("count"));
            var ions = new List<SelectedIon>();
            reader.ReadStartElement("selectedIonList"); // Throws exception if we are not at the "selectedIonList" tag.
            while (reader.ReadState == ReadState.Interactive)
            {
                // Handle exiting out properly at EndElement tags
                if (reader.NodeType != XmlNodeType.Element)
                {
                    reader.Read();
                    continue;
                }

                switch (reader.Name)
                {
                    case "selectedIon":
                        // Schema requirements: one to many instances of this element
                        var ion = new SelectedIon();
                        var innerReader = reader.ReadSubtree();
                        innerReader.MoveToContent();
                        innerReader.ReadStartElement("selectedIon"); // Throws exception if we are not at the "selectedIon" tag.
                        while (innerReader.ReadState == ReadState.Interactive)
                        {
                            // Handle exiting out properly at EndElement tags
                            if (innerReader.NodeType != XmlNodeType.Element)
                            {
                                innerReader.Read();
                                continue;
                            }
                            switch (innerReader.Name)
                            {
                                case "referenceableParamGroupRef":
                                    // Schema requirements: zero to many instances of this element
                                    innerReader.Skip();
                                    break;
                                case "cvParam":
                                    // Schema requirements: zero to many instances of this element
                                    /* MUST supply a *child* term of MS:1000455 (ion selection attribute) one or more times
                                     *   e.g.: MS:1000041 (charge state)
                                     *   e.g.: MS:1000042 (intensity)
                                     *   e.g.: MS:1000633 (possible charge state)
                                     *   e.g.: MS:1000744 (selected ion m/z)
                                     */
                                    switch (innerReader.GetAttribute("accession"))
                                    {
                                        case "MS:1000041":
                                            // name="charge state"
                                            ion.Charge = (int)Convert.ToDouble(innerReader.GetAttribute("value"));
                                            break;
                                        case "MS:1000744":
                                            // name="selected ion m/z"
                                            ion.SelectedIonMz = Convert.ToDouble(innerReader.GetAttribute("value"));
                                            break;
                                    }
                                    innerReader.Read(); // Consume the cvParam element (no child nodes)
                                    break;
                                case "userParam":
                                    // Schema requirements: zero to many instances of this element
                                    innerReader.Skip();
                                    break;
                                default:
                                    innerReader.Skip();
                                    break;
                            }
                        }
                        innerReader.Close();
                        ions.Add(ion);
                        reader.Read(); // "selectedIon" might not have child nodes
                        // We will either consume the EndElement, or the same element that was passed to ReadSpectrum (in case of no child nodes)
                        break;
                    default:
                        reader.Skip();
                        break;
                }
            }
            reader.Close();
            return ions;
        }

        /// <summary>
        /// Handle a single binaryDataArrayList element and child nodes
        /// Called by ReadSpectrum (xml hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single binaryDataArrayList element</param>
        /// <param name="defaultArrayLength">Default array length, coming from spectrum attribute</param>
        /// <returns></returns>
        private List<BinaryDataArray> ReadBinaryDataArrayList(XmlReader reader, int defaultArrayLength)
        {
            reader.MoveToContent();
            int bdArrays = Convert.ToInt32(reader.GetAttribute("count"));
            List<BinaryDataArray> bdaList = new List<BinaryDataArray>();
            reader.ReadStartElement("binaryDataArrayList"); // Throws exception if we are not at the "binaryDataArrayList" tag.
            while (reader.ReadState == ReadState.Interactive)
            {
                // Handle exiting out properly at EndElement tags
                if (reader.NodeType != XmlNodeType.Element)
                {
                    reader.Read();
                    continue;
                }

                switch (reader.Name)
                {
                    case "binaryDataArray":
                        // Schema requirements: two to many instances of this element
                        bdaList.Add(ReadBinaryDataArray(reader.ReadSubtree(), defaultArrayLength));
                        reader.ReadEndElement(); // "SpectrumIdentificationItem" must have child nodes
                        break;
                    default:
                        reader.Skip();
                        break;
                }
            }
            reader.Close();
            return bdaList;
        }

        /// <summary>
        /// Handle a single binaryDataArray element and child nodes
        /// Called by ReadBinaryDataArrayList (xml hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single binaryDataArray element</param>
        /// <param name="defaultLength">Default array length, coming from spectrum attribute</param>
        /// <returns></returns>
        private BinaryDataArray ReadBinaryDataArray(XmlReader reader, int defaultLength)
        {
            reader.MoveToContent();
            BinaryDataArray bda = new BinaryDataArray();
            bda.ArrayLength = defaultLength;
            int encLength = Convert.ToInt32(reader.GetAttribute("encodedLength"));
            int arrLength = Convert.ToInt32(reader.GetAttribute("arrayLength")); // Override the default; if non-existent, should get 0
            if (arrLength > 0)
            {
                bda.ArrayLength = arrLength;
            }
            bool compressed = false;
            reader.ReadStartElement("binaryDataArray"); // Throws exception if we are not at the "spectrum" tag.
            List<Param> paramList = new List<Param>();
            while (reader.ReadState == ReadState.Interactive)
            {
                // Handle exiting out properly at EndElement tags
                if (reader.NodeType != XmlNodeType.Element)
                {
                    reader.Read();
                    continue;
                }

                switch (reader.Name)
                {
                    case "referenceableParamGroupRef":
                        // Schema requirements: zero to many instances of this element
                        string rpgRef = reader.GetAttribute("ref");
                        paramList.AddRange(_referenceableParamGroups[rpgRef]);
                        reader.Read();
                        break;
                    case "cvParam":
                        // Schema requirements: zero to many instances of this element
                        paramList.Add(ReadCvParam(reader.ReadSubtree()));
                        reader.Read(); // Consume the cvParam element (no child nodes)
                        break;
                    case "userParam":
                        // Schema requirements: zero to many instances of this element
                        paramList.Add(ReadUserParam(reader.ReadSubtree()));
                        reader.Read();
                        break;
                    case "binary":
                        // Schema requirements: zero to many instances of this element
                        // Process the ParamList first.
                        foreach (Param param in paramList)
                        {
                            /*
                             * MUST supply a *child* term of MS:1000572 (binary data compression type) only once
                             *   e.g.: MS:1000574 (zlib compression)
                             *   e.g.: MS:1000576 (no compression)
                             * MUST supply a *child* term of MS:1000513 (binary data array) only once
                             *   e.g.: MS:1000514 (m/z array)
                             *   e.g.: MS:1000515 (intensity array)
                             *   e.g.: MS:1000516 (charge array)
                             *   e.g.: MS:1000517 (signal to noise array)
                             *   e.g.: MS:1000595 (time array)
                             *   e.g.: MS:1000617 (wavelength array)
                             *   e.g.: MS:1000786 (non-standard data array)
                             *   e.g.: MS:1000820 (flow rate array)
                             *   e.g.: MS:1000821 (pressure array)
                             *   e.g.: MS:1000822 (temperature array)
                             * MUST supply a *child* term of MS:1000518 (binary data type) only once
                             *   e.g.: MS:1000521 (32-bit float)
                             *   e.g.: MS:1000523 (64-bit float)
                             */
                            switch (param.Accession)
                            {
                                    // MUST supply a *child* term of MS:1000572 (binary data compression type) only once
                                case "MS:1000574":
                                    //   e.g.: MS:1000574 (zlib compression)
                                    compressed = true;
                                    break;
                                case "MS:1000576":
                                    //   e.g.: MS:1000576 (no compression)
                                    compressed = false;
                                    break;
                                    // MUST supply a *child* term of MS:1000513 (binary data array) only once
                                case "MS:1000514":
                                    //   e.g.: MS:1000514 (m/z array)
                                    bda.ArrayType = ArrayType.m_z_array;
                                    break;
                                case "MS:1000515":
                                    //   e.g.: MS:1000515 (intensity array)
                                    bda.ArrayType = ArrayType.intensity_array;
                                    break;
                                case "MS:1000516":
                                    //   e.g.: MS:1000516 (charge array)
                                    bda.ArrayType = ArrayType.charge_array;
                                    break;
                                case "MS:1000517":
                                    //   e.g.: MS:1000517 (signal to noise array)
                                    bda.ArrayType = ArrayType.signal_to_noise_array;
                                    break;
                                case "MS:1000595":
                                    //   e.g.: MS:1000595 (time array)
                                    bda.ArrayType = ArrayType.time_array;
                                    break;
                                case "MS:1000617":
                                    //   e.g.: MS:1000617 (wavelength array)
                                    bda.ArrayType = ArrayType.wavelength_array;
                                    break;
                                case "MS:1000786":
                                    //   e.g.: MS:1000786 (non-standard data array)
                                    bda.ArrayType = ArrayType.non_standard_data_array;
                                    break;
                                case "MS:1000820":
                                    //   e.g.: MS:1000820 (flow rate array)
                                    bda.ArrayType = ArrayType.flow_rate_array;
                                    break;
                                case "MS:1000821":
                                    //   e.g.: MS:1000821 (pressure array)
                                    bda.ArrayType = ArrayType.pressure_array;
                                    break;
                                case "MS:1000822":
                                    //   e.g.: MS:1000822 (temperature array)
                                    bda.ArrayType = ArrayType.temperature_array;
                                    break;
                                    // MUST supply a *child* term of MS:1000518 (binary data type) only once
                                case "MS:1000521":
                                    //   e.g.: MS:1000521 (32-bit float)
                                    bda.Precision = Precision.Precision32;
                                    break;
                                case "MS:1000523":
                                    //   e.g.: MS:1000523 (64-bit float)
                                    bda.Precision = Precision.Precision64;
                                    break;
                            }
                        }
                        byte[] bytes = Convert.FromBase64String(reader.ReadElementContentAsString()); // Consumes the start and end elements.
                        int dataSize = 8;
                        if (bda.Precision == Precision.Precision32)
                        {
                            dataSize = 4;
                        }
                        if (compressed)
                        {
                            bytes = DecompressZLib(bytes, bda.ArrayLength * dataSize);
                        }
                        if (bytes.Length % dataSize != 0 || bytes.Length / dataSize != bda.ArrayLength)
                        {
                            // We need to fail out.
                        }
                        //byte[] oneNumber = new byte[dataSize];
                        //bool swapBytes = true;
                        bda.Data = new double[bda.ArrayLength];
                        for (int i = 0; i < bytes.Length; i += dataSize)
                        {
                            // mzML binary data should always be Little Endian. Some other data formats may use Big Endian, which would require a byte swap
                            //Array.Copy(bytes, i, oneNumber, 0, dataSize);
                            //if (swapBytes)
                            //{
                            //  Array.Reverse(oneNumber);
                            //}
                            if (dataSize == 4)
                            {
                                //bda.Data[i / dataSize] = BitConverter.ToSingle(oneNumber, 0);
                                bda.Data[i / dataSize] = BitConverter.ToSingle(bytes, i);
                            }
                            else if (dataSize == 8)
                            {
                                //bda.Data[i / dataSize] = BitConverter.ToDouble(oneNumber, 0);
                                bda.Data[i / dataSize] = BitConverter.ToDouble(bytes, i);
                            }
                        }
                        break;
                    default:
                        reader.Skip();
                        break;
                }
            }
            reader.Close();
            return bda;
        }

        /*********************************************************************************************************************************************
         * TODO: Flesh out the algorithm/double check it, etc.
         * Do some more work here.
         * 
         ********************************************************************************************************************************************/
        private byte[] DecompressZLib(byte[] compressedBytes, int expectedBytes)
        {
            var msCompressed = new MemoryStream(compressedBytes);
            // We must skip the first two bytes
            // See http://george.chiramattel.com/blog/2007/09/deflatestream-block-length-does-not-match.html
            // EAT the zlib headers, the rest is a normal 'deflate'd stream
            msCompressed.ReadByte();
            msCompressed.ReadByte();
            //var msInflated = new MemoryStream((int)(msCompressed.Length * 2));
            //var newBytes = new byte[msCompressed.Length * 2];
            byte[] newBytes = new byte[expectedBytes];
            // The last 32 bits (4 bytes) are supposed to be an Adler-32 checksum. Might need to remove them as well.
            using (var inflater = new DeflateStream(msCompressed, CompressionMode.Decompress))
            {
                var bytesRead = inflater.Read(newBytes, 0, expectedBytes);
                if (bytesRead != expectedBytes)
                {
                    throw new XmlException("Fail decompressing data...");
                }
                //while (inflater.CanRead)
                //{
                //  var readBytes = new byte[4095];
                //  // Should be able to change to just this.
                //  var bytesRead = inflater.Read(readBytes, 0, readBytes.Length);
                //  if (bytesRead != 0)
                //  {
                //      msInflated.Write(readBytes, 0, bytesRead);
                //  }
                //}
            }
            //newBytes = new byte[msInflated.Length];
            //msInflated.Read(newBytes, 0, (int)msInflated.Length);
            return newBytes;
        }
    }
}