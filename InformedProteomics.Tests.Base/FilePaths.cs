using System.IO;

namespace InformedProteomics.Tests.Base
{
    public static class FilePaths
    {
        public static readonly string TestRawFilePath = Path.Combine(Utils.DEFAULT_TEST_FILE_FOLDER, @"TopDown\ProductionQCShew\QC_Shew_12_02_2_1Aug12_Cougar_12-06-11.raw");
        public static readonly string TestPbfFilePath = Path.Combine(Utils.DEFAULT_TEST_FILE_FOLDER, @"TopDown\ProductionQCShew\QC_Shew_12_02_2_1Aug12_Cougar_12-06-11.pbf");
        public static readonly string TestTopDownRawFilePathEtd = Path.Combine(Utils.DEFAULT_TEST_FILE_FOLDER, @"TopDown\E_coli_iscU_60_mock.raw");
        public static readonly string TestTopDownRawFilePathCid = Path.Combine(Utils.DEFAULT_TEST_FILE_FOLDER, @"TopDown\SBEP_STM_001_02272012_Aragon.raw");
        public static readonly string TestQExactiveRawFilePath = Path.Combine(Utils.DEFAULT_TEST_FILE_FOLDER, @"TopDown\ProductionQCShew\QC_Shew_13_04_A_17Feb14_Samwise_13-07-28.raw");
        public static readonly string TestQExactivePbfFilePath = Path.Combine(Utils.DEFAULT_TEST_FILE_FOLDER, @"TopDown\ProductionQCShew\QC_Shew_13_04_A_17Feb14_Samwise_13-07-28.pbf");
    }
}