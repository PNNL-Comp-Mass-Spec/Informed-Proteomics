set FilterClause=--where "cat != PNL_Domain && cat != Local_Testing"
set NUnitExe="C:\Program Files (x86)\NUnit.org\nunit-console\nunit3-console.exe"

%NUnitExe% ..\InformedProteomics.Tests.UnitTests\bin\Debug\InformedProteomics.Tests.UnitTests.dll %FilterClause% --out=NUnitConsole_UnitTests.txt --result:NUnit_UnitTests.xml --timeout=90000
%NUnitExe% ..\InformedProteomics.Tests.FunctionalTests\bin\Debug\InformedProteomics.Tests.FunctionalTests.dll %FilterClause% --out=NUnitConsole_FunctionalTests.txt --result:NUnit_FunctionalTests.xml --timeout=90000
%NUnitExe% ..\InformedProteomics.Tests.DevTests\bin\Debug\InformedProteomics.Tests.DevTests.dll %FilterClause% --out=NUnitConsole_DevTests.txt --result:NUnit_DevTests.xml --timeout=90000
%NUnitExe% ..\InformedProteomics.Test\bin\Debug\InformedProteomics.Test.dll %FilterClause% --out=NUnitConsole_Test.txt --result:NUnit_Test.xml --timeout=90000

pause
