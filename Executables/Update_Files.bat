if not exist PbfGen mkdir PbfGen
if not exist PbfGen mkdir ProMex
if not exist PbfGen mkdir MSPathFinder

xcopy ..\PbfGen\bin\Release\*.exe                                   PbfGen\ /D /Y
xcopy ..\PbfGen\bin\Release\*.dll                                   PbfGen\ /D /Y
xcopy ..\PbfGen\bin\Release\*.pdb                                   PbfGen\ /D /Y
xcopy ..\PbfGen\bin\Release\PbfGen.exe.config                       PbfGen\ /D /Y 
xcopy ..\PbfGen\bin\Release\ThermoRawFileReader.xml                 PbfGen\ /D /Y 
xcopy ..\PbfGen\bin\Release\InformedProteomics.Backend.xml          PbfGen\ /D /Y 

xcopy ..\ProMex\bin\Release\*.exe                                   ProMex\ /D /Y
xcopy ..\ProMex\bin\Release\*.dll                                   ProMex\ /D /Y
xcopy ..\ProMex\bin\Release\*.pdb                                   ProMex\ /D /Y
xcopy ..\ProMex\bin\Release\InformedProteomics.Backend.xml          ProMex\ /D /Y  
xcopy ..\ProMex\bin\Release\ProMex.exe.config                       ProMex\ /D /Y  
xcopy ..\ProMex\bin\Release\ThermoRawFileReader.xml                 ProMex\ /D /Y  

xcopy ..\MSPathFinderT\bin\Release\*.exe                            MSPathFinder\ /D /Y
xcopy ..\MSPathFinderT\bin\Release\*.dll                            MSPathFinder\ /D /Y
xcopy ..\MSPathFinderT\bin\Release\*.pdb                            MSPathFinder\ /D /Y
xcopy ..\MSPathFinderT\bin\Release\MSPathFinderT.exe.config         MSPathFinder\ /D /Y
xcopy ..\MSPathFinderT\bin\Release\ThermoRawFileReader.xml          MSPathFinder\ /D /Y
xcopy ..\MSPathFinderT\bin\Release\InformedProteomics.Backend.xml   MSPathFinder\ /D /Y
xcopy ..\MSPathFinderT\bin\Release\InformedProteomics.Backend.Database.xml MSPathFinder\ /D /Y

pause
