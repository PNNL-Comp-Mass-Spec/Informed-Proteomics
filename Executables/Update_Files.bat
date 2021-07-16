@echo off
echo.
echo Be sure to build PbfGen, ProMex, and MSPathFinderT in Release Mode
echo.
pause

if not exist PbfGen mkdir PbfGen
if not exist PbfGen mkdir ProMex
if not exist PbfGen mkdir ProMexAlign
if not exist PbfGen mkdir MSPathFinder

echo.
echo Copying PbgGen files
xcopy ..\PbfGen\bin\Release\*.exe                                   PbfGen\ /D /Y
xcopy ..\PbfGen\bin\Release\*.dll                                   PbfGen\ /D /Y
xcopy ..\PbfGen\bin\Release\*.pdb                                   PbfGen\ /D /Y
xcopy ..\PbfGen\bin\Release\PbfGen.exe.config                       PbfGen\ /D /Y 
xcopy ..\PbfGen\bin\Release\InformedProteomics.Backend.xml          PbfGen\ /D /Y 

echo.
echo Copying ProMex files
xcopy ..\ProMex\bin\Release\*.exe                                   ProMex\ /D /Y
xcopy ..\ProMex\bin\Release\*.dll                                   ProMex\ /D /Y
xcopy ..\ProMex\bin\Release\*.pdb                                   ProMex\ /D /Y
xcopy ..\ProMex\bin\Release\InformedProteomics.Backend.xml          ProMex\ /D /Y  
xcopy ..\ProMex\bin\Release\ProMex.exe.config                       ProMex\ /D /Y  

echo.
echo Copying ProMexAlign files
xcopy ..\PromexAlign\bin\Release\*.exe                              ProMexAlign\ /D /Y
xcopy ..\PromexAlign\bin\Release\*.dll                              ProMexAlign\ /D /Y
xcopy ..\PromexAlign\bin\Release\*.pdb                              ProMexAlign\ /D /Y
xcopy ..\PromexAlign\bin\Release\InformedProteomics.Backend.xml     ProMexAlign\ /D /Y  
xcopy ..\PromexAlign\bin\Release\ProMexAlign.exe.config             ProMexAlign\ /D /Y  

echo.
echo Copying MSPathFinderT files
xcopy ..\MSPathFinderT\bin\Release\*.exe                            MSPathFinder\ /D /Y
xcopy ..\MSPathFinderT\bin\Release\*.dll                            MSPathFinder\ /D /Y
xcopy ..\MSPathFinderT\bin\Release\*.pdb                            MSPathFinder\ /D /Y
xcopy ..\MSPathFinderT\bin\Release\MSPathFinderT.exe.config         MSPathFinder\ /D /Y
xcopy ..\MSPathFinderT\bin\Release\InformedProteomics.Backend.xml   MSPathFinder\ /D /Y
xcopy ..\MSPathFinderT\bin\Release\InformedProteomics.Backend.Database.xml MSPathFinder\ /D /Y

pause
