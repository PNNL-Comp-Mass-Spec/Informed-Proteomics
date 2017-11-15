if not exist PbfGen mkdir PbfGen
if not exist PbfGen mkdir ProMex
if not exist PbfGen mkdir MSPathFinder

xcopy ..\PbfGen\bin\Release\*.exe PbfGen\ /D /Y
xcopy ..\PbfGen\bin\Release\*.dll PbfGen\ /D /Y
xcopy ..\PbfGen\bin\Release\*.pdb PbfGen\ /D /Y

xcopy ..\ProMex\bin\Release\*.exe ProMex\ /D /Y
xcopy ..\ProMex\bin\Release\*.dll ProMex\ /D /Y
xcopy ..\ProMex\bin\Release\*.pdb ProMex\ /D /Y

xcopy ..\MSPathFinderT\bin\Release\*.exe MSPathFinder\ /D /Y
xcopy ..\MSPathFinderT\bin\Release\*.dll MSPathFinder\ /D /Y
xcopy ..\MSPathFinderT\bin\Release\*.pdb MSPathFinder\ /D /Y

pause
