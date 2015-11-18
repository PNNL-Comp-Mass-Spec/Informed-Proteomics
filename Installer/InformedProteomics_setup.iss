; This is an Inno Setup configuration file
; http://www.jrsoftware.org/isinfo.php

#define ApplicationVersion GetFileVersion('..\InformedProteomics.Backend\bin\x64\Release\InformedProteomics.Backend.dll')

[CustomMessages]
AppName=DeconTools
[Messages]
WelcomeLabel2=This will install [name/ver] on your computer.%n%n%n%nNOTICE:%nSome source files require access to a 64-bit ProteoWizard installation. Please install 64-bit ProteoWizard before using the program to avoid errors.%n%n
; Example with multiple lines:
; WelcomeLabel2=Welcome message%n%nAdditional sentence
[Files]
; InformedProteomics.Backend
Source: InformedProteomics.Backend\bin\x64\Release\InformedProteomics.Backend.dll                                                                 ; DestDir: {app}
Source: InformedProteomics.Backend\bin\x64\Release\InformedProteomics.Backend.dll.config                                                          ; DestDir: {app}

; InformedProteomics.Backend Nuget libraries
Source: InformedProteomics.Backend\bin\x64\Release\MathNet.Numerics.dll                                                                           ; DestDir: {app}

; Manually managed libraries
Source: lib\alglibnet2.dll                                                                                                                        ; DestDir: {app}
Source: lib\PNNLOmics.dll                                                                                                                         ; DestDir: {app}
Source: lib\ProteinFileReader.dll                                                                                                                 ; DestDir: {app}
Source: lib\x64\SAIS.dll                                                                                                                          ; DestDir: {app}
Source: lib\ThermoRawFileReaderDLL.dll                                                                                                            ; DestDir: {app}

; InformedProteomics.BottomUp
Source: InformedProteomics.BottomUp\bin\x64\Release\InformedProteomics.BottomUp.dll                                                               ; DestDir: {app}

; InformedProteomics.Graphics
Source: InformedProteomics.Graphics\bin\x64\Release\InformedProteomics.Graphics.dll                                                               ; DestDir: {app}

; InformedProteomics.Graphics Nuget libraries
Source: InformedProteomics.Graphics\bin\x64\Release\OxyPlot.dll                                                                                   ; DestDir: {app}
Source: InformedProteomics.Graphics\bin\x64\Release\OxyPlot.Wpf.dll                                                                               ; DestDir: {app}

; InformedProteomics.Scoring
Source: InformedProteomics.Scoring\bin\x64\Release\InformedProteomics.Scoring.dll                                                                 ; DestDir: {app}

; InformedProteomics.TopDown
Source: InformedProteomics.TopDown\bin\x64\Release\InformedProteomics.TopDown.dll                                                                 ; DestDir: {app}

; InformedProteomics.TopDown Nuget libraries
;Source: InformedProteomics.TopDown\bin\x64\Release\MathNet.Numerics.dll                                                                                   ; DestDir: {app}

; MSPathFinder (bottom up)
;Source: MSPathFinder\bin\x64\Release\MSPathFinder.exe                                                                                             ; DestDir: {app}
;Source: MSPathFinder\bin\x64\Release\MSPathFinder.exe.config                                                                                      ; DestDir: {app}

; MSPathFinder (top down)
Source: TopDownConsole\bin\x64\Release\MSPathFinderT.exe                                                                                           ; DestDir: {app}
Source: TopDownConsole\bin\x64\Release\MSPathFinderT.exe.config                                                                                    ; DestDir: {app}

; PbfGen
Source: PbfGen\bin\x64\Release\PbfGen.exe                                                                                                         ; DestDir: {app}
Source: PbfGen\bin\x64\Release\PbfGen.exe.config                                                                                                  ; DestDir: {app}

; ProMex
Source: MS1FeatureExtractor\bin\x64\Release\ProMex.exe                                                                                                         ; DestDir: {app}
Source: MS1FeatureExtractor\bin\x64\Release\ProMex.exe.config                                                                                                  ; DestDir: {app}

Source: README.md                                                                                                                                 ; DestDir: {app}\Readme.txt
;Source: RevisionHistory.txt                                                                                                                       ; DestDir: {app}
Source: Example_Files\Mod_Examples.txt                                                                                                            ; DestDir: {app}
Source: Example_Files\MSPathFinder_Mods.txt                                                                                                       ; DestDir: {app}
Source: Example_Files\MSPathFinder_Mods_Phospho.txt                                                                                               ; DestDir: {app}


[Dirs]
Name: {commonappdata}\InformedProteomics; Flags: uninsalwaysuninstall

;[Tasks]
;Name: desktopicon; Description: {cm:CreateDesktopIcon}; GroupDescription: {cm:AdditionalIcons}; Flags: unchecked
; Name: quicklaunchicon; Description: {cm:CreateQuickLaunchIcon}; GroupDescription: {cm:AdditionalIcons}; Flags: unchecked

[Icons]
Name: {group}\Command Line; Filename: {cmd}; Parameters: "/K ""set PATH=%PATH%;{app}"""; WorkingDir: {userdocs}; Comment: Command Prompt with direct access to InformedProteomics executables
Name: {group}\ReadMe; Filename: {app}\Readme.txt; Comment: InformedProteomics ReadMe
Name: {group}\Mod Examples; Filename: {app}\Mod_Examples.txt; Comment: Example Modifications
Name: {group}\MSPathFinder Mods; Filename: {app}\MSPathFinder_Mods.txt; Comment: Example Modifications file for MSPathFinder
Name: {group}\MSPathFinder Mods Phospho; Filename: {app}\MSPathFinder_Mods_Phospho.txt; Comment: Example Modifications file for MSPathFinder
Name: {group}\Uninstall InformedProteomics; Filename: {uninstallexe}

;C:\Windows\System32\cmd.exe /K "set PATH=%PATH%;%ProgramFiles%\InformedProteomics"

[Setup]
ArchitecturesInstallIn64BitMode=x64
ArchitecturesAllowed=x64
AppName=InformedProteomics
AppVersion={#ApplicationVersion}
;AppVerName=InformedProteomics
AppID=InformedProteomicsId
AppPublisher=Pacific Northwest National Laboratory
AppPublisherURL=http://omics.pnl.gov/software
AppSupportURL=http://omics.pnl.gov/software
AppUpdatesURL=http://omics.pnl.gov/software
DefaultDirName={pf}\InformedProteomics
DefaultGroupName=InformedProteomics
AppCopyright=© PNNL
;LicenseFile=.\License.rtf
PrivilegesRequired=poweruser
OutputBaseFilename=InformedProteomics_Installer
;VersionInfoVersion=1.57
VersionInfoVersion={#ApplicationVersion}
VersionInfoCompany=PNNL
VersionInfoDescription=InformedProteomics
VersionInfoCopyright=PNNL
DisableFinishedPage=true
ShowLanguageDialog=no
;SetupIconFile=..\MageFileProcessor\wand.ico
;InfoBeforeFile=.\readme.rtf
ChangesAssociations=false
;WizardImageFile=..\Deploy\Images\MageSetupSideImage.bmp
;WizardSmallImageFile=..\Deploy\Images\MageSetupSmallImage.bmp
;InfoAfterFile=.\postinstall.rtf
EnableDirDoesntExistWarning=false
AlwaysShowDirOnReadyPage=true
UninstallDisplayIcon={app}\delete_16x.ico
ShowTasksTreeLines=true
OutputDir=Installer\Output
SourceDir=..\
Compression=lzma
SolidCompression=yes
[Registry]
;Root: HKCR; Subkey: MageFile; ValueType: string; ValueName: ; ValueData:Mage File; Flags: uninsdeletekey
;Root: HKCR; Subkey: MageSetting\DefaultIcon; ValueType: string; ValueData: {app}\wand.ico,0; Flags: uninsdeletevalue
[UninstallDelete]
Name: {app}; Type: filesandordirs
