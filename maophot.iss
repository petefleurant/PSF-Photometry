; Script generated by the Inno Setup Script Wizard.
; SEE THE DOCUMENTATION FOR DETAILS ON CREATING INNO SETUP SCRIPT FILES!

#define MyAppName "MAOPhot"
#define MyAppVersion "1.1.0"
#define MyAppPublisher "Pete Fleurant"
#define MyAppURL "https://github.com/petefleurant/PSF-Photometry"
#define MyAppExeName "MAOPhot.exe"
#define MyAppAssocName MyAppName + " File"
#define MyAppAssocExt ".exe"
#define MyAppAssocKey StringChange(MyAppAssocName, " ", "") + MyAppAssocExt

[Setup]
; NOTE: The value of AppId uniquely identifies this application. Do not use the same AppId value in installers for other applications.
; (To generate a new GUID, click Tools | Generate GUID inside the IDE.)
AppId={{8983D6B5-1CDF-45A7-B19C-4BD911C107B2}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
;AppVerName={#MyAppName} {#MyAppVersion}
AppPublisher={#MyAppPublisher}
AppPublisherURL={#MyAppURL}
AppSupportURL={#MyAppURL}
AppUpdatesURL={#MyAppURL}
DefaultDirName={localappdata}\{#MyAppName}
; "ArchitecturesAllowed=x64compatible" specifies that Setup cannot run
; on anything but x64 and Windows 11 on Arm.
ArchitecturesAllowed=x64compatible
; "ArchitecturesInstallIn64BitMode=x64compatible" requests that the
; install be done in "64-bit mode" on x64 or Windows 11 on Arm,
; meaning it should use the native 64-bit Program Files directory and
; the 64-bit view of the registry.
ArchitecturesInstallIn64BitMode=x64compatible
ChangesAssociations=yes
DisableProgramGroupPage=yes
LicenseFile=E:\Astronomy\AAVSO\PSF-Photometry\PSF-Photometry\LICENSE
; Uncomment the following line to run in non administrative install mode (install for current user only.)
;PrivilegesRequired=lowest
OutputBaseFilename=MAOPhot_SETUP
Compression=lzma
SolidCompression=yes
WizardStyle=modern

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"; Flags: unchecked

[Dirs]
Name: "{app}\logs"; Permissions: users-modify
Name: "{app}\2022 8 1 V1117 Her (Example)"; Permissions: users-modify
Name: "{app}\2022 12 21 Z Tau (Example)"; Permissions: users-modify

[Files]
Source: "E:\Astronomy\AAVSO\PSF-Photometry\PSF-Photometry\dist\{#MyAppExeName}"; DestDir: "{app}"; Flags: ignoreversion
Source: "E:\Astronomy\AAVSO\PSF-Photometry\PSF-Photometry\2022 8 1 V1117 Her (Example)\*"; DestDir: "{app}\2022 8 1 V1117 Her (Example)"; Flags: ignoreversion recursesubdirs createallsubdirs; Permissions: users-modify
Source: "E:\Astronomy\AAVSO\PSF-Photometry\PSF-Photometry\2022 12 21 Z Tau (Example)\*"; DestDir: "{app}\2022 12 21 Z Tau (Example)"; Flags: ignoreversion recursesubdirs createallsubdirs; Permissions: users-modify
Source: "E:\Astronomy\AAVSO\PSF-Photometry\PSF-Photometry\HelpFile.pdf"; DestDir: "{app}"; Flags: ignoreversion
Source: "E:\Astronomy\AAVSO\PSF-Photometry\PSF-Photometry\LICENSE"; DestDir: "{app}"; Flags: ignoreversion
Source: "E:\Astronomy\AAVSO\PSF-Photometry\PSF-Photometry\README.md"; DestDir: "{app}"; Flags: ignoreversion
; NOTE: Don't use "Flags: ignoreversion" on any shared system files

[Registry]
Root: HKA; Subkey: "Software\Classes\{#MyAppAssocExt}\OpenWithProgids"; ValueType: string; ValueName: "{#MyAppAssocKey}"; ValueData: ""; Flags: uninsdeletevalue
Root: HKA; Subkey: "Software\Classes\{#MyAppAssocKey}"; ValueType: string; ValueName: ""; ValueData: "{#MyAppAssocName}"; Flags: uninsdeletekey
Root: HKA; Subkey: "Software\Classes\{#MyAppAssocKey}\DefaultIcon"; ValueType: string; ValueName: ""; ValueData: "{app}\{#MyAppExeName},0"
Root: HKA; Subkey: "Software\Classes\{#MyAppAssocKey}\shell\open\command"; ValueType: string; ValueName: ""; ValueData: """{app}\{#MyAppExeName}"" ""%1"""
Root: HKA; Subkey: "Software\Classes\Applications\{#MyAppExeName}\SupportedTypes"; ValueType: string; ValueName: ".myp"; ValueData: ""

[Icons]
Name: "{autoprograms}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}"
Name: "{autodesktop}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}"; Tasks: desktopicon

[Run]
Filename: "{app}\{#MyAppExeName}"; Description: "{cm:LaunchProgram,{#StringChange(MyAppName, '&', '&&')}}"; Flags: nowait postinstall skipifsilent

