; KIT-167 (section 8.6, pool S6): the Microsoft Visual C++ x64 runtime ("Visual C++ v14") the bundled DIAMOND needs on
; Windows. The per-user installer does not elevate, so it does not install the runtime: when the runtime is absent or
; older than the build gate F1 recorded, it says so in words and offers Microsoft's download page, and the installation
; continues. Without the runtime the application refuses a new mapping in words (reference_mapping_engine_runtime_missing);
; saved mappings and the demo tables still work. The registry is read in the 64-bit view. No runtime DLL is shipped
; (FSF GPL FAQ, "WindowsRuntimeAndGPL").
;
; MSR_VC_MIN_BLD: the oldest runtime build with which gate F1's search ran (run 45353725258, Server Core ltsc2025):
; without a runtime and with the oldest one Microsoft serves (14.0.24215) the program stops with 0xC0000135 (DLL not
; found: it imports MSVCP140_1.dll and VCRUNTIME140_1.dll, which that runtime lacks); with the current one (build 35211)
; it runs. No build between the two was tried, so the check asks for 35211 or later.
!ifndef MSR_VC_MIN_BLD
  !define MSR_VC_MIN_BLD 35211
!endif

!macro customInit
  IfSilent msr_vc_done
  SetRegView 64
  ClearErrors
  ReadRegDWORD $R8 HKLM "SOFTWARE\Microsoft\VisualStudio\14.0\VC\Runtimes\x64" "Installed"
  ReadRegDWORD $R9 HKLM "SOFTWARE\Microsoft\VisualStudio\14.0\VC\Runtimes\x64" "Bld"
  SetRegView 32
  StrCmp $R8 "1" 0 msr_vc_missing
  IntCmp $R9 ${MSR_VC_MIN_BLD} msr_vc_done msr_vc_missing msr_vc_done
  msr_vc_missing:
    MessageBox MB_YESNO|MB_ICONINFORMATION "MultiScholaR maps strain proteins to a reference organism with the program DIAMOND, which needs the Microsoft Visual C++ runtime (Visual C++ v14, 64-bit). It is missing or too old on this computer. The installation continues; everything else works without it.$\r$\n$\r$\nOpen Microsoft's download page for the runtime now?" IDNO msr_vc_done
    ExecShell "open" "https://learn.microsoft.com/cpp/windows/latest-supported-vc-redist"
  msr_vc_done:
!macroend
