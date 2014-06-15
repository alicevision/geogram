:: Get the Visual Studio installation dir
set MSVCDIR=%VS110COMNTOOLS%..\..\VC

:: Configure environment for Visual Studio
call "%MSVCDIR%\BIN\VCVARS32.BAT"

:: Set the generator to use
set CMAKE_VS_GENERATOR=Visual Studio 11
set CMAKE_VS_GENERATOR_TOOLSET=v110_xp

