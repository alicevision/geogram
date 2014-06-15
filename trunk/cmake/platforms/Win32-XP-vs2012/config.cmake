set(VORPALINE_ARCH_32 true)
include(${CMAKE_SOURCE_DIR}/cmake/platforms/Windows-vs.cmake)

# Target Windows XP SP2 operating system (and later)
# see http://msdn.microsoft.com/fr-fr/library/aa383745.aspx
add_definitions(-D_WIN32_WINNT=0x0502)
