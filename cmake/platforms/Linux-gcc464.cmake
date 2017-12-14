
# This configuration is for a compiler configured as follows:
# 
# cd gcc-4.6.4
# mkdir build
# cd build
# ../configure \
#        --prefix=/opt/gcc464 \
#        --program-suffix=464 \
#        --enable-languages=c,c++ \
#        --enable-shared --enable-threads=posix --disable-checking \
#        --with-system-zlib --enable-__cxa_atexit --disable-libunwind-exceptions \
#        --disable-multilib

set(CMAKE_CXX_COMPILER:FILEPATH /opt/gcc464/bin/g++464)
set(CMAKE_C_COMPILER:FILEPATH /opt/gcc464/bin/gcc464)

include(${GEOGRAM_SOURCE_DIR}/cmake/platforms/Linux-gcc.cmake)

