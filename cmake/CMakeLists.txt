include (ExternalProject)

set_property (DIRECTORY PROPERTY EP_BASE Dependencies)

set (DEPENDENCIES)
set (EXTRA_CMAKE_ARGS)

set( boost_download "https://boostorg.jfrog.io/artifactory/main/release/1.83.0/source/boost_1_83_0.tar.gz")
set( boost_SHA256 "c0685b68dd44cc46574cce86c4e17c0f611b15e195be9848dfd0769a0a207628" )
set( boost_INSTALL ${CMAKE_CURRENT_BINARY_DIR}/third_party/boost )
set( boost_INCLUDE_DIR ${boost_INSTALL}/include )
set( boost_LIB_DIR ${boost_INSTALL}/lib )

unset(Boost_DIR CACHE)
unset(Boost_INCLUDE_DIR CACHE)
unset(Boost_LIBRARY_DIR_DEBUG CACHE)
unset(Boost_LIBRARY_DIR_RELEASE CACHE)
unset(Boost_PROGRAM_OPTIONS_LIBRARY_DEBUG CACHE)


# Use static linking to avoid issues with system-wide installations of Boost.
list (APPEND DEPENDENCIES ep_boost)
ExternalProject_Add (ep_boost
  URL ${boost_download}
  URL_HASH SHA256=${boost_SHA256}
  CONFIGURE_COMMAND ./bootstrap.sh --with-libraries=program_options
  BUILD_COMMAND ./b2 link=static
  BUILD_IN_SOURCE 1
  INSTALL_COMMAND ""
  )
list (APPEND EXTRA_CMAKE_ARGS
  -DBOOST_ROOT=${CMAKE_CURRENT_BINARY_DIR}/Dependencies/Source/ep_boost
  -DBoost_NO_SYSTEM_PATHS=ON)

ExternalProject_Add (TimeEvolver
  DEPENDS ${DEPENDENCIES}
  SOURCE_DIR ${PROJECT_SOURCE_DIR}
  CMAKE_ARGS -DDOWNLOAD_BOOST=OFF -DNOSYSTEMBOOST=ON ${EXTRA_CMAKE_ARGS}
  INSTALL_COMMAND ""
  BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})
