project(func_test LANGUAGES CXX)

#TV -> FD : J'ai effectué des modifications afin de pouvoir compiler les tests sur ma version d'ubuntu 18.04
#set(Boost_USE_STATIC_LIBS ON)
#TODO : a supprimer parce que plus Boost (a vérifier)
if (WIN32)
  set(Boost_USE_STATIC_LIBS ON)
else()
  add_compile_definitions(BOOST_TEST_DYN_LINK)
endif()

find_package(Catch2 CONFIG REQUIRED)
include(Catch)

# include_directories(test_file)

aux_source_directory(. SOURCES_LIST)
foreach(filename ${SOURCES_LIST})
  get_filename_component(target ${filename} NAME_WE)
  add_executable(${target} ${filename})
  target_link_libraries(${target} libsource Catch2::Catch2)
  target_compile_features(${target} PRIVATE cxx_std_17)
  # detail each sub test
  catch_discover_tests(${target} WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/func_test)
endforeach(filename)