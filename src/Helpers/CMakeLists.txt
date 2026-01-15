add_library(helper_lib STATIC Chem.cpp Postprocess.cpp YamlHelpers.cpp Preprocess.cpp)

if (SOAR_BUILD_API)
   target_sources(helper_lib PRIVATE InferenceAPI.cpp SysUtils.cpp)
endif()

target_link_libraries(helper_lib PUBLIC yaml-cpp::yaml-cpp)
target_include_directories(helper_lib PUBLIC ${INCLUDE_DIR})
