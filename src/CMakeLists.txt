add_subdirectory(Helpers)

add_executable(train train.cpp)
add_executable(infer inference.cpp)
add_executable(data_prep data_preparation.cpp)

target_link_libraries(train PRIVATE yaml-cpp::yaml-cpp CNum::CNum helper_lib)
target_link_libraries(infer PRIVATE CNum::CNum helper_lib)
target_link_libraries(data_prep PRIVATE CNum::CNum helper_lib)

if (SOAR_BUILD_API)
   add_executable(api_driver api_driver.cpp)
   target_link_libraries(api_driver PRIVATE CNum::CNum CNum::CNum_Deploy helper_lib)
endif()

if (SOAR_BUILD_FAST)
	target_compile_options(train PRIVATE -Ofast -fno-math-errno -fno-trapping-math -ffp-contract=fast -ffast-math)
else()
	target_compile_options(train PRIVATE -O2)
endif()

target_compile_options(data_prep PRIVATE -O2)
