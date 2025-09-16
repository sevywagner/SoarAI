add_subdirectory(Helpers)

add_executable(train train.cpp)
add_executable(infer inference.cpp)
add_executable(data_prep data_preparation.cpp)
target_link_libraries(train PRIVATE CNum::CNum helper_lib)
target_link_libraries(infer PRIVATE CNum::CNum helper_lib)
target_link_libraries(data_prep PRIVATE CNum::CNum helper_lib)
