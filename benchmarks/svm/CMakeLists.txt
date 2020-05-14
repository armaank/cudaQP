get_directory_property(headers
                        DIRECTORY ${PROJECT_SOURCE_DIR}/benchmarks
                        DEFINITION headers)

set(headers ${headers}
${CMAKE_CURRENT_SOURCE_DIR}/test_svm.h PARENT_SCOPE)

get_directory_property(codegen_headers
                        DIRECTORY ${PROJECT_SOURCE_DIR}/benchmarks
                        DEFINITION codegen_headers)

set(codegen_headers ${codegen_headers}
        ${CMAKE_CURRENT_SOURCE_DIR}/data.h PARENT_SCOPE)
