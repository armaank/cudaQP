get_directory_property(headers
                        DIRECTORY ${PROJECT_SOURCE_DIR}/svm
                        DEFINITION headers)

set(headers ${headers}
${CMAKE_CURRENT_SOURCE_DIR}/test_solve_linsys.h PARENT_SCOPE)

get_directory_property(codegen_headers
                        DIRECTORY ${PROJECT_SOURCE_DIR}/svm
                        DEFINITION codegen_headers)

set(codegen_headers ${codegen_headers}
        ${CMAKE_CURRENT_SOURCE_DIR}/data.h PARENT_SCOPE)
