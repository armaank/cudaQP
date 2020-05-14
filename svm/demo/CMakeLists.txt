get_directory_property(headers
                        DIRECTORY ${PROJECT_SOURCE_DIR}/svm
                        DEFINITION headers)
                        
set(headers ${headers}
    ${CMAKE_CURRENT_SOURCE_DIR}/demo.h PARENT_SCOPE)

