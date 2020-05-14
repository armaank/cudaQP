# Add subdirectories of all the tests
add_subdirectory(svm)

set(test_headers ${headers} PARENT_SCOPE)
set(codegen_test_headers ${codegen_headers} PARENT_SCOPE)
