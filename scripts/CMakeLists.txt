
# BUILD files in scripts
add_executable(coo_csr_gen_code coo_csr_gen_code.cc)
add_executable(bcsr_csr_copy bcsr_csr_copy.cc)
add_executable(csr_coo_gen_code synth_csr_coo.cc)

set_target_properties(coo_csr_gen_code 
	PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/scripts"
    COMPILE_FLAGS "-g"
)

set_target_properties(csr_coo_gen_code 
	PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/scripts"
    COMPILE_FLAGS "-g"
)
include_directories(./)
add_custom_target(scripts
	DEPENDS coo_csr_gen_code
                csr_coo_gen_code
	COMMENT "Building scripts from synthesis"
	)
