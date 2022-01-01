
# BUILD files in scripts
add_executable(coo_csr_gen_code coo_csr_gen_code.cc)
add_executable(bcsr_csr_copy bcsr_csr_copy.cc)

set_target_properties(coo_csr_gen_code 
	PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/scripts"
    COMPILE_FLAGS "-g"
)

add_custom_target(scripts
	DEPENDS coo_csr_gen_code
		bcsr_csr_copy
	COMMENT "Building scripts from synthesis"
	)
