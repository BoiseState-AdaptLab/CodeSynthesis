
# BUILD files in scripts
add_executable(new_coo_csr_gen_code new_coo_csr_gen_code.cc)
add_executable(bcsr_csr_copy bcsr_csr_copy.cc)
add_executable(csr_coo_gen_code synth_csr_coo.cc)
add_executable(non_substitute_coo_csr non_substitute_coo_csr.cc)
add_executable(coo_dia_trial coo_dia_trial_synthesis.cc)
add_executable(bench_harness bench_harness.cc)
add_executable(make_fake_mtx_files make_fake_mtx_files.cpp)

set_target_properties(csr_coo_gen_code 
	PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/scripts"
    COMPILE_FLAGS "-g"
)
set_target_properties(make_fake_mtx_files
		PROPERTIES
		RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/scripts"
		COMPILE_FLAGS "-g -Wall -Werror -Wextra"
		)

#set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
#set (CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")

set_target_properties(bench_harness
		PROPERTIES
		RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/scripts"
		COMPILE_FLAGS "-g -Wall -Wextra"
		)

set_target_properties(coo_dia_trial
		PROPERTIES
		RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/scripts"
		COMPILE_FLAGS "-g"
		)

set_target_properties(new_coo_csr_gen_code 
	PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/scripts"
    COMPILE_FLAGS "-g"
)


set_target_properties(non_substitute_coo_csr 
	PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/scripts"
    COMPILE_FLAGS "-g"
)


include_directories(./)
add_custom_target(scripts
	DEPENDS bench_harness 
                new_coo_csr_gen_code
		csr_coo_gen_code
		coo_dia_trial
		non_substitute_coo_csr
	COMMENT "Building scripts from synthesis"
	)
