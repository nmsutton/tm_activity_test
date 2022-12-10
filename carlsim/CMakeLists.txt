project(moradi_stp_compare)

# Targets  
add_executable(moradi_stp_compare
    src/main_moradi_stp_compare.cpp
)

# definitions     
if(CARLSIM_NO_CUDA)
    target_compile_definitions(moradi_stp_compare PUBLIC
        -D__NO_CUDA__ )
endif()

# Linking    
target_link_libraries(moradi_stp_compare PRIVATE
    carlsim
    carlsim-stopwatch
)
target_include_directories(moradi_stp_compare PRIVATE
    ../../tools/stopwatch
)
   
        
# Installation   
install(TARGETS moradi_stp_compare
    EXPORT carlsim
    RUNTIME 
        DESTINATION $<$<CONFIG:Debug>:debug/>samples 
        COMPONENT moradi_stp_compare
        )

#install(FILES run_moradi_stp_compare DESTINATION $<$<CONFIG:Debug>:debug/>samples)          