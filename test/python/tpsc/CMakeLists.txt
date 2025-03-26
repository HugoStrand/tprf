FILE(GLOB all_tgz_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.tar.gz)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${all_tgz_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set(PREFIX tpsc-)

add_python_test(test_tpsc_chi_sum_rule ${PREFIX})
add_python_test(test_tpsc_G_Sigma_sum_rule ${PREFIX})
add_python_test(test_tpsc_hubbard_square_lattice ${PREFIX})
add_python_test(test_tpsc_hubbard_triangular_lattice ${PREFIX})
add_python_test(test_tpsc_improved_bubble ${PREFIX})
add_python_test(test_tpsc_plus_Sigma ${PREFIX})
add_python_test(test_tpsc_solver_h5 ${PREFIX})
add_python_test(test_triqs_TBLattice ${PREFIX})