if(${BUILD_TESTS})
        set(link_tests "driver;auxiliary;core;sweepers;solvers;${HDF5_LIBRARIES};pugixml")

        add_unit_test(test_3x3 ${link_tests})
        copy_file_if_changed(${CMAKE_CURRENT_SOURCE_DIR}/3x3.xml
            ${CMAKE_CURRENT_BINARY_DIR}/3x3.xml test_3x3)
        copy_file_if_changed(${CMAKE_CURRENT_SOURCE_DIR}/c5g7.xsl
            ${CMAKE_CURRENT_BINARY_DIR}/c5g7.xsl test_3x3)

        add_unit_test(test_moc_sn_3d ${link_tests})
        copy_file_if_changed(${CMAKE_CURRENT_SOURCE_DIR}/3d_moc.xml
            ${CMAKE_CURRENT_BINARY_DIR}/3d_moc.xml test_moc_sn_3d)
        copy_file_if_changed(${CMAKE_CURRENT_SOURCE_DIR}/3d_sn.xml
            ${CMAKE_CURRENT_BINARY_DIR}/3d_sn.xml test_moc_sn_3d)
        copy_file_if_changed(${CMAKE_CURRENT_SOURCE_DIR}/c5g7.xsl
            ${CMAKE_CURRENT_BINARY_DIR}/c5g7.xsl test_moc_sn_3d)
endif()
