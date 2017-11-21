cmake_minimum_required(VERSION 2.8)

set(LIB_SOURCES
    Bed2SetPair.cpp
    Bed2SetPair.h
    BedCompare.cpp
    BedCompare.h
    Bedpe.cpp
    Bedpe.h
    CMakeLists.txt
    CutterByRule.cpp
    CutterByRule.h
    ExecuteByRule.cpp
    ExecuteByRule.h
    Gene2.cpp
    Gene2.h
    MyTypes.h
    Printer.cpp
    Printer.h
    SetPair.cpp
    SetPair.h
    Util.cpp
    Util.h
    main.cpp
)

add_library(fusionEvaluator ${LIB_SOURCES})
add_executable(fusionEvaluator-bin main.cpp)
target_link_libraries(fusionEvaluator-bin fusionEvaluator)
set_target_properties(fusionEvaluator-bin PROPERTIES OUTPUT_NAME fusionToolEvaluator)
