add_library(util OBJECT 
    CNFFormula.h
    ResourceLimits.h
    SolverTypes.h
    Stamp.h
    StreamBuffer.h
    UnionFind.cc
    CaptureDistribution.cc
)

add_subdirectory(pybind11)