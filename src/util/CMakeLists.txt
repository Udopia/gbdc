add_library(util OBJECT 
    CNFFormula.h
    ResourceLimits.h
    SolverTypes.h
    Stamp.h
    StreamBuffer.h
    UnionFind.cc
    CaptureDistribution.cc
)
set_property(TARGET util PROPERTY POSITION_INDEPENDENT_CODE ON)