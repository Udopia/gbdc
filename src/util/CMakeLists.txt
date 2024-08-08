add_library(util OBJECT 
    CNFFormula.h
    ResourceLimits.h
    SolverTypes.h
    Stamp.h
    StreamBuffer.h
)

add_subdirectory(threadpool)
add_subdirectory(pybind11)