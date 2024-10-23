add_library(transform OBJECT 
    IndependentSet.h
    Normalize.h
    cnf2bip.h
)
set_property(TARGET transform PROPERTY POSITION_INDEPENDENT_CODE ON)
