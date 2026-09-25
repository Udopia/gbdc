add_library(transform OBJECT 
    cnf2bip.cc
    cnf2cnf.cc
    cnf2kis.h
    cnf2lig.cc
)
set_property(TARGET transform PROPERTY POSITION_INDEPENDENT_CODE ON)
