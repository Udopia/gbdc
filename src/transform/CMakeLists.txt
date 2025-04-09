add_library(transform OBJECT 
    cnf2bip.cc
    cnf2cnf.h
    cnf2kis.h
)
set_property(TARGET transform PROPERTY POSITION_INDEPENDENT_CODE ON)
