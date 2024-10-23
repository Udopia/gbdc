add_library(extract OBJECT 
    CNFBaseFeatures.cc
    CNFGateFeatures.cc
    OPBBaseFeatures.cc
    WCNFBaseFeatures.cc
)
set_property(TARGET extract PROPERTY POSITION_INDEPENDENT_CODE ON)