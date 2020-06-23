# RStar Tree Library
set(rstar_sources
        RStarTree.c
        RSTInstDel.c
        RSTJoin.c
        RSTUtil.c
        RSTInOut.c
        RSTInterUtil.c
        RSTQuery.c)

set(rstar_headers
        RStarTree.h
        RSTInstDel.h
        RSTJoin.h
        RSTUtil.h
        RSTInOut.h
        RSTInterUtil.h
        RSTQuery.h
        RSTBase.h)


add_library(cSmRST)
target_sources(cSmRST
        PRIVATE ${rstar_sources}
        PUBLIC ${rstar_headers})
target_compile_definitions(cSmRST
        PUBLIC FILE_OFFSET_BITS=64)