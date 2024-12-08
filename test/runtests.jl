using TestItemRunner
using TestItems

# For CI, run all tests not marked with the :extended tag
@run_package_tests filter=ti -> !(:extended in ti.tags) verbose=true
