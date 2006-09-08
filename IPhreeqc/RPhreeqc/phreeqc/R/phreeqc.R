
# Global Package Variables

.packageNameDLL = paste(.packageName, .Platform$dynlib.ext, sep="")

# Package Initialization

.onLoad =
function(libname, pkgname)
{
    library.dynam(.packageNameDLL, pkgname, libname)

    # If we need to initialize structures do it in a phreeqc_init C function.
    # Might need to finalize them with a _finalize function in .onUnload too.
    #invisible(.Call("phreeqc_init", PACKAGE=.packageName))
}

# Package Functions

phrReadDB =
function(filename)
{
    invisible(.Call("read_db", as.character(filename), PACKAGE=.packageName))
}

phrReadString =
function(string)
{
    invisible(.Call("read", as.character(string), PACKAGE=.packageName))
}

phrRun =
function()
{
    invisible(.Call("run", PACKAGE=.packageName))
}

phrRunFile =
function(filename)
{
    invisible(.Call("runFile", as.character(filename), PACKAGE=.packageName))
}

#phrGetColumn =
#function(column_number)
#{
#    return(.Call("getCol", as.integer(column_number), PACKAGE=.packageName))
#}

phrGetSelectedOutput =
function()
{
    return(.Call("getSelOut", PACKAGE=.packageName))
}

phrGetLastErrorString =
function()
{
    return(.Call("getLastErrorString", PACKAGE=.packageName))
}
