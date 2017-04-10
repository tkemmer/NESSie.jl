#=
    Utility function to
    - Create a temporary file,
    - Invoke a given function that presumably writes to the same file,
    - Clean up the temporary file, and
    - Return the file's full content.

    @param f
        Function f(fh::IOStream) that fills the file represented by fh
    @return String
=#
function readback(f::Function)
    (fname, fh) = mktemp()
    try
        f(fh)
        seek(fh, 0)
        return join(readlines(fh), "\n")
    finally
        close(fh)
        rm(fname)
    end
end
