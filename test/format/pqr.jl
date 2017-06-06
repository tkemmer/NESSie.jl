using NESSie.Format

testfiles = ((mktemp()..., 0), # empty file
             (mktemp()..., 2)) # dummy file
write(testfiles[2][2], """
REMARK   1 Lorem ipsum dolor sit amet.
REMARK   1
ATOM      1  N   THR     1     -17.108  25.866  23.850  0.1812 1.8240
HETATM    2  O   HOH     2     -21.160  40.444  40.509 -0.8340 1.6612
ATOM      3  CA  THR     3     -16.775  27.193  23.310  0.0034 1.9080
TER
END""")

context("readpqr") do
    try
        for (fname, fh, len) in testfiles, T in testtypes
            seekstart(fh)
            charges = readpqr(fh, T)

            @fact typeof(charges) --> Vector{Charge{T}}
            @fact length(charges) --> len
            if fname == testfiles[2][1]
                @fact charges[1].pos --> T[-17.108, 25.866, 23.850]
                @fact charges[1].val --> T(0.1812)
                @fact charges[2].pos --> T[-16.775, 27.193, 23.310]
                @fact charges[2].val --> T(0.0034)
            end
        end
    finally
        [begin close(fh); rm(fn); end for (fn, fh, _) in testfiles]
    end
end
