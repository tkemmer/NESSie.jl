using NESSie.TestModel
using NESSie.TestModel: scalemodel, legendre, spherical_besseli, spherical_besselk

context("XieModel") do
    for T in testtypes
        charges = [
            Charge(T[0, 0, 0, 1]...),
            Charge(T[1, 0, 0, 2]...),
            Charge(T[1, 1, 0, 3]...),
            Charge(T[0, 1, 0, 4]...)
        ]

        model = XieModel(2one(T), charges)
        @fact typeof(model) --> XieModel{T}
        @fact model.radius --> 2one(T)
        for (q1, q2) in zip(model.charges, scalemodel(charges, model.radius))
            @fact q1.pos --> q2.pos
            @fact q1.val --> q2.val
        end
        @fact model.params --> defaultopt(T)
    end
end

context("scalemodel") do
    for T in testtypes
        charges = [
            Charge(T[0, 0, 0, 1]...),
            Charge(T[1, 0, 0, 2]...),
            Charge(T[1, 1, 0, 3]...),
            Charge(T[0, 1, 0, 4]...)
        ]

        # translation only
        scaled = scalemodel(charges, T(√(.5)/.8))
        @fact typeof(scaled) --> Vector{Charge{T}}
        @fact scaled[1].pos --> roughly(T[-.5, -.5, 0])
        @fact scaled[1].val --> one(T)
        @fact scaled[2].pos --> roughly(T[.5, -.5, 0])
        @fact scaled[2].val --> 2one(T)
        @fact scaled[3].pos --> roughly(T[.5, .5, 0])
        @fact scaled[3].val --> 3one(T)
        @fact scaled[4].pos --> roughly(T[-.5, .5, 0])
        @fact scaled[4].val --> 4one(T)

        # translation and scaling
        scaled = scalemodel(charges, T(1.25))
        @fact typeof(scaled) --> Vector{Charge{T}}
        @fact scaled[1].pos --> roughly(T[-√.5, -√.5, 0])
        @fact scaled[1].val --> one(T)
        @fact scaled[2].pos --> roughly(T[√.5, -√.5, 0])
        @fact scaled[2].val --> 2one(T)
        @fact scaled[3].pos --> roughly(T[√.5, √.5, 0])
        @fact scaled[3].val --> 3one(T)
        @fact scaled[4].pos --> roughly(T[-√.5, √.5, 0])
        @fact scaled[4].val --> 4one(T)
    end
end

context("legendre") do
    for T in testtypes
        # 0-length
        p = legendre(0, one(T))
        @fact isa(p, Function) --> true
        @fact_throws p(-1)
        @fact_throws p(0)
        @fact_throws p(1)

        # one element
        p = legendre(1, one(T))
        @fact isa(p, Function) --> true
        @fact typeof(p(0)) --> T
        @fact p(0) --> one(T)
        @fact_throws p(-1)
        @fact_throws p(1)

        # more elements
        p = legendre(5, one(T)/2)
        @fact isa(p, Function) --> true
        @fact typeof(p(0)) --> T
        @fact p(0) --> one(T)
        @fact typeof(p(1)) --> T
        @fact p(1) --> one(T)/2
        @fact typeof(p(2)) --> T
        @fact p(2) --> roughly(T(-0.125))
        @fact typeof(p(3)) --> T
        @fact p(3) --> roughly(T(-0.4375))
        @fact typeof(p(4)) --> T
        @fact p(4) --> roughly(T(-0.2890625))
        @fact_throws p(-1)
        @fact_throws p(5)
    end
end

context("spherical_besseli") do
    for T in testtypes
        # 0-length
        i = spherical_besseli(-2, one(T))
        @fact isa(i, Function) --> true
        @fact_throws i(-2)
        @fact_throws i(-1)
        @fact_throws i(0)

        # some elements
        i = spherical_besseli(0, one(T)/2)
        @fact isa(i, Function) --> true
        @fact typeof(i(-1)) --> T
        @fact i(-1) --> roughly(T(√π * 1.2723896474148502))
        @fact typeof(i(0)) --> T
        @fact i(0) --> roughly(T(√π * 0.5879930867904171))
        @fact_throws i(-2)
        @fact_throws i(1)
    end
end

context("spherical_besselk") do
    for T in testtypes
        # 0-length
        k = spherical_besselk(-2, one(T))
        @fact isa(k, Function) --> true
        @fact_throws k(-2)
        @fact_throws k(-1)
        @fact_throws k(0)

        # some elements
        k = spherical_besselk(0, one(T)/2)
        @fact isa(k, Function) --> true
        @fact typeof(k(-1)) --> T
        @fact k(-1) --> roughly(T(√π * 1.0750476034999203))
        @fact typeof(k(0)) --> T
        @fact k(0) --> roughly(T(√π * 1.0750476034999203))
        @fact_throws k(-2)
        @fact_throws k(1)
    end
end

@pending NonlocalXieModel1 --> :nothing
@pending coefficients --> :nothing
@pending φΣ --> :nothing
