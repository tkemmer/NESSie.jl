using NESSie.TestModel

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
        for (q1, q2) in zip(model.charges, TestModel.scalemodel(charges, model.radius))
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
        scaled = TestModel.scalemodel(charges, T(√(.5)/.8))
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
        scaled = TestModel.scalemodel(charges, T(1.25))
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
