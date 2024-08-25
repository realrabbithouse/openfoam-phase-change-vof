    scalar alphaCoNum = -GREAT;

    forAll(fluidRegions, regioni)
    {
        alphaCoNum = max
        (
            alphaCourantNo
            (
                fluidRegions[regioni],
                runTime,
                mixtureFluid[regioni],
                phiFluid[regioni]
            ),
            alphaCoNum
        );
    }
