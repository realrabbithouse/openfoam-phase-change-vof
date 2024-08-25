    scalar CoNum = -GREAT;

    forAll(fluidRegions, regioni)
    {
        CoNum = max
        (
            fluidCourantNo
            (
                fluidRegions[regioni],
                runTime,
                phiFluid[regioni]
            ),
            CoNum
        );
    }
