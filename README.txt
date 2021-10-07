
# Files to be produced

#--------------------------------------------

1) Files with radiation/smearing/acceptance cuts: Options 3                                                                            # default (option 4 / G2018 doesn't work)

2) Files without radiation but with smearing/acceptance cuts: Options 1 + 2                                                            # effect of radiation

3) Files with 1p0pi truth, with radiation/smearing/acceptance cuts: Options 3 + 4                                                      # closure test & acceptance correction

#--------------------------------------------

4) Files at truth 1p0pi level without radiation/smearing/acceptance cuts: Options 1 + 2                                                # absolute truth / data-mc comparisons

#--------------------------------------------

5) Files at truth 1p0pi level without radiation, with offset for G2018, with smearing, without acceptance cuts: Options 1 + 2          # acceptance correction

#--------------------------------------------


###########################################

root -b
.L OverlayLFGProxy.cpp
OverlayLFGProxy()
OverlayLFGProxy("","MC")
OverlayLFGProxy("xB_")
OverlayLFGProxy("xB_","MC")

root -b 
.L OverlayPL.cpp
OverlayPL()
OverlayPL("xB_")
OverlayPL("","MC")
OverlayPL("xB_","MC")

root -b
.L PlotTwoD.cpp
PlotTwoD()
PlotTwoD("xB_")
PlotTwoD("","MC")
PlotTwoD("xB_","MC")

root -b
.L OverlayGoodBad.cpp
OverlayGoodBad()
OverlayGoodBad("xB_")
OverlayGoodBad("","MC")
OverlayGoodBad("xB_","MC")

root -b 
.L OverlayPMiss2DSlices.cpp
OverlayPMiss2DSlices()
OverlayPMiss2DSlices("xB_")
OverlayPMiss2DSlices("","MC")
OverlayPMiss2DSlices("xB_","MC")

root -b
.L Overlay.cpp
Overlay()
Overlay("xB_")

###########################################
