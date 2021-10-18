
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
