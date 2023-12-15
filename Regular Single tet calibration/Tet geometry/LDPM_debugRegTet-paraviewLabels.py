
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get active source.
        
lDPMgeo000parasingleTetFacets000vtk = FindSource('LDPM_debugRegTet-para-facets.000.vtk')
                             
SetActiveSource(lDPMgeo000parasingleTetFacets000vtk)

# create a query selection
QuerySelect(QueryString='(id >= 0)', FieldType='CELL', InsideOut=0)

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# get display properties
lDPMgeo000parasingleTetFacets000vtkDisplay = GetDisplayProperties(lDPMgeo000parasingleTetFacets000vtk, view=renderView1)

# Properties modified on lDPMgeo000parasingleTetFacets000vtkDisplay
lDPMgeo000parasingleTetFacets000vtkDisplay.SelectionCellLabelVisibility = 1

# Properties modified on lDPMgeo000parasingleTetFacets000vtkDisplay
lDPMgeo000parasingleTetFacets000vtkDisplay.SelectionCellLabelColor = [0.0, 0.0, 0.0]
lDPMgeo000parasingleTetFacets000vtkDisplay.SelectionCellLabelFontSize = 14
lDPMgeo000parasingleTetFacets000vtkDisplay.SelectionPointLabelFontSize = 14


#--------------------------------------------
# uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
        
