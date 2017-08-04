--             pixel.lua
-- FR bibliotheque interne, incluse quand on appelle require("pixel")
-- FR fournit des fonctions simples pour afficher des (gros) pixels
-- EN internal library, included when one calls require("pixel")
-- EN defines easy-to-use functions for displaying (big) pixels

function GLUP.Pixel3D(x,y,z)
   GLUP.Vertex(x,y,z)
   GLUP.Vertex(x+1,y,z)
   GLUP.Vertex(x,y+1,z)
   GLUP.Vertex(x+1,y+1,z)      
   GLUP.Vertex(x,y,z+1)
   GLUP.Vertex(x+1,y,z+1)
   GLUP.Vertex(x,y+1,z+1)
   GLUP.Vertex(x+1,y+1,z+1) 
end

function GLUP.Pixel2D(x,y)
   GLUP.Pixel3D(x,y,0.0)
end

pix = GLUP.Pixel2D
pix3d = GLUP.Pixel3D
col = GLUP.Color

function pixBegin()
   GLUP.Enable(GLUP.VERTEX_COLORS)
   GLUP.Begin(GLUP.HEXAHEDRA)
end

function pixEnd()
   GLUP.End()
   GLUP.Disable(GLUP.VERTEX_COLORS)
end

function pixGrid()
   GLUP.SetColor(GLUP.FRONT_COLOR,0,0,0.5)
   GLUP.Begin(GLUP.LINES)
   local xm,ym,zm,xM,yM,zM
   xm,ym,zm,xM,yM,zM = GLUP.GetRegionOfInterest()
   for x=xm,xM,1 do
      GLUP.Vertex(x,ym,0)
      GLUP.Vertex(x,yM,0)
   end
   for y=ym,yM,1 do
      GLUP.Vertex(xm,y,0)
      GLUP.Vertex(xM,y,0)
   end
   GLUP.End()
end

function GLUP.init_graphics()
   GLUP.SetRegionOfInterest(1,1,1,21,21,1)
end
