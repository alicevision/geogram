#snapshot.sh: generates an image of a given mesh under different viewing angles,
#using vorpaview and imagemagick.

#******************************************************************************
#usage: snapshot.sh meshfile
#generated image is in output.png
#******************************************************************************

# Size of each snapshot
geometry=800x800 

# Views to generate. For each view: x,y,z,angle
# where x,y,z is the axis of rotation.
views="1,1,1,70 1,0,0,90 0,1,0,90 0,0,1,0"

# Simulated keypress events (zooms a little bit and hides mesh)
keys=zzzzzzzzm

for view in $views
do
    image=/tmp/snap_`echo $view|sed -e 's|,|_|g'`.ppm
    images="$images $image"
    vorpaview $1 gfx:geometry=$geometry gfx:background=false \
	      gfx:rotate=$view gfx:snapshot=$image \
	      gfx:keypress=$keys > /dev/null
#    convert $image -trim /tmp/tmp.ppm
#    mv /tmp/tmp.ppm $image
done

# Assemble the images using imagemagick
montage $images -adjoin output.png
